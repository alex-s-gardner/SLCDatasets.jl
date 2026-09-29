# The merged image itself: a matrix whose samples are read from the bursts when they are asked for.
#
# A subswath is a few gigabytes and a correlator reads it in windows, so nothing is copied up front. The
# array holds a [`BurstGrid`](@ref) and one raster per burst, and a read resolves each row to the burst
# covering it. Rows and samples outside every burst's valid region read as zero — the merged image is
# wider and taller than the union of the regions placed in it, and the corners it leaves are not data.
#
# Reading a window is what this is built for: `A[rows, cols]` groups the rows by the burst they come
# from and takes one strip from each, where scalar indexing pays a row lookup per sample. That is the
# operation `AutoRIFT.jl`'s blocked path uses.
#
# Zero is what the margin reads as because the correlator treats it as absent, and because a filter over
# a window straddling the edge must not see whatever the processor left in the file. It is not a claim
# that no sample there is zero: [`validmask`](@ref) is what says where the data is, and a caller that
# needs the distinction should use it rather than testing for zero.

"""
    BurstRaster(raster, row_offset = 0)

One burst's samples: the raster holding them, and where in it the burst begins.

A `.SAFE` stores a subswath's bursts stacked in one raster, so each burst begins at its own offset into
it; a burst delivered as its own file begins at the start. Pairing the two makes those the same thing to
a reader.
"""
struct BurstRaster{S<:AbstractMatrix}
    raster::S
    # Where in `raster` the burst begins: zero when it holds that burst alone, and
    # `(burst - 1) * lines_per_burst` when it holds the whole subswath stacked, which is how a `.SAFE`
    # measurement raster stores it. A placement's `burst_rows` count from the burst's own first line, so
    # this is what turns them into rows of the file.
    row_offset::Int
end

BurstRaster(raster::AbstractMatrix) = BurstRaster(raster, 0)

# The rows of the file a placement's rows come from.
_file_rows(s::BurstRaster, rows::AbstractUnitRange) = (first(rows) + s.row_offset):(last(rows) + s.row_offset)

"""
    ConcatenatedBursts{T,S} <: DiskArrays.AbstractDiskArray{T,2}

The bursts of one Sentinel-1 subswath as a single matrix, read on demand.

Rows come from whichever burst covers them under the array's [`BurstGrid`](@ref); rows and samples no
burst reaches read as `zero(T)`. Build one with [`pixels`](@ref) on a merged [`SLC`](@ref) rather than
calling this constructor.

Index with ranges — `A[rows, cols]` — wherever the shape of the read allows it: a window is taken as
one strip per burst it spans, while scalar indexing resolves a row per sample.

A `DiskArrays` array, so views, broadcasting and strided indexing behave as they do for every other raster
this package returns. `eachchunk` reports the rows each burst serves, which `burst_grid` makes
disjoint and consecutive.
"""
struct ConcatenatedBursts{T,S<:AbstractMatrix{T}} <: DiskArrays.AbstractDiskArray{T,2}
    # Parallel to `grid.placements`: `sources[k]` is what `grid.placements[k]` reads from. One entry per
    # placement rather than one per product, so bursts delivered as separate files and bursts sharing a
    # subswath's raster are the same array.
    sources::Vector{BurstRaster{S}}
    grid::BurstGrid
end

function ConcatenatedBursts(sources::AbstractVector{<:BurstRaster}, grid::BurstGrid)
    length(sources) == length(grid.placements) || throw(ArgumentError(
        "the grid places $(length(grid.placements)) bursts but $(length(sources)) rasters were " *
        "given; there must be one raster per placed burst"))
    for (p, src) in zip(grid.placements, sources)
        rows = _file_rows(src, p.burst_rows)
        # A raster too small for the rows and samples its placement takes would read outside it.
        (last(rows) <= size(src.raster, 1) && last(p.cols) <= size(src.raster, 2)) ||
            throw(ArgumentError(
                "burst $(p.burst) contributes rows $rows and samples $(p.cols) of a raster that is " *
                "$(size(src.raster, 1))x$(size(src.raster, 2))"))
    end
    held = collect(sources)
    return ConcatenatedBursts{eltype(eltype(held)),fieldtype(eltype(held), :raster)}(held, grid)
end

# One raster shared by every placement, each reading at its burst's own offset into it.
ConcatenatedBursts(raster::AbstractMatrix, grid::BurstGrid, row_offsets) =
    ConcatenatedBursts([BurstRaster(raster, o) for o in row_offsets], grid)

# One raster per placement, each holding that burst alone.
ConcatenatedBursts(rasters::AbstractVector{<:AbstractMatrix}, grid::BurstGrid) =
    ConcatenatedBursts(map(BurstRaster, rasters), grid)

Base.size(A::ConcatenatedBursts) = size(A.grid)
Base.IndexStyle(::Type{<:ConcatenatedBursts}) = IndexCartesian()

"""
    grid(A::ConcatenatedBursts) -> BurstGrid

Where each burst sits in the merged image.
"""
grid(A::ConcatenatedBursts) = A.grid

function Base.getindex(A::ConcatenatedBursts{T}, i::Int, j::Int) where {T}
    @boundscheck checkbounds(A, i, j)
    k = placement_index(A.grid, i)
    k == 0 && return zero(T)
    p = A.grid.placements[k]
    j in p.cols || return zero(T)
    src = A.sources[k]
    return src.raster[src.row_offset + first(p.burst_rows) + (i - first(p.grid_rows)), j]
end

DiskArrays.readblock!(A::ConcatenatedBursts, dest::AbstractArray,
                      rows::AbstractUnitRange{<:Integer},
                      cols::AbstractUnitRange{<:Integer}) = _read_window!(dest, A, rows, cols)

"""
    DiskArrays.eachchunk(A::ConcatenatedBursts)

The rows each burst serves, full width, with any rows no burst reaches as bands of their own.

`burst_grid` trims each placement to its valid lines, so the placements are disjoint and consecutive — on a
real subswath, `19:1421`, `1422:2764`, `2765:4106` and so on. That makes them a partition of the azimuth
axis, which is what a chunk grid has to be, and the granularity a row actually belongs to.

The leading and trailing rows no burst images — the first eighteen of that subswath, which read as zero —
become bands of their own, since a grid that skipped them would leave rows in no chunk at all.

A grid whose placements are *not* disjoint and ascending cannot be expressed this way, and reports one chunk
covering the array instead. `burst_grid` does not build one, but the constructor accepts a hand-made grid.
"""
function DiskArrays.eachchunk(A::ConcatenatedBursts)
    bands = UnitRange{Int}[]
    next = 1
    for p in A.grid.placements
        # Out of order or overlapping: not a partition, so there is no honest grid to report.
        first(p.grid_rows) < next && return DiskArrays.GridChunks(A, size(A))
        first(p.grid_rows) > next && push!(bands, next:(first(p.grid_rows) - 1))
        push!(bands, p.grid_rows)
        next = last(p.grid_rows) + 1
    end
    next <= size(A, 1) && push!(bands, next:size(A, 1))
    isempty(bands) && return DiskArrays.GridChunks(A, size(A))
    return DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = map(length, bands)),
                                 DiskArrays.RegularChunks(size(A, 2), 0, size(A, 2)))
end

# Reported rather than hidden behind `Unchunked`: a band is the real unit a row belongs to. Nothing forces one
# to be read whole — the sources are uncompressed strips in a mapped file — and declaring it costs nothing,
# measured at 2.3 ns a sample for a 512² window against 2.2 under `Unchunked`, because `DiskArrays`' batch
# strategy takes a dense sub-range directly rather than expanding it to chunk bounds.
DiskArrays.haschunks(::ConcatenatedBursts) = DiskArrays.Chunked()

# One strip per burst the window spans, taken from that burst's raster in a single read. The window is
# zeroed first, so the rows and samples no burst covers need no separate pass.
function _read_window!(out, A::ConcatenatedBursts{T}, rows, cols) where {T}
    fill!(out, zero(T))
    each_overlap(A.grid, rows, cols) do k, p, grid_rows, window_cols
        src = A.sources[k]
        shift = src.row_offset + first(p.burst_rows) - first(p.grid_rows)
        strip = src.raster[(first(grid_rows) + shift):(last(grid_rows) + shift), window_cols]
        copyto!(_window_view(out, rows, cols, grid_rows, window_cols), strip)
    end
    return out
end

# The part of a window-shaped result that a placement's rows and columns land in. Offsets come from `out`'s
# own axes, since `DiskArrays` hands `readblock!` whatever array the caller is filling and for a view that is
# not one-based.
function _window_view(out, rows, cols, grid_rows, window_cols)
    i0, j0 = first(axes(out, 1)), first(axes(out, 2))
    return view(out,
                (i0 + first(grid_rows) - first(rows)):(i0 + last(grid_rows) - first(rows)),
                (j0 + first(window_cols) - first(cols)):(j0 + last(window_cols) - first(cols)))
end

"""
    BurstValidMask <: AbstractMatrix{Bool}

Which samples of a merged image were imaged, derived from the grid rather than stored.

`true` exactly where [`ConcatenatedBursts`](@ref) reads from a burst's valid region. The image reads
zero outside, but a sample inside can be zero too, so this is what distinguishes absent from dark and
what a correlator should be given rather than a test against zero.
"""
struct BurstValidMask <: AbstractMatrix{Bool}
    grid::BurstGrid
end

Base.size(m::BurstValidMask) = size(m.grid)
Base.IndexStyle(::Type{BurstValidMask}) = IndexCartesian()

function Base.getindex(m::BurstValidMask, i::Int, j::Int)
    @boundscheck checkbounds(m, i, j)
    p = placement_at(m.grid, i)
    return p !== nothing && j in p.cols
end

function Base.getindex(m::BurstValidMask, rows::AbstractUnitRange{<:Integer},
                       cols::AbstractUnitRange{<:Integer})
    @boundscheck checkbounds(m, rows, cols)
    out = fill(false, length(rows), length(cols))
    each_overlap(m.grid, rows, cols) do _, _, grid_rows, window_cols
        fill!(_window_view(out, rows, cols, grid_rows, window_cols), true)
    end
    return out
end

"""
    validmask(A) -> AbstractMatrix{Bool}

Which samples of `A` were imaged.

For a merged subswath this is [`BurstValidMask`](@ref), computed from the burst placements and holding
nothing, so asking for it costs no read of the rasters.
"""
validmask(A::ConcatenatedBursts) = BurstValidMask(A.grid)
