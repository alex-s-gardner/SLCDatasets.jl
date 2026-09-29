# The NISAR backend: HDF5 products under `/science/{LSAR,SSAR}/<productType>`.
#
# Group paths are built the way isce3's `nisar.products.readers.Base` builds them — a root path found
# by probing for the sensor band, then the product type, then a fixed subgroup name. The product type
# is read from the identification group and selects the reader, which is also how isce3's
# `open_product` dispatches. A product whose type is `SLC` is an early-mission RSLC.

const SCIENCE_ROOT = "science"
const SENSOR_BANDS = ("LSAR", "SSAR")

# Speed of light, for converting the processed center frequency to a wavelength. This is the value
# isce3 uses (`isce3.core.speed_of_light`), and the conversion must match it to the bit.
const SPEED_OF_LIGHT = 299792458.0

"""
    NisarBackend <: AbstractSLCBackend

A NISAR-format HDF5 product: a file, the sensor band group inside it, and the product type.

`frequency` selects the sub-band (`"A"` or `"B"`); a product lists the ones it carries at
`identification/listOfFrequencies`.
"""
struct NisarBackend <: AbstractSLCBackend
    path::String
    band::String
    product_type::String
    frequency::String
end

product_path(b::NisarBackend) = string(SCIENCE_ROOT, "/", b.band, "/", b.product_type)
metadata_path(b::NisarBackend) = string(product_path(b), "/metadata")
swath_path(b::NisarBackend) = string(product_path(b), "/swaths")
frequency_path(b::NisarBackend) = string(swath_path(b), "/frequency", b.frequency)
identification_path(b::NisarBackend) = string(SCIENCE_ROOT, "/", b.band, "/identification")
orbit_path(b::NisarBackend) = string(metadata_path(b), "/orbit")

# A geocoded product stores its samples on a map grid under `grids` rather than `swaths`, so it has no
# slant-range axis for `RadarGeometry` to describe and is out of this package's scope. Recognized by name
# in order to say that, rather than failing on a missing group.
const GEOCODED_TYPES = ("GSLC", "GCOV", "GUNW", "GOFF")

"""
    nisar_band(h) -> String

The sensor band group of an open NISAR product: `"LSAR"` or `"SSAR"`.
"""
function nisar_band(h)::String
    haskey(h, SCIENCE_ROOT) || throw(ArgumentError(
        "`$(HDF5.filename(h))` has no `/$SCIENCE_ROOT` group, so it is not a NISAR-format product"))
    science = _group(h, SCIENCE_ROOT)
    for band in SENSOR_BANDS
        haskey(science, band) && return band
    end
    throw(ArgumentError("`$(HDF5.filename(h))` has none of $(join(SENSOR_BANDS, ", ")) under " *
                        "`/$SCIENCE_ROOT`, so it is not a NISAR-format product"))
end

nisar_band(path::AbstractString) = h5open(nisar_band, path, "r")

"""
    nisar_product_type(h, band) -> String

The product type recorded in an open NISAR product's identification group, as the group name that
holds it. Early-mission products name the group `SLC` where later ones name it `RSLC`; both report
`RSLC`.
"""
function nisar_product_type(h, band::AbstractString)::String
    g = _group(h, string(SCIENCE_ROOT, "/", band))
    haskey(g, "identification") || throw(ArgumentError(
        "`$(HDF5.filename(h))` has no `identification` group under `/$SCIENCE_ROOT/$band`"))
    declared = _text(g, "identification/productType")
    # The declared type names the group except for the early-mission spelling.
    declared == "SLC" && return haskey(g, "RSLC") ? "RSLC" : "SLC"
    return declared
end

nisar_product_type(path::AbstractString, band::AbstractString) =
    h5open(h -> nisar_product_type(h, band), path, "r")

# Indexing an HDF5 file or group gives back `Union{Attribute,Dataset,Datatype,Group}`, so `read` on the
# result infers as `Any` and every value derived from it dispatches at runtime. These accessors assert
# what the layout says a name is, which both keeps inference concrete and turns a product whose
# structure differs from the NISAR layout into an error naming the path rather than a `MethodError`
# deeper in.

@noinline _not_a(kind, parent, name) = throw(ArgumentError(
    "`$name` under `$(HDF5.name(parent))` is a $(nameof(typeof(parent))) member that is not a " *
    "$kind, so this is not a NISAR-format product"))

function _dataset(parent, name::AbstractString)
    d = parent[name]
    d isa HDF5.Dataset || _not_a("dataset", parent, name)
    return d
end

function _group(parent, name::AbstractString)
    g = parent[name]
    g isa HDF5.Group || _not_a("group", parent, name)
    return g
end

# `read` on a dataset infers as `Any`, because a dataset's element type is only known once the file is
# open. So each accessor below converts to the type the layout calls for and annotates that, which keeps
# inference concrete through the arithmetic that follows. The conversion is explicit rather than
# `read(d, T)`, which reinterprets the stored bytes as `T` and errors on a width mismatch instead of
# converting — a product storing a count as `Int32` or times as `Float32` is read correctly here.

# A dataset's value where only a scalar is wanted.
_scalar(parent, name::AbstractString)::Float64 = _only_number(_dataset(parent, name))

# A whole-number field, stored variously as an integer or a float across product generations.
function _scalar_int(parent, name::AbstractString)
    d = _dataset(parent, name)
    v = _only_number(d)
    isinteger(v) || throw(ArgumentError(
        "`$(HDF5.name(d))` holds $v where a whole number is expected, so this is not a " *
        "NISAR-format product"))
    return Int(v)
end

# A scalar dataset holds a number, a zero-dimensional array, or a one-element one depending on how the
# product was written; `only` covers the arrays and the bare number falls through.
_only_number(d)::Float64 = _as_number(read(d))
_as_number(x::Number) = Float64(x)
_as_number(x::AbstractArray) = Float64(only(x))

# A fixed-length string dataset, surfaced as a `String` or as bytes; `_string` narrows either and strips
# the NUL padding.
_text(parent, name::AbstractString)::String = _string(read(_dataset(parent, name)))

# The first and last element of a one-dimensional dataset, and its length. Reading the dataset whole
# would transfer the entire axis — 27,000 azimuth times for a real RSLC — where only the ends are used.
function _axis_bounds(parent, name::AbstractString)
    d = _dataset(parent, name)
    n = length(d)
    n > 0 || throw(ArgumentError("`$(HDF5.name(d))` is empty, so it has no bounds"))
    # Only the endpoints are read: the axis itself is tens of thousands of samples.
    return _as_number(d[1]), _as_number(d[n]), n
end

# `read_attribute` has no type-asserting form, so its `Any` is narrowed here rather than at each call.
function _units_epoch(parent, name::AbstractString)
    d = _dataset(parent, name)
    haskey(HDF5.attributes(d), "units") || throw(ArgumentError(
        "`$(HDF5.name(d))` has no `units` attribute, so the epoch its times are measured from is " *
        "unknown"))
    units = read_attribute(d, "units")
    return parse_cf_epoch(_string(units))
end

# Every read in the eager path goes through one open file handle. Opening a NISAR product costs about
# as much as the metadata reads themselves, so a reader that reopened per query would spend most of its
# time in `h5open`.
function read_identification(b::NisarBackend)
    return h5open(h -> read_identification(b, h), b.path, "r")
end

function read_identification(b::NisarBackend, h)
    g = _group(h, identification_path(b))
    # An absent optional field reads as empty rather than failing: a product missing one is still
    # usable, where one missing `absoluteOrbitNumber` is not.
    get_str(name)::String = haskey(g, name) ? _text(g, name) : ""
    return Identification(
        get_str("missionId"),
        get_str("productType"),
        _scalar_int(g, "absoluteOrbitNumber"),
        get_str("orbitPassDirection"),
        get_str("lookDirection"),
        get_str("zeroDopplerStartTime"),
        get_str("zeroDopplerEndTime"),
        get_str("boundingPolygon"),
    )
end

read_geometry(b::NisarBackend) = h5open(h -> read_geometry(b, h), b.path, "r")

function read_geometry(b::NisarBackend, h)
    b.product_type in GEOCODED_TYPES && throw(ArgumentError(
        "`$(b.path)` is a $(b.product_type) product, which stores its samples on a map grid and so " *
        "carries no slant-range/azimuth geometry. This package reads SLCs in radar geometry; open an " *
        "RSLC, or use a raster library for a geocoded product"))
    freq = _group(h, frequency_path(b))
    swaths = _group(h, swath_path(b))
    ident = _group(h, identification_path(b))

    near_range, far_range, nsamples = _axis_bounds(freq, "slantRange")
    sensing_start, sensing_stop, nlines = _axis_bounds(swaths, "zeroDopplerTime")

    look = _text(ident, "lookDirection")
    side = lowercase(look) == "left" ? LookLeft :
           lowercase(look) == "right" ? LookRight :
           throw(ArgumentError(
               "`$(b.path)` records lookDirection \"$look\"; expected \"Left\" or \"Right\""))

    return RadarGeometry(
        near_range,
        far_range,
        _scalar(freq, "slantRangeSpacing"),
        SPEED_OF_LIGHT / _scalar(freq, "processedCenterFrequency"),
        1 / _scalar(swaths, "zeroDopplerTimeSpacing"),
        sensing_start,
        sensing_stop,
        nlines,
        nsamples,
        side,
        # The azimuth times and the orbit's state vector times share one epoch in this format, so the
        # epoch is read once here and reported for both.
        _units_epoch(swaths, "zeroDopplerTime"),
    )
end

read_orbit(b::NisarBackend) = h5open(h -> read_orbit(b, h), b.path, "r")

function read_orbit(b::NisarBackend, h)
    g = _group(h, orbit_path(b))
    time = _read_vector(g, "time")
    # Stored (N, 3) row-major, which HDF5.jl presents as (3, N).
    pos = _read_matrix(g, "position")
    vel = _read_matrix(g, "velocity")
    size(pos, 1) == 3 || throw(ArgumentError(
        "orbit position in `$(b.path)` has leading dimension $(size(pos, 1)), expected 3"))
    axes(pos) == axes(vel) || throw(DimensionMismatch(
        "orbit position and velocity in `$(b.path)` have axes $(axes(pos)) and $(axes(vel))"))
    size(pos, 2) == length(time) || throw(DimensionMismatch(
        "orbit in `$(b.path)` has $(length(time)) times but $(size(pos, 2)) state vectors"))
    get_str(name, default)::String = haskey(g, name) ? _text(g, name) : default
    return StateVectors(
        time,
        _columns_as_svectors(pos),
        _columns_as_svectors(vel),
        _units_epoch(g, "time"),
        get_str("interpMethod", "Hermite"),
        get_str("orbitType", "Custom"),
    )
end

# A product may store these as `Float32`, so the element type is converted rather than propagated into
# `StateVectors`, whose fields are `Float64`. `convert` is a no-op when the file already holds `Float64`.
_read_vector(parent, name::AbstractString)::Vector{Float64} =
    convert(Vector{Float64}, read(_dataset(parent, name))::AbstractVector)

_read_matrix(parent, name::AbstractString)::Matrix{Float64} =
    convert(Matrix{Float64}, read(_dataset(parent, name))::AbstractMatrix)

# The (3, N) array HDF5.jl hands back, as N 3-vectors.
function _columns_as_svectors(a::AbstractMatrix{Float64})
    Base.require_one_based_indexing(a)
    return [SVector{3,Float64}(a[1, i], a[2, i], a[3, i]) for i in axes(a, 2)]
end

# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------

"""
    NisarRaster{T} <: AbstractMatrix{T}

A NISAR product's sample array, read a window at a time out of its HDF5 dataset.

**The file's layout is the transpose of Julia's.** A NISAR image is stored row-major with azimuth (or
northing) slowest, so `HDF5.jl` reports it as `(range, azimuth)`; this reverses that, so indexing is
`(line, sample)` for a radar-geometry product and `(northing, easting)` for a geocoded one — the order
every other array in this package uses.

Nothing is held: an RSLC band is 50511 x 57760 `ComplexF32` and a GSLC band 161280 x 80640, which are
23 GB and 104 GB of samples against files of 11 GiB and 10 GiB. Index with ranges, so that a window is
one hyperslab read rather than one read per element.
"""
struct NisarRaster{T} <: AbstractMatrix{T}
    path::String
    dataset::String
    dims::Tuple{Int,Int}
end

# **HDF5 is not thread-safe unless the library was built for it, and the shipped one is not.** Two tasks
# calling into the C library segfault rather than racing visibly, so every call a caller can reach from more
# than one task goes through this.
#
# It is therefore the ceiling on how much of a machine a blocked reader can use, and as little as possible
# belongs inside it. `H5Dread` does not qualify: on a filtered dataset it inflates the chunks itself, which
# is the bulk of a window read — hence [`_read_by_chunks!`](@ref), which keeps only the compressed byte read
# under the lock.
const HDF5_IO = ReentrantLock()

# **A window read opens the file, and keeping the handle or the dataset open instead does not work.** A
# `Dict` of live `HDF5.Dataset`s segfaults a blocked reader inside `H5CX_pop` under `H5Tget_member_name`:
# each read derives a compound datatype from the shared dataset, and `HDF5.jl` finalizes those from
# whichever thread runs the garbage collector, with no lock, while another thread is inside the library.
# Opening and closing within the locked region keeps every library call on one thread at a time.
#
# **Raising the chunk cache does not help either.** A NISAR band is stored in 512 x 512 chunks — 2 MiB of
# `ComplexF32`, shuffled and deflated — against a library default cache of 1 MiB, so no chunk is ever
# retained and overlapping windows re-inflate; but `H5Pset_chunk_cache` at 1 GiB reads back as 1 MiB from
# `H5Dget_access_plist` and leaves a re-read of ten overlapping 2048² windows at 1.45 s against 1.42 s. A
# reader that needs chunk reuse has to hold the samples itself.

"""
    nisar_samples_path(b::NisarBackend) -> String

Where `b`'s polarization's samples live: under `swaths` for a radar-geometry product, `grids` for a
geocoded one.

The polarization is the first this frequency lists, which is the only one a single-polarization product
carries and the co-polarized channel of a dual one.
"""
function nisar_samples_path(b::NisarBackend, polarization = nothing)
    group = b.product_type in GEOCODED_TYPES ? "grids" : "swaths"
    base = string(product_path(b), "/", group, "/frequency", b.frequency)
    pol = polarization === nothing ? nothing : String(polarization)
    return h5open(b.path, "r") do h
        g = _group(h, base)
        if pol === nothing
            haskey(g, "listOfPolarizations") || throw(ArgumentError(
                "`$(b.path)` lists no polarizations under $base"))
            pol = String(first(read(g["listOfPolarizations"])))
        end
        haskey(g, pol) || throw(ArgumentError(
            "`$(b.path)` has no $pol under $base; it lists " *
            join(String.(read(g["listOfPolarizations"])), ", ")))
        return string(base, "/", pol)
    end
end

function NisarRaster(b::NisarBackend, polarization = nothing)
    ds = nisar_samples_path(b, polarization)
    T, nx, ny = h5open(b.path, "r") do h
        d = h[ds]
        (eltype(d), size(d, 1), size(d, 2))
    end
    return NisarRaster{T}(b.path, ds, (ny, nx))
end

Base.size(r::NisarRaster) = r.dims

function Base.getindex(r::NisarRaster{T}, rows::AbstractUnitRange{<:Integer},
                       cols::AbstractUnitRange{<:Integer}) where {T}
    @boundscheck checkbounds(r, rows, cols)
    plan = _chunk_plan(r)
    isnothing(plan) || return _read_by_chunks!(similar(Matrix{T}, length(rows), length(cols)),
                                               r, plan, rows, cols)

    # The hyperslab under the lock, the transpose outside it: `permutedims` of a window is a third of the
    # cost and needs no library call, so holding `HDF5_IO` across it would serialize concurrent readers on
    # work that has nothing to serialize.
    #
    # **`::Matrix{T}` because `HDF5.Dataset`'s `getindex` is not inferrable**, and this method's return type
    # is what every consumer's loop is compiled against: without the annotation `Amplitude`'s window read
    # dispatches `_magnitude` once per sample, which measures 158 ns an element against 12 ns — a 2048²
    # window of a GSLC costing 0.71 s instead of 0.05 s. `T` is the dataset's own `eltype`, read when the
    # raster was constructed, so the assertion is the file's own type and not a coercion.
    raw::Matrix{T} = @lock HDF5_IO h5open(h -> h[r.dataset][cols, rows], r.path, "r")
    return permutedims(raw)
end

# ---------------------------------------------------------------------------
# Reading a window chunk by chunk
# ---------------------------------------------------------------------------

"""
    ChunkPlan{T}

How a dataset's samples are stored: the chunk shape, the value an unwritten chunk reads as, and the filters
to undo in read order.

`chunk` and `nchunks` are in the *file's* `HDF5.jl`-reported order — `(sample, line)` — because that is the
order a chunk's bytes arrive in and the order the offsets are computed from.
"""
struct ChunkPlan{T}
    chunk::Tuple{Int,Int}
    nchunks::Tuple{Int,Int}
    bytes::Int
    # **The dataset's own fill value, not zero.** A chunked dataset allocates a chunk when it is first
    # written, so a region never written has no bytes to read and reads as the fill value instead — which on
    # a NISAR GSLC is `NaN + NaN*im` over everything outside the imaged swath. Assuming zero turns a gap
    # into valid-looking black.
    fill::T
    # `(filter id, first client value)` in the order to *undo* them, which is the reverse of the order the
    # creation property list lists and therefore of the order they were applied on write.
    undo::Vector{Tuple{Int,Int}}
end

const H5Z_DEFLATE = 1
const H5Z_SHUFFLE = 2

# Resolved once per dataset, because it costs three library calls and never changes. `nothing` means the
# chunked path does not apply: an unchunked dataset, or one carrying a filter with no decoder here — in
# which case `H5Dread` is still correct, just serial.
const CHUNK_PLANS = Dict{Tuple{String,String},Union{ChunkPlan,Nothing}}()
const CHUNK_PLAN_LOCK = ReentrantLock()

function _chunk_plan(r::NisarRaster{T}) where {T}
    key = (r.path, r.dataset)
    @lock CHUNK_PLAN_LOCK begin
        haskey(CHUNK_PLANS, key) && return CHUNK_PLANS[key]
        plan = @lock HDF5_IO h5open(h -> _chunk_plan(h[r.dataset], T), r.path, "r")
        CHUNK_PLANS[key] = plan
        return plan
    end
end

function _chunk_plan(d::HDF5.Dataset, ::Type{T}) where {T}
    dcpl = HDF5.get_create_properties(d)
    # `dcpl.chunk` errors rather than returning `nothing` when there are no chunks, so the layout decides
    # whether to ask.
    dcpl.layout === :chunked || return nothing
    chunk = Tuple(Int.(dcpl.chunk))
    length(chunk) == 2 || return nothing
    nf = Int(HDF5.API.h5p_get_nfilters(dcpl))
    undo = Tuple{Int,Int}[]
    for i in (nf - 1):-1:0
        cd = Vector{Cuint}(undef, 8)
        ne = Ref{Csize_t}(length(cd))
        name = Vector{UInt8}(undef, 64)
        flags = Ref{Cuint}()
        config = Ref{Cuint}()
        id = Int(HDF5.API.h5p_get_filter(dcpl, i, flags, ne, cd, length(name), name, config))
        id in (H5Z_DEFLATE, H5Z_SHUFFLE) || return nothing
        push!(undo, (id, ne[] > 0 ? Int(cd[1]) : 0))
    end
    # The value an unwritten chunk reads as. `H5Pget_fill_value` converts into whatever type it is handed,
    # and `HDF5.datatype(d)` is the dataset's own — so a compound `ComplexF32` arrives as one.
    fill = Ref{T}(zero(T))
    HDF5.API.h5p_get_fill_value(dcpl, HDF5.datatype(d), fill)
    return ChunkPlan{T}(chunk, Tuple(cld.(size(d), chunk)), prod(chunk) * sizeof(T), fill[], undo)
end

"""
    _read_by_chunks!(out, r, plan, rows, cols) -> out

`r[rows, cols]` read a chunk at a time, with the decompression outside `HDF5_IO`.

**This is what lets a blocked reader use more than one core.** `H5Dread` runs the filter pipeline itself, so
a window read that goes through it inflates inside the lock and every other task waits — measured at 2.6 of
10 threads on a blocked correlation, with nine threads parked in `__psynch_cvwait` and one in
`inflate_fast`. `H5Dread_chunk` instead hands back a chunk's stored bytes with no filter applied, so the
serialized part is a byte read and the inflate is ordinary Julia work on any thread.

**Buffers are pooled and the stored bytes land in one array**, because the allocation is otherwise what
limits the concurrency this exists to deliver: a byte vector per chunk is 24 MiB per 2048² window and two
scratch buffers per task another 48 MiB, and ten block reads in flight turn that into half a gigabyte of
garbage per round.

Measured on a NISAR GSLC band, ten overlapping 2048² amplitude windows on ten threads: 1.312 s through
`H5Dread` — slower than reading them one at a time, which is what a contended lock looks like — against
0.126 s here, with concurrency now worth 1.60x rather than 0.66x. One window alone went from 0.707 s to
0.029 s, the larger part of that from making this method's return type inferrable.
"""
function _read_by_chunks!(out::Matrix{T}, r::NisarRaster{T}, plan::ChunkPlan{T},
                          rows::AbstractUnitRange{<:Integer},
                          cols::AbstractUnitRange{<:Integer}) where {T}
    cs, cl = plan.chunk                  # file order: samples, then lines
    # `(is, il)`: the chunk's index along samples and along lines. Flat, because the work is handed out by
    # linear index and a two-dimensional comprehension would make `pairs` yield `CartesianIndex` keys.
    jobs = vec([(is, il)
                for is in (fld(first(cols) - 1, cs) + 1):(fld(last(cols) - 1, cs) + 1),
                    il in (fld(first(rows) - 1, cl) + 1):(fld(last(rows) - 1, cl) + 1)])

    # Every chunk's stored bytes end to end, with `span[k]` naming chunk `k`'s slice and its filter mask. An
    # unallocated chunk gets an empty span, which is how the fill value is signalled.
    slack = _scratch_bytes(plan)
    bytes = _take_scratch(length(jobs) * slack)
    span = Vector{UnitRange{Int}}(undef, length(jobs))
    mask = Vector{UInt32}(undef, length(jobs))
    try
        # The locked region is the open, the chunk lookups and the byte reads — nothing else. One open for
        # the whole window rather than one per chunk.
        @lock HDF5_IO h5open(r.path, "r") do h
            d = h[r.dataset]
            pos = 1
            m = Ref{UInt32}(0)
            for (k, (is, il)) in pairs(jobs)
                off = HDF5.API.hsize_t[(il - 1) * cl, (is - 1) * cs]
                info = HDF5.API.h5d_get_chunk_info_by_coord(d, off)
                if info.size == 0
                    span[k] = pos:(pos - 1)          # empty: never written, so the fill value applies
                    continue
                end
                n = Int(info.size)
                # `h5d_read_chunk` writes from the start of whatever buffer it is given, so each chunk is
                # read into its own slot with a view rather than into a shared front.
                dest = view(bytes, pos:(pos + n - 1))
                HDF5.API.h5d_read_chunk(d, HDF5.API.H5P_DEFAULT, off, m, dest)
                mask[k] = m[]
                span[k] = pos:(pos + n - 1)
                pos += n
            end
        end

        ntasks = max(1, min(Threads.nthreads(), length(jobs)))
        @sync for t in 1:ntasks
            Threads.@spawn begin
                infl = _take_scratch(slack)
                plain = _take_scratch(slack)
                try
                    for k in t:ntasks:length(jobs)
                        _place_chunk!(out, plan, jobs[k], view(bytes, span[k]), mask[k],
                                      infl, plain, rows, cols)
                    end
                finally
                    _give_scratch(infl)
                    _give_scratch(plain)
                end
            end
        end
    finally
        _give_scratch(bytes)
    end
    return out
end

# **A pool, because these buffers are megabytes and a blocked reader asks for them thousands of times.** A
# 2048² window of a NISAR band needs 24 MiB for the stored bytes and 2.4 MiB per decoding task; allocating
# that per read made the garbage collector, not the inflate, the limit on how much of the machine a
# concurrent read could use.
#
# Bounded by how many reads are in flight at once rather than by a byte budget: a buffer is only ever held
# for the duration of one read, so the pool settles at the concurrency the caller actually uses. Buffers are
# handed out by "big enough", and a request larger than anything free allocates rather than growing one in
# place, so a single outsized read does not permanently inflate every slot.
const SCRATCH_LOCK = ReentrantLock()
const SCRATCH_FREE = Vector{Vector{UInt8}}()

function _take_scratch(nbytes::Integer)
    n = Int(nbytes)
    got = @lock SCRATCH_LOCK begin
        i = findfirst(b -> length(b) >= n, SCRATCH_FREE)
        isnothing(i) ? nothing : popat!(SCRATCH_FREE, i)
    end
    return isnothing(got) ? Vector{UInt8}(undef, n) : got
end

_give_scratch(b::Vector{UInt8}) = (@lock SCRATCH_LOCK push!(SCRATCH_FREE, b); nothing)

# Room for a chunk that did not compress. zlib's worst case on incompressible input is the input plus about
# a thousandth, and a stored chunk is never larger than that, so the whole chunk plus an eighth is slack
# enough to hold any intermediate the filter chain produces.
_scratch_bytes(plan::ChunkPlan) = plan.bytes + plan.bytes ÷ 8 + 1024

# One decoded chunk's contribution to `out`, which is indexed `(line, sample)` where the chunk's bytes are
# `(sample, line)` — so this is a transposing copy of the rectangle the two have in common.
function _place_chunk!(out::Matrix{T}, plan::ChunkPlan{T}, job::Tuple{Int,Int},
                       raw::AbstractVector{UInt8}, mask::UInt32,
                       infl::Vector{UInt8}, plain::Vector{UInt8},
                       rows::AbstractUnitRange{<:Integer},
                       cols::AbstractUnitRange{<:Integer}) where {T}
    is, il = job
    cs, cl = plan.chunk
    srows = intersect(rows, ((il - 1) * cl + 1):(il * cl))
    scols = intersect(cols, ((is - 1) * cs + 1):(is * cs))
    (isempty(srows) || isempty(scols)) && return out

    if isempty(raw)
        for s in scols, l in srows
            out[l - first(rows) + 1, s - first(cols) + 1] = plan.fill
        end
        return out
    end

    n = length(raw)
    copyto!(plain, 1, raw, 1, n)
    for (i, (id, cd)) in pairs(plan.undo)
        (mask >> (length(plan.undo) - i)) & 0x1 == 1 && continue
        if id == H5Z_DEFLATE
            n = _inflate!(infl, plain, n)
            plain, infl = infl, plain
        else                                      # H5Z_SHUFFLE
            # Whatever length it is handed, as the filter itself does: shuffle is not necessarily the first
            # filter a writer applied, so its input is not necessarily a whole chunk.
            _unshuffle!(infl, plain, cd == 0 ? sizeof(T) : cd, n)
            plain, infl = infl, plain
        end
    end
    n == plan.bytes || throw(ArgumentError(
        "a chunk decoded to $n bytes where the chunk is $(plan.bytes)"))

    got = reshape(reinterpret(T, view(plain, 1:plan.bytes)), cs, cl)
    # `out` is contiguous down its first axis, so the destination is walked in order and `got` strided.
    for s in scols
        j = s - first(cols) + 1
        gs = s - (is - 1) * cs
        for l in srows
            out[l - first(rows) + 1, j] = got[gs, l - (il - 1) * cl]
        end
    end
    return out
end

# zlib `uncompress` on the first `n` bytes of `src`, which is the stream HDF5's deflate filter writes —
# `compress2`, so zlib-wrapped with an adler32 rather than a raw deflate block.
function _inflate!(dst::Vector{UInt8}, src::Vector{UInt8}, n::Integer)
    len = Ref{Csize_t}(length(dst))
    rc = ccall((:uncompress, Zlib_jll.libz), Cint,
               (Ptr{UInt8}, Ptr{Csize_t}, Ptr{UInt8}, Csize_t),
               dst, len, src, Csize_t(n))
    rc == 0 || throw(ArgumentError("zlib could not inflate a chunk: uncompress returned $rc"))
    return Int(len[])
end

# HDF5's shuffle filter groups the first byte of every element, then the second, and so on. The trailing
# `n % elsize` bytes are left where they are, which is what the filter does with them.
function _unshuffle!(dst::Vector{UInt8}, src::Vector{UInt8}, elsize::Integer, n::Integer)
    es = Int(elsize)
    nelem = Int(n) ÷ es
    k = 1
    for b in 1:es
        j = b
        for _ in 1:nelem
            dst[j] = src[k]
            j += es
            k += 1
        end
    end
    for i in (nelem * es + 1):Int(n)
        dst[i] = src[i]
    end
    return dst
end

Base.getindex(r::NisarRaster, i::Int, j::Int) = r[i:i, j:j][1, 1]
Base.getindex(r::NisarRaster, rows::AbstractUnitRange{<:Integer}, j::Int) = r[rows, j:j][:, 1]
Base.getindex(r::NisarRaster, i::Int, cols::AbstractUnitRange{<:Integer}) = r[i:i, cols][1, :]

read_pixels(b::NisarBackend) = NisarRaster(b)

# ---------------------------------------------------------------------------
# Geocoded products
# ---------------------------------------------------------------------------

"""
    GeocodedGrid

The map grid a geocoded product's samples lie on: `origin` at the outer corner of the first pixel,
signed `spacing`, `size` as `(rows, columns)`, and the EPSG code of the projection.

**`origin` is the corner, not the first pixel's center.** The product stores `xCoordinates` and
`yCoordinates` as centers; a geotransform names the corner, and a consumer intersecting two grids or
handing one to GDAL wants the latter.
"""
struct GeocodedGrid
    origin::Tuple{Float64,Float64}
    spacing::Tuple{Float64,Float64}
    size::Tuple{Int,Int}
    epsg::Int
end

"""
    GeocodedProduct <: AbstractSLC

A geocoded NISAR product: its identification, the map grid its samples lie on, and the samples.

Separate from [`SLC`](@ref) because it has no slant-range geometry to describe — `read_geometry` refuses
one — while still carrying an acquisition's identification and samples. [`pixels`](@ref) reads it.
"""
struct GeocodedProduct{B<:AbstractSLCBackend} <: AbstractSLC
    backend::B
    identification::Identification
    grid::GeocodedGrid
end

"""
    geocoded_grid(b::NisarBackend) -> GeocodedGrid

The map grid of a geocoded product's frequency group.
"""
function geocoded_grid(b::NisarBackend)
    b.product_type in GEOCODED_TYPES || throw(ArgumentError(
        "`$(b.path)` is a $(b.product_type) product, which lies in radar geometry and so has no map " *
        "grid; `read_geometry` describes it instead"))
    base = string(product_path(b), "/grids/frequency", b.frequency)
    return h5open(b.path, "r") do h
        g = _group(h, base)
        x = read(_dataset(g, "xCoordinates"))
        y = read(_dataset(g, "yCoordinates"))
        dx = _scalar(g, "xCoordinateSpacing")
        dy = _scalar(g, "yCoordinateSpacing")
        return GeocodedGrid((first(x) - dx / 2, first(y) - dy / 2), (dx, dy),
                            (length(y), length(x)), _scalar_int(g, "projection"))
    end
end

"""
    open_geocoded(path; frequency = nothing) -> GeocodedProduct

Open a geocoded NISAR product — a GSLC, GCOV, GUNW or GOFF.

[`open_slc`](@ref) refuses these: they carry no slant-range axis, so there is no `RadarGeometry` to read.
What they do carry is a map grid and samples, which is what this returns.
"""
function open_geocoded(path::AbstractString; frequency = nothing)
    ispath(path) || throw(ArgumentError("`$path` is not a readable file"))
    ishdf5(path) || throw(ArgumentError("`$path` is not an HDF5 file"))
    return h5open(path, "r") do h
        band = nisar_band(h)
        product_type = nisar_product_type(h, band)
        product_type in GEOCODED_TYPES || throw(ArgumentError(
            "`$path` is a $product_type product, which lies in radar geometry; use `open_slc`"))
        freq = frequency === nothing ? default_frequency(h, band) : String(frequency)
        b = NisarBackend(path, band, product_type, freq)
        return GeocodedProduct(b, read_identification(b, h), geocoded_grid(b))
    end
end

grid(g::GeocodedProduct) = g.grid
nlines(g::GeocodedProduct) = g.grid.size[1]
nsamples(g::GeocodedProduct) = g.grid.size[2]
