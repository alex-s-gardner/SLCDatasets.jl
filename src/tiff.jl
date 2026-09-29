# Reading the measurement raster of a Sentinel-1 product.
#
# A measurement TIFF is uncompressed, one strip per line, one complex-int16 sample per pixel, and its
# strip offsets are tabulated. So a line is a contiguous run of bytes at a known offset and the file
# can be memory-mapped and indexed rather than decoded: no TIFF library is needed, and no more of the
# raster is touched than is asked for. The same layout holds for a `.SAFE` measurement file and for a
# single-burst raster from ASF's burst extractor.
#
# The layout is not assumed. Every property the reader relies on is a tag in the file, so each is read
# and checked, and a raster written some other way is refused by name rather than misread: a tiled or
# compressed file has no per-line offsets to index, and the values would be garbage rather than absent.

const TIFF_LITTLE_ENDIAN = 0x4949  # "II"
const TIFF_BIG_ENDIAN = 0x4d4d     # "MM"

const TIFF_MAGIC = 42
const TIFF_BIGTIFF_MAGIC = 43

# The tags read below. A measurement file carries others — image description, datetime, the GDAL
# metadata a burst raster picks up — and those are irrelevant to addressing pixels.
const TIFF_IMAGE_WIDTH = 256
const TIFF_IMAGE_LENGTH = 257
const TIFF_BITS_PER_SAMPLE = 258
const TIFF_COMPRESSION = 259
const TIFF_STRIP_OFFSETS = 273
const TIFF_SAMPLES_PER_PIXEL = 277
const TIFF_ROWS_PER_STRIP = 278
const TIFF_STRIP_BYTE_COUNTS = 279
const TIFF_PLANAR_CONFIGURATION = 284
const TIFF_TILE_WIDTH = 322
const TIFF_SAMPLE_FORMAT = 339

const TIFF_COMPRESSION_NONE = 1
const TIFF_PLANAR_CONTIGUOUS = 1
const TIFF_SAMPLE_FORMAT_COMPLEX_INT = 5

# Field types, by the code stored in a directory entry. Only the widths needed to read the tags above
# are listed; an entry of any other type is refused where it is read.
const TIFF_TYPE_SIZES = Dict{Int,Int}(1 => 1, 2 => 1, 3 => 2, 4 => 4, 5 => 8, 7 => 1,
                                      8 => 2, 9 => 4, 10 => 8, 11 => 4, 12 => 8, 16 => 8, 17 => 8)

"""
    TiffField

One entry of a TIFF image file directory: what the tag is, and where its values are.

An entry's last four bytes are either the values themselves or the offset to them, depending on their
total size, so both readings are kept: `payload` is those bytes as an integer, and `valueoffset` is
where they sit in the file. Which one applies is resolved when the values are read.
"""
struct TiffField
    tag::Int
    typ::Int
    count::Int
    payload::UInt32
    valueoffset::Int
end

"""
    StripedTiff{T} <: DiskArrays.AbstractDiskArray{T,2}

The raster of an uncompressed, one-line-per-strip TIFF, memory-mapped and addressed by line.

Indexing reads from the mapping, so a window costs the lines it covers and nothing else. Build one
with [`open_tiff`](@ref), which checks the file has this layout.

Indexing a range of lines and samples — `t[rows, cols]` — is the operation to prefer: it copies each
line's run in one go, where scalar indexing pays the address arithmetic per pixel.

A `DiskArrays` array, so views, broadcasting and strided indexing behave as they do for every other raster
this package returns. `eachchunk` reports the strips — one line each, which is the file's real granularity.
"""
struct StripedTiff{T,D<:AbstractVector{UInt8}} <: DiskArrays.AbstractDiskArray{T,2}
    path::String
    data::D
    # Byte offset and length of each line's samples, in file order. Parallel to `1:nlines`.
    offsets::Vector{Int64}
    bytecounts::Vector{Int64}
    nlines::Int
    nsamples::Int
    swapped::Bool
end

"""
    open_tiff(path) -> StripedTiff

Memory-map the TIFF at `path` and read its directory.

Throws unless the file is a classic TIFF whose raster is uncompressed, stored one line per strip, with
a single complex 16-bit integer sample per pixel — the layout a Sentinel-1 measurement file has. The
raster itself is not read.
"""
function open_tiff(path::AbstractString)
    isfile(path) || throw(ArgumentError("`$path` is not a file, so it holds no TIFF raster"))
    data = open(io -> mmap(io, Vector{UInt8}, filesize(io)), path, "r")
    return StripedTiff(String(path), data)
end

function StripedTiff(path::AbstractString, data::AbstractVector{UInt8})
    length(data) >= 8 || throw(ArgumentError(
        "`$path` is $(length(data)) bytes, too short to hold a TIFF header"))

    byteorder = (UInt16(data[1]) << 8) | UInt16(data[2])
    swapped = if byteorder == TIFF_LITTLE_ENDIAN
        !(ENDIAN_BOM == 0x04030201)
    elseif byteorder == TIFF_BIG_ENDIAN
        ENDIAN_BOM == 0x04030201
    else
        throw(ArgumentError("`$path` does not start with a TIFF byte-order mark (`II` or `MM`), " *
                            "so it is not a TIFF file"))
    end

    magic = _tiff_load(UInt16, data, 3, swapped)
    if magic == TIFF_BIGTIFF_MAGIC
        throw(ArgumentError(
            "`$path` is a BigTIFF, whose directory this reader does not read. Sentinel-1 " *
            "measurement rasters are classic TIFF; a BigTIFF is written by something else"))
    end
    magic == TIFF_MAGIC || throw(ArgumentError(
        "`$path` has TIFF version $magic rather than $TIFF_MAGIC, so it is not a TIFF file"))

    fields = _read_directory(path, data, swapped)
    return _striped_tiff(path, data, fields, swapped)
end

# One entry is 12 bytes: tag, type, count, then four bytes of values or an offset to them.
function _read_directory(path, data, swapped)
    offset = Int(_tiff_load(UInt32, data, 5, swapped)) + 1
    _tiff_inbounds(path, data, offset, 2, "the image file directory")
    n = Int(_tiff_load(UInt16, data, offset, swapped))
    n > 0 || throw(ArgumentError("`$path` has an empty TIFF directory, so it names no raster"))
    _tiff_inbounds(path, data, offset + 2, 12 * n, "the image file directory")

    fields = Dict{Int,TiffField}()
    for i in 0:(n - 1)
        e = offset + 2 + 12 * i
        tag = Int(_tiff_load(UInt16, data, e, swapped))
        fields[tag] = TiffField(tag,
                                Int(_tiff_load(UInt16, data, e + 2, swapped)),
                                Int(_tiff_load(UInt32, data, e + 4, swapped)),
                                _tiff_load(UInt32, data, e + 8, swapped),
                                e + 8)
    end
    return fields
end

function _striped_tiff(path, data, fields, swapped)
    # Tiling and striping are alternatives, so a tiled file has no strip offsets at all. Naming the
    # tiling is a better error than reporting the strip tags missing.
    haskey(fields, TIFF_TILE_WIDTH) && throw(ArgumentError(
        "`$path` stores its raster in tiles rather than strips, which this reader does not read. " *
        "Sentinel-1 measurement rasters are striped one line per strip"))

    nsamples = _tiff_scalar(path, data, fields, TIFF_IMAGE_WIDTH, swapped)
    nlines = _tiff_scalar(path, data, fields, TIFF_IMAGE_LENGTH, swapped)

    _tiff_require(path, data, fields, TIFF_COMPRESSION, TIFF_COMPRESSION_NONE, swapped,
                  "is compressed, so its lines are not addressable without decoding it")
    _tiff_require(path, data, fields, TIFF_ROWS_PER_STRIP, 1, swapped,
                  "does not store one line per strip, so a line is not one addressable run")
    _tiff_require(path, data, fields, TIFF_SAMPLES_PER_PIXEL, 1, swapped,
                  "does not carry exactly one sample per pixel")
    _tiff_require(path, data, fields, TIFF_BITS_PER_SAMPLE, 32, swapped,
                  "does not carry 32-bit samples, so they are not complex 16-bit integers")
    _tiff_require(path, data, fields, TIFF_SAMPLE_FORMAT, TIFF_SAMPLE_FORMAT_COMPLEX_INT, swapped,
                  "does not carry complex-integer samples")
    # Meaningless with one sample per pixel, and absent from some writers' output; checked only when
    # present, since a value naming a planar layout would contradict the sample count.
    if haskey(fields, TIFF_PLANAR_CONFIGURATION)
        _tiff_require(path, data, fields, TIFF_PLANAR_CONFIGURATION, TIFF_PLANAR_CONTIGUOUS, swapped,
                      "does not store its samples contiguously")
    end

    offsets = _tiff_values(path, data, fields, TIFF_STRIP_OFFSETS, swapped)
    bytecounts = _tiff_values(path, data, fields, TIFF_STRIP_BYTE_COUNTS, swapped)

    length(offsets) == nlines || throw(ArgumentError(
        "`$path` tabulates $(length(offsets)) strip offsets for $nlines lines; with one line per " *
        "strip there must be one offset per line"))
    length(bytecounts) == nlines || throw(ArgumentError(
        "`$path` tabulates $(length(bytecounts)) strip byte counts for $nlines lines; with one " *
        "line per strip there must be one per line"))

    # A line is `nsamples` complex-int16 samples. A file whose byte counts say otherwise is laid out
    # differently than its tags claim, and indexing it would silently read across line boundaries.
    #
    # The loop reports nothing itself: building a message per line would cost more than the check, and
    # a raster has as many lines as it is tall. It finds the first line that is wrong and leaves the
    # complaining to `_bad_strip`.
    expected = 4 * nsamples
    limit = length(data) - expected + 1
    bad = 0
    for i in eachindex(bytecounts, offsets)
        if bytecounts[i] != expected || !(1 <= offsets[i] + 1 <= limit)
            bad = i
            break
        end
    end
    bad == 0 || _bad_strip(path, data, offsets, bytecounts, bad, nsamples, expected)

    return StripedTiff{Complex{Int16},typeof(data)}(String(path), data, offsets, bytecounts,
                                                    nlines, nsamples, swapped)
end

@noinline _tiff_truncated(path, what) = throw(ArgumentError(
    "`$path` ends before $what, so the file is truncated"))

# Which of the two things is wrong with the line the scan stopped on.
@noinline function _bad_strip(path, data, offsets, bytecounts, i, nsamples, expected)
    bytecounts[i] == expected || throw(ArgumentError(
        "`$path` gives line $i a strip of $(bytecounts[i]) bytes, but $nsamples complex 16-bit " *
        "integer samples occupy $expected"))
    _tiff_inbounds(path, data, Int(offsets[i]) + 1, expected, "line $i")
end

function _tiff_inbounds(path, data, offset, len, what)
    (offset >= 1 && offset + len - 1 <= length(data)) || _tiff_truncated(path, what)
    return nothing
end

# `offset` is a 1-based byte index into `data`. The bytes are assembled rather than loaded as a `T`:
# nothing in a TIFF is required to sit at a multiple of its own width, and a real measurement raster puts
# its strips at offsets that are not.
function _tiff_load(::Type{T}, data, offset, swapped) where {T}
    x = zero(T)
    if swapped
        for k in 0:(sizeof(T) - 1)
            x = (x << 8) | T(data[offset + k])
        end
    else
        for k in (sizeof(T) - 1):-1:0
            x = (x << 8) | T(data[offset + k])
        end
    end
    return x
end

@noinline _tiff_no_tag(path, tag) = throw(ArgumentError(
    "`$path` has no TIFF tag $tag, so it does not describe a readable raster"))

function _tiff_scalar(path, data, fields, tag, swapped)
    haskey(fields, tag) || _tiff_no_tag(path, tag)
    vals = _tiff_values(path, data, fields, tag, swapped)
    length(vals) == 1 || throw(ArgumentError(
        "`$path` gives TIFF tag $tag $(length(vals)) values where one is expected"))
    return Int(only(vals))
end

function _tiff_require(path, data, fields, tag, expected, swapped, complaint)
    haskey(fields, tag) || _tiff_no_tag(path, tag)
    got = _tiff_scalar(path, data, fields, tag, swapped)
    got == expected || throw(ArgumentError(
        "`$path` $complaint: TIFF tag $tag is $got rather than $expected"))
    return nothing
end

# The values of one field, as `Int64`. Values totalling four bytes or fewer are stored in the entry
# itself; anything longer is at the offset the entry holds.
function _tiff_values(path, data, fields, tag, swapped)
    haskey(fields, tag) || _tiff_no_tag(path, tag)
    f = fields[tag]
    size = get(TIFF_TYPE_SIZES, f.typ, 0)
    size == 0 && throw(ArgumentError(
        "`$path` gives TIFF tag $tag field type $(f.typ), which this reader does not read"))
    total = size * f.count

    # Values fitting in four bytes sit in the entry itself, at its last four bytes; longer runs sit at
    # the offset those bytes hold. Either way each value is read at its own offset, since a value
    # narrower than the payload occupies the *first* bytes of it and byte-swapping the payload as a
    # whole would move it.
    base = total <= 4 ? f.valueoffset : Int(f.payload) + 1
    _tiff_inbounds(path, data, base, total, "the values of TIFF tag $tag")

    out = Vector{Int64}(undef, f.count)
    for i in 1:f.count
        out[i] = _tiff_field_value(data, f.typ, base + (i - 1) * size, swapped)
    end
    return out
end

_tiff_field_value(data, typ, offset, swapped) =
    typ == 1 || typ == 2 || typ == 7 ? Int64(data[offset]) :
    typ == 3 ? Int64(_tiff_load(UInt16, data, offset, swapped)) :
    typ == 8 ? Int64(_tiff_load(Int16, data, offset, swapped)) :
    typ == 4 ? Int64(_tiff_load(UInt32, data, offset, swapped)) :
    typ == 9 ? Int64(_tiff_load(Int32, data, offset, swapped)) :
    typ == 16 ? Int64(_tiff_load(UInt64, data, offset, swapped)) :
    typ == 17 ? _tiff_load(Int64, data, offset, swapped) :
    throw(ArgumentError("a TIFF field of type $typ is not readable as an integer"))

Base.size(t::StripedTiff) = (t.nlines, t.nsamples)
Base.IndexStyle(::Type{<:StripedTiff}) = IndexCartesian()

"""
    path(t::StripedTiff) -> String

The file the raster is mapped from.
"""
path(t::StripedTiff) = t.path

# The byte index of sample `j` of line `i`.
_sample_offset(t::StripedTiff, i::Int, j::Int) = Int(t.offsets[i]) + 4 * (j - 1) + 1

function Base.getindex(t::StripedTiff{T}, i::Int, j::Int) where {T}
    @boundscheck checkbounds(t, i, j)
    off = _sample_offset(t, i, j)
    re = _tiff_load(Int16, t.data, off, t.swapped)
    im = _tiff_load(Int16, t.data, off + 2, t.swapped)
    return T(re, im)
end

# A window, copied line by line. Each line's samples are contiguous, so this is one `copyto!` per line
# rather than the address arithmetic scalar indexing pays per pixel — the difference between a
# reasonable read and an unreasonable one when a caller pulls blocks out of a multi-gigabyte raster.
DiskArrays.readblock!(t::StripedTiff, dest::AbstractArray, rows::AbstractUnitRange{<:Integer},
                      cols::AbstractUnitRange{<:Integer}) = _copy_window!(dest, t, rows, cols)

# **The strips, which are the file's real granularity: one line each.** Reported rather than hidden behind
# `Unchunked`, even though a partial line here costs proportionally less than a whole one — the strips are
# uncompressed and the file is memory-mapped, so nothing forces a chunk to be read whole. Declaring them
# costs nothing: measured on a 19,717-sample subswath, a 512² window reads at 1.7 ns a sample under
# `Chunked` against 1.8 under `Unchunked`, because `DiskArrays`' batch strategy takes a dense sub-range
# directly rather than expanding it to chunk bounds.
DiskArrays.haschunks(::StripedTiff) = DiskArrays.Chunked()
DiskArrays.eachchunk(t::StripedTiff) = DiskArrays.GridChunks(t, (1, size(t, 2)))

# A line's samples are contiguous in the file but a row of the result is not, since the result is
# column-major. So the window is gathered line by line into a transposed buffer — where each line *is*
# contiguous — and transposed once at the end, which lets both halves run at memory bandwidth instead of
# writing every sample to a strided address. Measured on a real subswath, that is a fifth to a third off
# the read for any window wide enough to matter, and the buffer is the same size as the result.
#
# The bytes go through a buffer rather than being read from the mapping in place because a strip begins
# wherever the writer put it, which in a real product is not a multiple of the sample size; a sample
# straddling that boundary cannot be loaded directly.
function _copy_window!(out, t::StripedTiff{T}, rows, cols) where {T}
    (isempty(rows) || isempty(cols)) && return out
    ncols = length(cols)
    first_col = Int(first(cols))
    data = t.data

    gathered = Matrix{T}(undef, ncols, length(rows))
    bytes = reinterpret(UInt8, gathered)
    linebytes = 4 * ncols
    for (di, i) in enumerate(rows)
        copyto!(bytes, (di - 1) * linebytes + 1, data, _sample_offset(t, Int(i), first_col), linebytes)
    end
    if t.swapped
        for k in eachindex(gathered)
            v = gathered[k]
            gathered[k] = T(bswap(real(v)), bswap(imag(v)))
        end
    end
    permutedims!(out, gathered, (2, 1))
    return out
end
