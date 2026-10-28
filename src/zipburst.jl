# Bursts read straight out of a Sentinel-1 SLC zip, the way ASF's burst extractor cuts them.
#
# A subswath's measurement raster stacks its bursts in one uncompressed TIFF, one strip per line, the
# samples of every line contiguous from the first line to the last. Burst `k` of a subswath with `L` lines
# and `S` samples per burst is therefore the `4 L S` bytes from `D + 4 (k - 1) L S`, where `D` is where the
# samples begin. The zip stores the raster deflated, so those bytes are reached by inflating the member
# from its start, and the inflater stops once the last burst asked for has gone by
# (`ChunkManifests.inflate_ranges!`): a burst early in a subswath costs only the bytes before it.
#
# Each burst is written to a file of its own, a TIFF of the layout the extractor serves, so a product
# assembled from these bursts reads them exactly as it reads extracted ones.

"""
    ZipBurst(zip, swath, polarization, burst; dir = nothing) -> ZipBurst

Burst `burst` (1-based) of subswath `swath` of the Sentinel-1 IW SLC whose `.zip` is at `zip` — a local
path, an `https://` URL or an `s3://` URI — in `polarization`.

The raster is inflated out of the zip into `dir` when it is first read, a temporary directory by
default; a raster already there is not read again. [`zip_bursts`](@ref) lists every burst of a zip.
"""
struct ZipBurst <: AbstractSLCSource
    zip::String
    slc::String
    swath::Int
    polarization::String
    burst::Int
    dir::Union{Nothing,String}
end

function ZipBurst(zip::AbstractString, swath::Integer, polarization::AbstractString, burst::Integer;
                  dir = nothing)
    swath in 1:3 || throw(ArgumentError("Sentinel-1 IW subswaths are 1, 2 and 3; got $swath"))
    burst >= 1 || throw(ArgumentError("bursts are numbered from 1, so burst $burst does not exist"))
    pol = uppercase(String(polarization))
    pol in ("HH", "HV", "VH", "VV") || throw(ArgumentError(
        "a Sentinel-1 polarization is HH, HV, VH or VV; got `$polarization`"))
    return ZipBurst(String(zip), _zip_product_name(zip), Int(swath), pol, Int(burst),
                    dir === nothing ? nothing : String(dir))
end

const _ZIP_BURST_SPEC = r"^(.*\.zip)#IW([1-3])/(HH|HV|VH|VV)/(\d+)$"

"""
    ZipBurst(spec; dir = nothing) -> ZipBurst
    zip_burst_spec(b::ZipBurst) -> String

A burst named as `<zip>#IW<n>/<POL>/<index>`: the zip's location, then the subswath, polarization and
index the extractor's URL for the same burst ends with — counted from zero, as there.
"""
function ZipBurst(spec::AbstractString; dir = nothing)
    m = match(_ZIP_BURST_SPEC, spec)
    m === nothing && throw(ArgumentError(
        "`$spec` does not name a burst of an SLC zip; one reads `<SLC zip>#IW<1-3>/<HH|HV|VH|VV>/<index>`"))
    return ZipBurst(something(m.captures[1]), parse(Int, something(m.captures[2])), something(m.captures[3]),
                    parse(Int, something(m.captures[4])) + 1; dir)
end

zip_burst_spec(b::ZipBurst) = b.zip * "#IW" * string(b.swath) * "/" * b.polarization * "/" * string(b.burst - 1)

"""
    is_zip_burst_spec(s) -> Bool

Whether `s` names a burst of an SLC zip as [`ZipBurst`](@ref) reads one.
"""
is_zip_burst_spec(s::AbstractString) = occursin(_ZIP_BURST_SPEC, s)

# The product a zip holds, from its file name: `S1A_IW_SLC__1SDV_…_D6DE.zip` holds `S1A_IW_SLC__1SDV_…_D6DE`.
function _zip_product_name(zip::AbstractString)
    base = String(zip[(_last_slash(zip) + 1):end])
    m = match(r"^(S1[A-D]_IW_SLC__[0-9A-Z_]+)\.zip$", base)
    m === nothing && throw(ArgumentError(
        "`$zip` is not a Sentinel-1 IW SLC zip, which is named `<product>.zip` after the product it holds"))
    return String(something(m.captures[1]))
end

_zip_burst_label(b::ZipBurst) = "burst " * string(b.burst) * " of IW" * string(b.swath) * " " * b.polarization *
                                " of `" * b.zip * "`"

# Named as an extracted burst's raster is (`_asf_local`), with the burst counted from zero.
function _zip_local(b::ZipBurst)
    file = lowercase(b.polarization) * "-" * string(b.burst - 1) * ".tiff"
    return joinpath(b.dir === nothing ? _scratch_dir() : b.dir, b.slc * "-iw" * string(b.swath) * "-" * file)
end

"""
    zip_rasters(bursts::AbstractVector{ZipBurst}) -> Vector{String}

The local path of each burst's raster, inflating those not already in their directories. Each subswath
member is inflated once for all its bursts, from its start to the end of the last of them, and the members
are inflated concurrently.
"""
function zip_rasters(bursts::AbstractVector{ZipBurst})
    paths = [_zip_local(b) for b in bursts]
    todo = Int[]
    waits = Base.Event[]
    @lock _DOWNLOADING_LOCK for k in eachindex(bursts)
        _downloaded(paths[k]) && continue
        if haskey(_DOWNLOADING, paths[k])
            push!(waits, _DOWNLOADING[paths[k]])
        else
            _DOWNLOADING[paths[k]] = Base.Event()
            push!(todo, k)
        end
    end
    try
        groups = Dict{Tuple{String,Int,String},Vector{Int}}()
        for k in todo
            push!(get!(Vector{Int}, groups, (bursts[k].zip, bursts[k].swath, bursts[k].polarization)), k)
        end
        tasks = [Threads.@spawn _inflate_bursts(bursts[ks], paths[ks]) for ks in values(groups)]
        # Every task has stopped writing before an error is reported.
        failed = nothing
        for t in tasks
            try
                wait(t)
            catch err
                failed === nothing && (failed = err)
            end
        end
        failed === nothing || throw(failed)
    finally
        @lock _DOWNLOADING_LOCK for k in todo
            notify(pop!(_DOWNLOADING, paths[k]))
        end
    end
    foreach(wait, waits)
    for p in paths
        _downloaded(p) || throw(ErrorException("inflating `$p` failed in the task that started it"))
    end
    return paths
end

# How much of a measurement member's head is read for its layout: its first directory and the strip
# tables that directory's values point to, a few hundred kilobytes at most.
const _ZIP_TIFF_HEAD = 1 << 20

# Bursts of one subswath member, written to `dests`.
function _inflate_bursts(bs::Vector{ZipBurst}, dests::Vector{String})
    b = first(bs)
    r = ZipArchive(b.zip)
    m = _member(r, _zip_measurement_member(r, b.swath, b.polarization))
    a = zip_subswath(b)
    L, S, N = a.lines_per_burst, a.samples_per_burst, nbursts(a)
    head = read_member(r.source, m, 0, min(_ZIP_TIFF_HEAD, m.uncompressed); transport = r.transport)
    D = _measurement_data_offset(head, L, S, N, b.zip * " " * m.name)
    order = sortperm([x.burst for x in bs])
    nb = 4 * L * S
    for k in order
        1 <= bs[k].burst <= N || throw(ArgumentError(
            "subswath IW$(b.swath) of `$(b.zip)` has $N bursts, so burst $(bs[k].burst) does not exist"))
    end
    ranges = [(D + (bs[k].burst - 1) * nb, nb) for k in order]
    parts = [dests[k] * ".part" for k in order]
    ios = IOStream[]
    try
        for p in parts
            mkpath(dirname(p))
            io = open(p, "w")
            push!(ios, io)
            _write_burst_tiff_header(io, L, S)
        end
        inflate_ranges!(ios, r.source, m, ranges; transport = r.transport)
    finally
        foreach(close, ios)
    end
    for (j, k) in enumerate(order)
        mv(parts[j], dests[k]; force = true)
    end
    return dests
end

function _zip_measurement_member(r::ZipArchive, swath::Integer, polarization::AbstractString)
    id_str = _swath_id(swath, polarization)
    hits = sort!([m.name for m in r.members
                  if occursin("/measurement/", m.name) && occursin(id_str, m.name) && endswith(m.name, ".tiff")])
    isempty(hits) && throw(ArgumentError(
        "`$(r.uri)` has no measurement raster for subswath IW$swath polarization $(uppercase(polarization))"))
    return only(hits)
end

_le32(b::Vector{UInt8}, i::Int) = Int(b[i]) | Int(b[i + 1]) << 8 | Int(b[i + 2]) << 16 | Int(b[i + 3]) << 24
_le16(b::Vector{UInt8}, i::Int) = Int(b[i]) | Int(b[i + 1]) << 8

"""
    _measurement_data_offset(head, L, S, N, context) -> Int

Where the samples of a measurement raster begin, from the first `head` bytes of the file, for a subswath
of `N` bursts of `L` lines of `S` samples.

The file's header points to its directory, which follows the samples, a byte of padding at most past
their end. Products also carry a copy of the directory at byte 8, followed by the out-of-line values it
points to: the samples begin at its first strip offset when that is set, and where those values end when
the copy's strip table is zeroed. Either way the result must end the samples where the header's
directory begins, which is checked, as are the copy's width and height against the annotation's.
"""
function _measurement_data_offset(head::Vector{UInt8}, L::Int, S::Int, N::Int, context::AbstractString)
    length(head) >= 16 && head[1] == UInt8('I') && head[2] == UInt8('I') && _le16(head, 3) == 42 ||
        throw(ArgumentError("$context is not a little-endian classic TIFF"))
    ifd = _le32(head, 5)
    n = _le16(head, 9)
    fail() = throw(ArgumentError(
        "$context is not laid out as a Sentinel-1 measurement raster of $N bursts of $L lines of $S samples: " *
        "its samples could not be located without reading its directory, which follows them"))
    9 + 2 + 12n + 4 <= length(head) || fail()
    sizes = (1, 1, 2, 4, 8, 1, 1, 2, 4, 8, 4, 8)
    width = height = first_strip = 0
    values_end = 8 + 2 + 12n + 4
    for i in 0:(n - 1)
        e = 9 + 2 + 12i
        tag, typ, count, value = _le16(head, e), _le16(head, e + 2), _le32(head, e + 4), _le32(head, e + 8)
        1 <= typ <= length(sizes) || fail()
        bytes = sizes[typ] * count
        if bytes > 4
            value + bytes <= length(head) || fail()
            values_end = max(values_end, value + bytes)
        end
        if tag == 256
            width = typ == 3 ? value & 0xffff : value
        elseif tag == 257
            height = typ == 3 ? value & 0xffff : value
        elseif tag == 273
            typ == 4 || fail()
            first_strip = count == 1 ? value : _le32(head, value + 1)
        end
    end
    (width == S && height == N * L) || fail()
    D = first_strip > 0 ? first_strip : values_end
    gap = ifd - (D + 4 * N * L * S)
    0 <= gap <= 1 || fail()
    return D
end

# The header a burst's TIFF starts with: one directory of the tags `open_measurement` reads, a strip a
# line, the strip tables, then the `L` lines of `S` complex 16-bit samples that follow it.
function _write_burst_tiff_header(io::IO, L::Int, S::Int)
    nentries = 11
    dir_at = 8
    offsets_at = dir_at + 2 + 12 * nentries + 4
    counts_at = offsets_at + 4L
    data_at = counts_at + 4L
    w16(x) = write(io, htol(UInt16(x)))
    w32(x) = write(io, htol(UInt32(x)))
    entry(tag, typ, count, value) = (w16(tag); w16(typ); w32(count); w32(value))
    write(io, UInt8('I'), UInt8('I'))
    w16(42)
    w32(dir_at)
    w16(nentries)
    entry(256, 4, 1, S)              # ImageWidth
    entry(257, 4, 1, L)              # ImageLength
    entry(258, 3, 1, 32)             # BitsPerSample
    entry(259, 3, 1, 1)              # Compression: none
    entry(262, 3, 1, 1)              # PhotometricInterpretation
    entry(273, 4, L, offsets_at)     # StripOffsets
    entry(277, 3, 1, 1)              # SamplesPerPixel
    entry(278, 4, 1, 1)              # RowsPerStrip
    entry(279, 4, L, counts_at)      # StripByteCounts
    entry(284, 3, 1, 1)              # PlanarConfiguration
    entry(339, 3, 1, 5)              # SampleFormat: complex signed integer
    w32(0)
    for i in 0:(L - 1)
        w32(data_at + 4S * i)
    end
    for _ in 1:L
        w32(4S)
    end
    return io
end

