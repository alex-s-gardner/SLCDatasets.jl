# ASF's burst extractor: where a burst's two files are, and fetching them.
#
# The extractor serves a Sentinel-1 SLC as its individual bursts: `{i}.tiff` is one burst's raster and
# `{i}.xml` its metadata, both indexed from zero at `sentinel1-burst.asf.alaska.edu/{slc}/{swath}/{pol}/`.
# A request is redirected through Earthdata Login, authenticated from `~/.netrc`, to a pre-signed S3 URL
# valid for an hour that serves `Range` requests. While the extractor is still generating a burst it
# answers `202` instead.
#
# The metadata file is not a slice of the annotation: it wraps the whole parent product's annotation for
# every subswath, and is the same file for every burst of that product and polarization. So it is fetched
# once per product and polarization and kept in memory. A raster, a couple of hundred megabytes read
# whole by any consumer, is downloaded to a directory and read from there.

"""
    ASF_BURST_HOST

The host ASF's burst extractor serves from.
"""
const ASF_BURST_HOST = "sentinel1-burst.asf.alaska.edu"

"""
    AsfBurst(slc, swath, polarization, burst; dir = nothing) -> AsfBurst
    AsfBurst(url; dir = nothing) -> AsfBurst

One burst of a Sentinel-1 SLC, from ASF's burst extractor.

`slc` is the parent granule's name, `burst` its **1-based** index within the subswath — the extractor
addresses bursts from zero and the conversion happens here, so that a burst has the same number it has
everywhere else in this package. `url` is the extractor's address of either of the burst's files, as
ASF's search lists it: `https://sentinel1-burst.asf.alaska.edu/<slc>/IW<n>/<POL>/<i>.tiff`.

The raster is downloaded into `dir` when it is first read, a temporary directory by default; pass one
to keep rasters across sessions, and a raster already there is not fetched again. The metadata is held
in memory.

Authentication is from `~/.netrc`, so a `urs.earthdata.nasa.gov` entry is the only setup. The extractor
generates a burst on first request and answers `202` until it is ready; this asks again every
`ASF_BURST_WAIT` seconds, up to `ASF_BURST_ATTEMPTS` times.

# Examples

```julia
b = AsfBurst("S1A_IW_SLC__1SSH_20151120T080202_...", 2, "HH", 3)
b = AsfBurst("https://sentinel1-burst.asf.alaska.edu/S1A_IW_SLC__1SSH_20151120T080202_.../IW2/HH/2.tiff")
s = open_slc(b; orbit = eof)
```
"""
struct AsfBurst <: AbstractSLCSource
    slc::String
    swath::Int
    polarization::String
    burst::Int
    dir::Union{Nothing,String}
end

function AsfBurst(slc::AbstractString, swath::Integer, polarization::AbstractString,
                  burst::Integer; dir = nothing)
    return AsfBurst(String(slc), Int(swath), _check_burst(swath, burst, polarization), Int(burst),
                    dir === nothing ? nothing : String(dir))
end

# A burst's subswath and index checked, and its polarization in upper case.
function _check_burst(swath::Integer, burst::Integer, polarization::AbstractString)
    swath in 1:3 || throw(ArgumentError(
        "Sentinel-1 IW subswaths are 1, 2 and 3; got $swath"))
    burst >= 1 || throw(ArgumentError(
        "bursts are numbered from 1, so burst $burst does not exist"))
    pol = uppercase(String(polarization))
    pol in ("HH", "HV", "VH", "VV") || throw(ArgumentError(
        "a Sentinel-1 polarization is HH, HV, VH or VV; got `$polarization`"))
    return pol
end

const _ASF_URL = r"^https://sentinel1-burst\.asf\.alaska\.edu/(S1[A-D]_IW_SLC__[0-9A-Z_]+)/IW([1-3])/(HH|HV|VH|VV)/(\d+)\.(tiff|xml)$"

"""
    is_asf_burst_url(url) -> Bool

Whether `url` is the burst extractor's address of a burst's raster or metadata.
"""
is_asf_burst_url(url::AbstractString) = occursin(_ASF_URL, url)

function AsfBurst(url::AbstractString; dir = nothing)
    m = match(_ASF_URL, url)
    m === nothing && throw(ArgumentError(
        "`$url` is not a burst extractor URL; one reads \
         `https://sentinel1-burst.asf.alaska.edu/<SLC granule>/IW<1-3>/<HH|HV|VH|VV>/<index>.tiff`"))
    return AsfBurst(something(m.captures[1]), parse(Int, something(m.captures[2])), something(m.captures[3]),
                    parse(Int, something(m.captures[4])) + 1; dir)
end

# Built a few pieces at a time: a concatenation of many arguments is compiled for an argument tuple of
# unknown length, which a trimmed executable cannot dispatch.
"""
    asf_burst_url(b::AsfBurst, extension) -> String

The URL of one of a burst's two files, `"tiff"` or `"xml"`.
"""
function asf_burst_url(b::AsfBurst, extension::AbstractString)
    dir = "https://" * ASF_BURST_HOST * "/" * b.slc
    file = string(b.burst - 1) * "." * extension
    return dir * "/IW" * string(b.swath) * "/" * b.polarization * "/" * file
end

# Where a burst's raster is kept, however it is fetched: named after its product, subswath, polarization
# and zero-based index.
function _burst_local(b)
    file = lowercase(b.polarization) * "-" * string(b.burst - 1) * ".tiff"
    return joinpath(b.dir === nothing ? _scratch_dir() : b.dir, b.slc * "-iw" * string(b.swath) * "-" * file)
end

# One directory per session for rasters downloaded or inflated with no directory named.
const _SCRATCH_DIR = Ref{Union{Nothing,String}}(nothing)
function _scratch_dir()
    d = _SCRATCH_DIR[]
    d === nothing || return d
    d = mktempdir(; cleanup = true)
    _SCRATCH_DIR[] = d
    return d
end

localpath(b::AsfBurst) = only(asf_rasters([b]))

"""
    with_claimed(fetch!, paths, what) -> paths

The files `paths`, each fetched by this call or by another already fetching it. `fetch!(todo)` fetches
`paths[todo]`, the ones neither local nor being fetched; a later caller asking for one of them waits for
it rather than fetching it again. `what` names the fetch in the error for a file that did not arrive.
"""
function with_claimed(fetch!::F, paths::Vector{String}, what::String) where {F}
    todo = Int[]
    waits = Base.Event[]
    @lock _DOWNLOADING_LOCK for k in eachindex(paths)
        _downloaded(paths[k]) && continue
        if haskey(_DOWNLOADING, paths[k])
            push!(waits, _DOWNLOADING[paths[k]])
        else
            _DOWNLOADING[paths[k]] = Base.Event()
            push!(todo, k)
        end
    end
    try
        isempty(todo) || fetch!(todo)
    finally
        @lock _DOWNLOADING_LOCK for k in todo
            notify(pop!(_DOWNLOADING, paths[k]))
        end
    end
    foreach(wait, waits)
    for p in paths
        _downloaded(p) || throw(ErrorException(what * " `" * p * "` failed in the task that started it"))
    end
    return paths
end

# How long to keep asking. The extractor answers `202` while it builds a burst, which takes a minute or
# two; a request that never succeeds should say so rather than retry forever.
const ASF_BURST_ATTEMPTS = 20
const ASF_BURST_WAIT = 15.0

"""
    asf_resolve(urls) -> Vector{String}

The pre-signed URL each extractor URL redirects to once its burst has been generated. All are asked at
once, so the extractor generates them together; those still answering `202` are asked again every
`ASF_BURST_WAIT` seconds, up to `ASF_BURST_ATTEMPTS` rounds.
"""
function asf_resolve(urls::Vector{String})
    out = fill("", length(urls))
    pending = collect(eachindex(urls))
    for attempt in 1:ASF_BURST_ATTEMPTS
        again = Int[]
        answers = http_resolve(urls[pending])
        for (k, (status, final)) in zip(pending, answers)
            if status == 202
                push!(again, k)
            elseif status == 200 || status == 206
                # A login page served with 200 is where the chain ends when Earthdata Login did not
                # authenticate the request.
                occursin(EDL_HOST, final) && _asf_unauthenticated(urls[k])
                out[k] = final
            else
                _asf_refused(urls[k], status)
            end
        end
        isempty(again) && return out
        pending = again
        attempt < ASF_BURST_ATTEMPTS && sleep(ASF_BURST_WAIT)
    end
    _asf_not_ready(urls[first(pending)], length(pending) - 1)
end

@noinline _asf_unauthenticated(url::String) = throw(ArgumentError(
    "`$url` ended at Earthdata Login rather than at the burst; ~/.netrc needs a working \
     `machine urs.earthdata.nasa.gov` entry"))

@noinline _asf_refused(url::String, status::Int) = throw(ArgumentError(
    "`$url` returned HTTP $status. A 401 means ~/.netrc has no working `machine urs.earthdata.nasa.gov` \
     entry; a 404 means the granule, subswath, polarization or burst index does not exist"))

@noinline _asf_not_ready(url::String, others::Int) = throw(ArgumentError(
    "`$url` and $others other burst files were still being generated after \
     $(round(Int, ASF_BURST_ATTEMPTS * ASF_BURST_WAIT)) s. ASF builds a burst on first request; try again \
     once it has had longer"))

# The metadata file of each product and polarization, as fetched.
const _ASF_METADATA = Dict{Tuple{String,String},String}()
const _ASF_METADATA_LOCK = ReentrantLock()

_asf_key(b::AsfBurst) = (b.slc, b.polarization)

"""
    asf_metadata(bursts) -> Vector{String}

The metadata file of each burst's parent product and polarization, fetched once per pair in the
session and shared by every burst of it.
"""
function asf_metadata(bursts::AbstractVector{AsfBurst})
    @lock _ASF_METADATA_LOCK begin
        missing_ = AsfBurst[]
        for b in bursts
            haskey(_ASF_METADATA, _asf_key(b)) || any(m -> _asf_key(m) == _asf_key(b), missing_) ||
                push!(missing_, b)
        end
        if !isempty(missing_)
            signed = asf_resolve([asf_burst_url(b, "xml") for b in missing_])
            t = HTTPTransport()
            for (b, url) in zip(missing_, signed)
                _ASF_METADATA[_asf_key(b)] = String(fetchall(t, url))
            end
        end
        return String[_ASF_METADATA[_asf_key(b)] for b in bursts]
    end
end

"""
    asf_rasters(bursts) -> Vector{String}

The local path of each burst's raster, downloading those not already in their directories. The
downloads are resolved together, then read 16 MiB ranges at a time, several ranges in flight.
"""
function asf_rasters(bursts::AbstractVector{AsfBurst})
    paths = [_burst_local(b) for b in bursts]
    return with_claimed(paths, "downloading") do todo
        signed = asf_resolve([asf_burst_url(bursts[k], "tiff") for k in todo])
        for (k, url) in zip(todo, signed)
            mkpath(dirname(paths[k]))
            _download(url, paths[k])
        end
    end
end

# The rasters being downloaded, which a second caller waits for rather than downloads again.
const _DOWNLOADING = Dict{String,Base.Event}()
const _DOWNLOADING_LOCK = ReentrantLock()

_downloaded(path::String) = isfile(path) && filesize(path) > 0

# A download is fetched a batch of ranges at a time, so a raster is never held in memory whole.
const _BLOCK = 16 * 1024 * 1024
const _BLOCKS_IN_FLIGHT = 8

# `url` into `dest`, through a `.part` file so an interrupted download is never mistaken for a raster.
function _download(url::String, dest::String)
    t = HTTPTransport()
    n = Int(objectsize(t, url))
    part = dest * ".part"
    io = open(part, "w")
    try
        for start in 0:(_BLOCK * _BLOCKS_IN_FLIGHT):(n - 1)
            stop = min(n, start + _BLOCK * _BLOCKS_IN_FLIGHT)
            ranges = [ByteRange(UInt64(o), UInt64(min(_BLOCK, stop - o))) for o in start:_BLOCK:(stop - 1)]
            for block in fetchranges(t, url, ranges)
                write(io, block)
            end
        end
    finally
        close(io)
    end
    filesize(part) == n || throw(ErrorException("downloading `$dest` wrote $(filesize(part)) of its $n bytes"))
    mv(part, dest; force = true)
    return dest
end
