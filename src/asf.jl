# Bursts as ASF distributes them, one file per burst.
#
# ASF's burst extractor serves a Sentinel-1 SLC as its individual bursts: `{i}.tiff` is one burst's
# raster and `{i}.xml` its metadata, both indexed from zero at
# `sentinel1-burst.asf.alaska.edu/{slc}/{swath}/{pol}/`. So a caller wanting a few bursts of a subswath
# transfers those bursts rather than the several gigabytes of the whole product.
#
# The raster has the same layout a `.SAFE` measurement file does — uncompressed, one strip per line,
# complex 16-bit samples — so `open_measurement` reads it unchanged. Unlike a `.SAFE`, each is its own file
# holding one burst, which is why a merge places them at row offset zero rather than at a burst's offset
# into a stacked raster.
#
# The XML is not a slice of the annotation. It wraps the *whole* product's annotation for every subswath:
# `<burst><metadata><product><swath>IW2</swath><content>` holds a complete product annotation listing
# every burst of the parent SLC, not just the one fetched. So `read_annotation` parses the `<content>`
# subtree as it stands, and which burst the file holds is known from its index rather than from the XML.

"""
    ASF_BURST_HOST

The host ASF's burst extractor serves from.
"""
const ASF_BURST_HOST = "sentinel1-burst.asf.alaska.edu"

"""
    AsfBurst(slc, swath, polarization, burst; dir = nothing) -> AsfBurst

One burst of a Sentinel-1 SLC, from ASF's burst extractor.

`slc` is the parent granule's name, `burst` its **1-based** index within the subswath — the extractor
addresses bursts from zero and the conversion happens here, so that a burst has the same number it has
everywhere else in this package.

Both files are fetched whole rather than by range: a burst raster is a couple of hundred megabytes and is
wanted entire. `dir` is where they are kept, defaulting to a temporary directory; pass one to keep them
across sessions, and an already-downloaded file there is not fetched again.

Authentication is from `~/.netrc`, so a `urs.earthdata.nasa.gov` entry is the only setup. The extractor generates a burst on first request and answers `202` until it is ready;
this retries until it is.

# Examples

```julia
b = AsfBurst("S1A_IW_SLC__1SSH_20151120T080202_...", 2, "HH", 3)
s = open_slc(b)
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
    swath in S1_SUBSWATHS || throw(ArgumentError(
        "Sentinel-1 IW subswaths are 1, 2 and 3; got $swath"))
    burst >= 1 || throw(ArgumentError(
        "bursts are numbered from 1, so burst $burst does not exist"))
    pol = uppercase(String(polarization))
    pol in ("HH", "HV", "VH", "VV") || throw(ArgumentError(
        "a Sentinel-1 polarization is HH, HV, VH or VV; got `$polarization`"))
    return AsfBurst(String(slc), Int(swath), pol, Int(burst),
                    dir === nothing ? nothing : String(dir))
end

# The extractor numbers bursts from zero.
_asf_stem(b::AsfBurst) = string(ASF_BURST_HOST, "/", b.slc, "/IW", b.swath, "/", b.polarization,
                                "/", b.burst - 1)

"""
    asf_burst_url(b::AsfBurst, extension) -> String

The URL of one of a burst's two files, `"tiff"` or `"xml"`.
"""
asf_burst_url(b::AsfBurst, extension::AbstractString) =
    string("https://", _asf_stem(b), ".", extension)

_asf_local(b::AsfBurst, extension::AbstractString) =
    joinpath(b.dir === nothing ? _default_asf_dir() : b.dir,
             string(b.slc, "-iw", b.swath, "-", lowercase(b.polarization), "-",
                    b.burst - 1, ".", extension))

# One directory per session, so the two files of a burst and the bursts of a subswath share it.
const _ASF_DIR = Ref{Union{Nothing,String}}(nothing)
function _default_asf_dir()
    d = _ASF_DIR[]
    d === nothing || return d
    d = mktempdir(; cleanup = true)
    _ASF_DIR[] = d
    return d
end

localpath(b::AsfBurst) = _asf_fetch(b, "tiff")

"""
    asf_annotation_path(b::AsfBurst) -> String

The local path of a burst's metadata file, fetching it if it is not already there.
"""
asf_annotation_path(b::AsfBurst) = _asf_fetch(b, "xml")

# How long to keep asking. The extractor answers `202` while it builds a product, and a burst raster
# takes a couple of minutes; a request that never succeeds should say so rather than retry forever.
const ASF_BURST_ATTEMPTS = 20
const ASF_BURST_WAIT = 15.0

function _asf_fetch(b::AsfBurst, extension::AbstractString)
    dest = _asf_local(b, extension)
    isfile(dest) && filesize(dest) > 0 && return dest
    mkpath(dirname(dest))
    url = asf_burst_url(b, extension)

    partial = string(dest, ".part")
    for attempt in 1:ASF_BURST_ATTEMPTS
        code = _asf_curl(url, partial)
        if code == 200
            # A `202` body is a few bytes of JSON, and so is an error page; a burst file is not.
            filesize(partial) > 0 || throw(ArgumentError(
                "`$url` returned an empty file"))
            mv(partial, dest; force = true)
            return dest
        elseif code == 202
            # Still being generated. Nothing is wrong yet, so wait and ask again.
            rm(partial; force = true)
            attempt == ASF_BURST_ATTEMPTS && break
            sleep(ASF_BURST_WAIT)
        else
            rm(partial; force = true)
            throw(ArgumentError(
                "`$url` returned HTTP $code. A `401` means `~/.netrc` has no working " *
                "`machine urs.earthdata.nasa.gov` entry; a `404` means the granule, subswath, " *
                "polarization or burst index does not exist"))
        end
    end
    rm(partial; force = true)
    throw(ArgumentError(
        "`$url` was still being generated after " *
        "$(round(Int, ASF_BURST_ATTEMPTS * ASF_BURST_WAIT)) s. ASF builds a burst on first request; " *
        "try again once it has had longer"))
end

# The cookie jar is what makes this work: the request is redirected to Earthdata to authenticate and
# back, and without somewhere to keep the session cookie the redirect lands on the login service and
# returns a page of HTML with status 200.
function _asf_curl(url::AbstractString, dest::AbstractString)
    jar = joinpath(dirname(dest), ".asf-cookies")
    out = IOBuffer()
    cmd = `curl --silent --show-error --location --netrc
                --cookie-jar $jar --cookie $jar
                --max-time 900 --output $dest --write-out "%{http_code}" $url`
    try
        run(pipeline(cmd; stdout = out, stderr = devnull))
    catch e
        throw(ArgumentError("fetching `$url` failed: $(sprint(showerror, e))"))
    end
    text = strip(String(take!(out)))
    code = tryparse(Int, text)
    code === nothing && throw(ArgumentError(
        "fetching `$url` returned no HTTP status; curl said `$text`"))
    return code
end

"""
    asf_annotation(b::AsfBurst) -> Node

The product annotation of the subswath a burst belongs to, from the burst's metadata file.

The file wraps a full annotation per subswath and polarization of the parent granule, each under a
`<content>` element whose children are the sections a `.SAFE` annotation has at its root. So that
element *is* the annotation, and what is returned lists every burst of the subswath rather than only
this one.
"""
function asf_annotation(b::AsfBurst)
    # The document is held by the returned node, so it must stay reachable: `readxml`'s result is kept
    # alive by the node's own reference to its document.
    r = root(readxml(asf_annotation_path(b)))
    want = string("IW", b.swath)
    for node in findall("metadata/product", r)
        _findtext(node, "swath") == want || continue
        uppercase(_findtext(node, "polarisation")) == b.polarization || continue
        content = findfirst("content", node)
        content === nothing && break
        return content
    end
    throw(ArgumentError(
        "`$(asf_annotation_path(b))` carries no product annotation for subswath IW$(b.swath) " *
        "polarization $(b.polarization)"))
end

"""
    AsfBurstBackend <: AbstractSLCBackend

One ASF-delivered burst: its own raster, and the annotation of the subswath it belongs to.
"""
struct AsfBurstBackend <: AbstractBurstBackend
    source::AsfBurst
    annotation::SubswathAnnotation
    orbit::String
end

_path(b::AsfBurstBackend) = asf_burst_url(b.source, "tiff")
_leading_annotation(b::AsfBurstBackend) = b.annotation

burst_index(b::AsfBurstBackend) = b.source.burst
burst_swath(b::AsfBurstBackend) = b.source.swath
burst_polarization(b::AsfBurstBackend) = lowercase(b.source.polarization)
burst_source(b::AsfBurstBackend) = b.source.slc
orbit_path(b::AsfBurstBackend) = b.orbit

function _asf_check_burst(b::AsfBurstBackend)
    n = nbursts(b.annotation)
    i = b.source.burst
    1 <= i <= n || throw(ArgumentError(
        "subswath IW$(b.source.swath) of `$(b.source.slc)` has $n bursts, so burst $i does not exist"))
    return i
end

_anchor(b::AsfBurstBackend) = b.annotation.burst_start[_asf_check_burst(b)]

# A burst's geometry and identification are its own, exactly as for a burst read from a `.SAFE`: the
# annotation an ASF burst carries is the whole subswath's, so the same fields describe it and the same
# helpers build the records.
read_geometry(b::AsfBurstBackend) =
    _s1_geometry(b.annotation, _anchor(b), b.annotation.lines_per_burst,
                 b.annotation.samples_per_burst)

read_identification(b::AsfBurstBackend) =
    _s1_identification(b.annotation, _anchor(b), _burst_stop(b.annotation, _anchor(b)))

function read_orbit(b::AsfBurstBackend)
    start = _anchor(b)
    return _s1_orbit(b.orbit, start, _burst_stop(b.annotation, start),
                     "burst $(b.source.burst) of subswath IW$(b.source.swath)")
end

# The file holds this burst alone, so its rows are the burst's own.
function read_pixels(b::AsfBurstBackend)
    a = b.annotation
    file = localpath(b.source)
    raster = open_measurement(file)
    size(raster, 1) == a.lines_per_burst || throw(ArgumentError(
        "`$file` has $(size(raster, 1)) lines but subswath IW$(b.source.swath) records " *
        "$(a.lines_per_burst) per burst"))
    size(raster, 2) == a.samples_per_burst || throw(ArgumentError(
        "`$file` is $(size(raster, 2)) samples wide but subswath IW$(b.source.swath) " *
        "records $(a.samples_per_burst)"))
    return raster
end

"""
    open_slc(b::AsfBurst; orbit) -> SLC

One ASF-delivered burst as an acquisition.

`orbit` is a POEORB or RESORB `.EOF` file, required for the same reason it is for a `.SAFE`: a
Sentinel-1 product carries no state vectors.
"""
function open_slc(b::AsfBurst; orbit = nothing)
    orbit_path = _require_orbit(orbit, string("burst ", b.burst, " of `", b.slc, "`"))
    a = read_annotation(asf_annotation(b), b.swath)
    return SLC(AsfBurstBackend(b, a, orbit_path))
end

"""
    asf_bursts(slc, swath, polarization, bursts; orbit, dir = nothing) -> Vector{SLC}

Several bursts of one subswath from ASF's burst extractor.

`bursts` is a range or vector of 1-based burst indices. The metadata file is fetched and parsed once
rather than per burst, since each of them carries the whole subswath's annotation; the rasters are
fetched as they are read.

The result merges with [`merge_bursts`](@ref).

# Examples

```julia
b = asf_bursts("S1A_IW_SLC__1SSH_20151120T080202_...", 2, "HH", 3:7; orbit = eof)
ref = merge_bursts(b)
```
"""
function asf_bursts(slc::AbstractString, swath::Integer, polarization::AbstractString,
                    bursts; orbit = nothing, dir = nothing)
    idx = collect(Int, bursts)
    isempty(idx) && throw(ArgumentError("no bursts were named"))
    orbit_path = _require_orbit(orbit, string("bursts of `", slc, "`"))

    # One fetch and one parse: every burst's file holds the same subswath annotation.
    a = read_annotation(asf_annotation(AsfBurst(slc, swath, polarization, first(idx); dir)), swath)
    return [SLC(AsfBurstBackend(AsfBurst(slc, swath, polarization, i; dir), a, orbit_path))
            for i in idx]
end

# How a merge of these reads its samples is in `merge.jl`, where `MergedBurstBackend` is defined.
