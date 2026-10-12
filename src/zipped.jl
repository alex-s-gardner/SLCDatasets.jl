# A zipped Sentinel-1 product — a local zip, or one over HTTP(S) or in S3 — through ChunkManifests.jl's zip
# reader: the central directory and the annotation members are read in place, and a measurement raster,
# which the zip stores deflated, is inflated whole into a directory and read from there. Inflating a
# member once is what any read of it costs: a deflate stream is not addressable without inflating what
# precedes it.

"""
    ZipArchive

A zip archive's location, the location its bytes are read from, the transport they are read through, and
its central directory.
"""
struct ZipArchive
    uri::String
    source::String
    transport::LeafTransport
    members::Vector{ZipMember}
end

# Each remote archive opened this session, with when, so reopening one — a product's annotations, then each
# of its subswaths' rasters — reads its central directory once. Its pre-signed URL is resolved again once it
# is `_ZIP_FRESH` seconds old, well inside the hour such a URL is valid for. A local archive is read again,
# which costs nothing and sees the file as it is.
const _ZIPS = Dict{String,Tuple{Float64,ZipArchive}}()
const _ZIPS_LOCK = ReentrantLock()
const _ZIP_FRESH = 600.0

function ZipArchive(uri::AbstractString)
    u = String(uri)
    _is_remote(u) || return (t = LocalTransport(); ZipArchive(u, u, t, zip_members(u; transport = t)))
    @lock _ZIPS_LOCK begin
        hit = get(_ZIPS, u, nothing)
        hit === nothing || time() - hit[1] > _ZIP_FRESH || return hit[2]
    end
    source, t = _zip_source(u)
    r = ZipArchive(u, source, t, zip_members(source; transport = t))
    @lock _ZIPS_LOCK _ZIPS[u] = (time(), r)
    return r
end

# Where a product location's bytes are read from, and through what. An HTTP(S) location is followed once to
# where it redirects — through Earthdata Login with the `~/.netrc` credentials, for an Earthdata host such
# as ASF's, which answers a bearer token for a zip with `401` — so each range request after goes straight
# to the pre-signed URL rather than through the login chain.
function _zip_source(uri::String)
    startswith(uri, "s3://") && return (uri, S3Transport(first(split(uri[6:end], '/'; limit = 2))))
    (startswith(uri, "http://") || startswith(uri, "https://")) || return (uri, LocalTransport())
    status, final = only(http_resolve([uri]))
    (status == 200 || status == 206) || throw(ArgumentError(
        "`$uri` returned HTTP $status; an Earthdata host needs a `machine urs.earthdata.nasa.gov` entry in \
         ~/.netrc"))
    occursin(EDL_HOST, final) && throw(ArgumentError(
        "`$uri` ended at Earthdata Login rather than at the archive; ~/.netrc needs a working \
         `machine urs.earthdata.nasa.gov` entry"))
    return (final, HTTPTransport())
end

_is_remote(uri::AbstractString) = occursin(r"^[a-zA-Z][a-zA-Z0-9+.\-]*://", uri)

with_zip(f, path::AbstractString) = f(ZipArchive(path))

zip_names(r::ZipArchive) = [m.name for m in r.members]

function _member(r::ZipArchive, name::AbstractString)
    for m in r.members
        m.name == name && return m
    end
    throw(ArgumentError("`" * r.uri * "` has no member " * String(name)))
end

read_zip_entry(r::ZipArchive, name::AbstractString) =
    String(read_member(r.source, _member(r, name); transport = r.transport))

# Each measurement member inflated so far, by archive and member name, and where it was written; and those
# being inflated, which a second caller waits for rather than inflates again.
const _EXTRACTED = Dict{Tuple{String,String},String}()
const _EXTRACTING = Dict{Tuple{String,String},Base.Event}()
const _EXTRACTED_LOCK = ReentrantLock()

"""
    extract_measurements(r::ZipArchive, names, dir) -> Vector{String}

The local file each of the members `names` of `r` is inflated into under `dir`, inflating those not yet
inflated in this session concurrently: a remote archive is then read through several streams at once, and
each member's compressed bytes are transferred exactly once. Two archives are inflated at once when two
tasks ask.
"""
function extract_measurements(r::ZipArchive, names::Vector{String}, dir::String)
    want = [(r.uri, n) for n in names]
    mine = Tuple{String,String}[]
    waits = Base.Event[]
    @lock _EXTRACTED_LOCK for k in want
        haskey(_EXTRACTED, k) && continue
        if haskey(_EXTRACTING, k)
            push!(waits, _EXTRACTING[k])
        else
            _EXTRACTING[k] = Base.Event()
            push!(mine, k)
        end
    end
    if !isempty(mine)
        mkpath(dir)
        tasks = [Threads.@spawn extract_member(r.source, _member(r, k[2]), joinpath(dir, basename(k[2]));
                                               transport = r.transport) for k in mine]
        try
            for (k, t) in zip(mine, tasks)
                file = fetch(t)::String
                @lock _EXTRACTED_LOCK _EXTRACTED[k] = file
            end
        finally
            @lock _EXTRACTED_LOCK for k in mine
                notify(pop!(_EXTRACTING, k))
            end
        end
    end
    foreach(wait, waits)
    return @lock _EXTRACTED_LOCK begin
        for k in want
            haskey(_EXTRACTED, k) || throw(ErrorException(
                "inflating " * k[2] * " of " * k[1] * " failed in the task that started it"))
        end
        String[_EXTRACTED[k] for k in want]
    end
end
