# Reading a remote product in place.
#
# A NISAR RSLC is tens of gigabytes, essentially all of it image samples, and its metadata is tens of
# kilobytes within the first few megabytes. ChunkManifests.jl scans the product with byte-range requests
# for that metadata alone, and the samples are then read a chunk at a time as they are indexed, so neither
# opening a product nor reading a window of it transfers the rest.

"""
    RemoteHTTP(url; token = get(ENV, "EARTHDATA_TOKEN", nothing), netrc = ...)

A product reachable over HTTP(S) from a server honoring range requests, read in place.

`token` is an Earthdata Login bearer token, sent as `Authorization: Bearer <token>`. An Earthdata host
answers it with a redirect to a pre-signed URL, and the header is not sent on to that other host. The
redirect is followed once and the product read from the pre-signed URL, which is followed again when it
expires.

Without a token, an Earthdata host — one under `nasa.gov` or `alaska.edu` — is read with a token
obtained from the `urs.earthdata.nasa.gov` entry of `netrc`, `~/.netrc` (or `ENV["NETRC"]`) by default,
when that file has one; see [`earthdata_token`](@ref). Any other URL is read with no `Authorization`
header, and `netrc = nothing` turns the lookup off.
"""
struct RemoteHTTP <: AbstractSLCSource
    url::String
    token::Union{Nothing,String}
    # A `.netrc` holding Earthdata Login credentials, read for a token when the product is opened.
    netrc::Union{Nothing,String}
end

function RemoteHTTP(url::AbstractString; token = get(ENV, "EARTHDATA_TOKEN", nothing),
                    netrc = _default_netrc(url, token))
    startswith(url, "http://") || startswith(url, "https://") || throw(ArgumentError(
        "an HTTP source must be addressed as `https://host/path`, got `$url`"))
    token === nothing || !isempty(token) || throw(ArgumentError(
        "`token` is empty; pass a token, or `nothing` for a public URL"))
    return RemoteHTTP(String(url), token === nothing ? nothing : String(token),
                      netrc === nothing ? nothing : String(netrc))
end

# The hosts Earthdata Login guards: NASA's DAACs and ASF's. A token is a credential, so it is sent only
# where it can be meant for.
function _earthdata_host(url::AbstractString)
    m = match(r"^https?://([^/:?#]+)", url)
    m === nothing && return false
    host = lowercase(m[1])
    return any(d -> host == d || endswith(host, "." * d), ("nasa.gov", "alaska.edu"))
end

function _default_netrc(url::AbstractString, token)
    token === nothing && _earthdata_host(url) || return nothing
    path = get(ENV, "NETRC", joinpath(homedir(), ".netrc"))
    isfile(path) && _netrc_login(path, EDL_HOST) !== nothing || return nothing
    return path
end

"""
    RemoteS3(uri)

A product in S3, addressed as `s3://bucket/key`, read in place.

Reading needs `using AWSS3`, which ChunkManifests.jl's S3 transport is built on, and credentials as
the AWS tools resolve them.
"""
struct RemoteS3 <: AbstractSLCSource
    uri::String

    function RemoteS3(uri::AbstractString)
        startswith(uri, "s3://") || throw(ArgumentError(
            "an S3 source must be addressed as `s3://bucket/key`, got `$uri`"))
        return new(String(uri))
    end
end

source_uri(src::RemoteHTTP) = src.url
source_uri(src::RemoteS3) = src.uri

function source_transport(src::RemoteHTTP)
    token = src.token !== nothing ? src.token :
            src.netrc !== nothing ? earthdata_token(src.netrc) : nothing
    token === nothing && return nothing
    return HTTPTransport(; headers = ["Authorization" => "Bearer $token"], resolve = true)
end

# Only what is read in place reaches a remote product; a Sentinel-1 raster needs an unpacked `.SAFE`.
localpath(src::Union{RemoteHTTP,RemoteS3}) = throw(ArgumentError(
    "`$(source_uri(src))` is read in place and has no local path. A remote product is opened with " *
    "`open_slc` as a NISAR product; a Sentinel-1 product has to be downloaded and unpacked first"))

# ---------------------------------------------------------------------------
# Earthdata Login
# ---------------------------------------------------------------------------

const EDL_HOST = "urs.earthdata.nasa.gov"

"""
    EDL_TOKEN_URL

Earthdata Login's `find_or_create_token` endpoint: given a user's credentials it returns one of their
bearer tokens, creating one only if they hold none.
"""
const EDL_TOKEN_URL = "https://urs.earthdata.nasa.gov/api/users/find_or_create_token"

# One token per `.netrc` per session: Earthdata Login issues the same token until it expires, so
# asking again on every open would only add a round trip.
const EDL_TOKENS = Dict{Tuple{String,String},String}()
const EDL_TOKENS_LOCK = ReentrantLock()

# The `login` and `password` of `machine` in a `.netrc`, which may spread an entry across lines.
function _netrc_login(path::AbstractString, machine::AbstractString)
    words = split(read(path, String))
    i = 1
    while i < length(words)
        if words[i] == "machine" && words[i + 1] == machine
            fields = Dict{String,String}()
            j = i + 2
            while j < length(words) && words[j] ∉ ("machine", "default")
                fields[words[j]] = words[j + 1]
                j += 2
            end
            haskey(fields, "login") && haskey(fields, "password") || return nothing
            return (fields["login"], fields["password"])
        end
        i += 1
    end
    return nothing
end

"""
    earthdata_token(netrc = joinpath(homedir(), ".netrc"); url = EDL_TOKEN_URL) -> String

An Earthdata Login bearer token for the `urs.earthdata.nasa.gov` credentials in `netrc`, from
[`EDL_TOKEN_URL`](@ref). The user's existing token is returned where there is one; otherwise Earthdata
Login creates one in their account. Fetched once per session for each file.
"""
function earthdata_token(netrc::AbstractString = joinpath(homedir(), ".netrc"); url::AbstractString = EDL_TOKEN_URL)
    key = (String(netrc), String(url))
    return @lock EDL_TOKENS_LOCK get!(EDL_TOKENS, key) do
        creds = isfile(netrc) ? _netrc_login(netrc, EDL_HOST) : nothing
        creds === nothing && throw(ArgumentError(
            "`$netrc` has no `machine $EDL_HOST login … password …` entry, so no Earthdata Login " *
            "token can be obtained from it. Add one, or pass `token` or set EARTHDATA_TOKEN"))
        status, body = http_request(url; method = "POST",
                                    headers = ["Authorization: Basic " * _base64(creds[1] * ":" * creds[2])])
        status == 200 || throw(ErrorException(
            "Earthdata Login refused the credentials in `$netrc` for a token: HTTP $status"))
        m = match(r"\"access_token\"\s*:\s*\"([^\"]+)\"", String(body))
        m === nothing && throw(ErrorException(
            "Earthdata Login answered without an access token, so none could be obtained"))
        String(m[1])
    end
end

const _BASE64_ALPHABET = codeunits("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/")

# Standard base64 with padding (RFC 4648 §4), for the Basic credentials header.
function _base64(s::String)
    src = codeunits(s)
    out = UInt8[]
    sizehint!(out, 4 * cld(length(src), 3))
    for i in 1:3:length(src)
        n = length(src) - i + 1
        b = UInt32(src[i]) << 16 | (n > 1 ? UInt32(src[i + 1]) << 8 : 0x0) | (n > 2 ? UInt32(src[i + 2]) : 0x0)
        for k in 0:3
            push!(out, k <= min(n, 3) ? _BASE64_ALPHABET[((b >> (18 - 6k)) & 0x3f) + 1] : UInt8('='))
        end
    end
    return String(out)
end
