# Reading a remote product in place.
#
# A NISAR RSLC is tens of gigabytes, essentially all of it image samples, and its metadata is tens of
# kilobytes within the first few megabytes. ChunkManifests.jl scans the product with byte-range requests
# for that metadata alone, and the samples are then read a chunk at a time as they are indexed, so neither
# opening a product nor reading a window of it transfers the rest.

"""
    RemoteHTTP(url; token = get(ENV, "EARTHDATA_TOKEN", nothing))

A product reachable over HTTP(S) from a server honoring range requests, read in place.

`token` is an Earthdata Login bearer token, sent as `Authorization: Bearer <token>`. An Earthdata host
answers it with a redirect to a pre-signed URL, and the header is not sent on to that other host. A
public URL needs no token; without one, no `Authorization` header is sent.
"""
struct RemoteHTTP <: AbstractSLCSource
    url::String
    token::Union{Nothing,String}
end

function RemoteHTTP(url::AbstractString; token = get(ENV, "EARTHDATA_TOKEN", nothing))
    startswith(url, "http://") || startswith(url, "https://") || throw(ArgumentError(
        "an HTTP source must be addressed as `https://host/path`, got `$url`"))
    token === nothing || !isempty(token) || throw(ArgumentError(
        "`token` is empty; pass a token, or `nothing` for a public URL"))
    return RemoteHTTP(String(url), token === nothing ? nothing : String(token))
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

source_transport(src::RemoteHTTP) = src.token === nothing ? nothing :
    HTTPTransport(; headers = ["Authorization" => "Bearer $(src.token)"])

# Only what is read in place reaches a remote product; a Sentinel-1 raster needs an unpacked `.SAFE`.
localpath(src::Union{RemoteHTTP,RemoteS3}) = throw(ArgumentError(
    "`$(source_uri(src))` is read in place and has no local path. A remote product is opened with " *
    "`open_slc` as a NISAR product; a Sentinel-1 product has to be downloaded and unpacked first"))
