# The remote sources, exercised against local HTTP servers.
#
# An Earthdata host takes a bearer token and answers with a redirect to a pre-signed URL on another
# host, which refuses a request carrying an `Authorization` header. The two servers here play those two
# parts over the fixture product, so the whole path — scanning in place, the redirect, the header kept
# from the second host, and reading samples by range — is checked without a network. The live granule is
# `live_nisar.jl`.

using SLCDatasets
using SLCDatasets: localpath, source_for, source_uri, source_transport, _netrc_login, _earthdata_host
using Base64: base64encode
using ChunkManifests: HTTPTransport
using HTTP
using Test

@testset "source_for dispatches on the spelling" begin
    @test source_for("/data/x.h5") isa LocalFile
    @test source_for("https://example.org/x.h5") isa RemoteHTTP
    @test source_for("http://example.org/x.h5") isa RemoteHTTP
    @test source_for("s3://bucket/key.h5") isa RemoteS3
    @test localpath(LocalFile("/data/x.h5")) == "/data/x.h5"
    @test source_uri(RemoteS3("s3://bucket/key.h5")) == "s3://bucket/key.h5"
end

@testset "constructors reject what cannot work" begin
    @test_throws "must be addressed as `s3://bucket/key`" RemoteS3("https://example.org/x.h5")
    @test_throws "must be addressed as `https://host/path`" RemoteHTTP("s3://bucket/key.h5")
    @test_throws "`token` is empty" RemoteHTTP("https://example.org/x.h5"; token = "")
    withenv("EARTHDATA_TOKEN" => "from-env") do
        @test RemoteHTTP("https://example.org/x.h5").token == "from-env"
    end
    withenv("EARTHDATA_TOKEN" => nothing) do
        @test RemoteHTTP("https://example.org/x.h5").token === nothing
        @test source_transport(RemoteHTTP("https://example.org/x.h5")) === nothing
    end
    @test source_transport(RemoteHTTP("https://example.org/x.h5"; token = nothing, netrc = nothing)) === nothing
    @test source_transport(RemoteHTTP("https://example.org/x.h5"; token = "t")) isa HTTPTransport
    # Read from where the download redirects, not through the redirect chain on every request.
    @test source_transport(RemoteHTTP("https://example.org/x.h5"; token = "t")).redirects !== nothing
    @test_throws "is read in place and has no local path" localpath(RemoteHTTP("https://example.org/x.h5"))
end

# Serves `bytes` by range, refusing any request that carries an `Authorization` header, as a pre-signed
# URL does.
function _data_handler(bytes)
    return function (req)
        isempty(HTTP.header(req, "Authorization", "")) ||
            return HTTP.Response(400, "only one auth mechanism allowed")
        m = match(r"bytes=(\d+)-(\d+)", HTTP.header(req, "Range", ""))
        m === nothing && return HTTP.Response(200, bytes)
        a, b = parse(Int, m[1]), min(parse(Int, m[2]), length(bytes) - 1)
        return HTTP.Response(206, ["Content-Range" => "bytes $a-$b/$(length(bytes))"], bytes[(a + 1):(b + 1)])
    end
end

# Answers a request carrying the token with a redirect to `target`, and any other with 401.
_auth_handler(token, target) = req ->
    HTTP.header(req, "Authorization", "") == "Bearer $token" ?
    HTTP.Response(307, ["Location" => target]) : HTTP.Response(401, "unauthorized")

@testset "an Earthdata-style product reads in place through the redirect" begin
    A = ComplexF32[complex(i, -j) for i in 1:37, j in 1:23]
    mktempdir() do dir
        path = write_fixture_product(joinpath(dir, "remote.h5"); samples = A, chunk = (5, 8),
                                     filters = (:shuffle, :deflate))
        reference = open_slc(path)
        data = HTTP.serve!(_data_handler(read(path)), "127.0.0.1", 0)
        # Named as `localhost` so the redirect crosses to another host, as an Earthdata one does.
        target = "http://localhost:$(HTTP.port(data))/remote.h5"
        auth = HTTP.serve!(_auth_handler("secret", target), "127.0.0.1", 0)
        url = "http://127.0.0.1:$(HTTP.port(auth))/remote.h5"
        try
            s = open_slc(RemoteHTTP(url; token = "secret"))
            @test s.geometry.starting_range === reference.geometry.starting_range
            @test s.geometry.sensing_start === reference.geometry.sensing_start
            @test s.identification.absolute_orbit == reference.identification.absolute_orbit
            @test orbit(s).position == orbit(reference).position
            p = pixels(s)
            @test p[:, :] == A
            @test p[6:20, 3:18] == A[6:20, 3:18]

            @test_throws Exception open_slc(RemoteHTTP(url; token = "wrong"))
        finally
            close(auth)
            close(data)
        end
    end
end

@testset "a .netrc entry is read wherever it is spelled" begin
    mktempdir() do dir
        one = joinpath(dir, "one")
        write(one, "machine other.host login a password b\nmachine urs.earthdata.nasa.gov login me password s3cret\n")
        @test _netrc_login(one, "urs.earthdata.nasa.gov") == ("me", "s3cret")
        spread = joinpath(dir, "spread")
        write(spread, "machine urs.earthdata.nasa.gov\n    password p\n    login u\ndefault login x password y\n")
        @test _netrc_login(spread, "urs.earthdata.nasa.gov") == ("u", "p")
        @test _netrc_login(spread, "nowhere") === nothing
        @test_throws "has no `machine urs.earthdata.nasa.gov login" earthdata_token(spread * "-absent")
    end
end

@testset "only an Earthdata host is sent a token from .netrc" begin
    @test _earthdata_host("https://nisar.asf.earthdatacloud.nasa.gov/NISAR/x.h5")
    @test _earthdata_host("https://sentinel1.asf.alaska.edu/SLC/x.zip")
    @test !_earthdata_host("https://example.org/x.h5")
    @test !_earthdata_host("https://nasa.gov.example.org/x.h5")
    mktempdir() do dir
        netrc = joinpath(dir, "netrc")
        write(netrc, "machine urs.earthdata.nasa.gov login me password s3cret\n")
        withenv("EARTHDATA_TOKEN" => nothing, "NETRC" => netrc) do
            @test RemoteHTTP("https://data.asf.alaska.edu/x.h5").netrc == netrc
            @test RemoteHTTP("https://example.org/x.h5").netrc === nothing
            @test RemoteHTTP("https://data.asf.alaska.edu/x.h5"; token = "t").netrc === nothing
        end
    end
end

@testset "a token is asked of Earthdata Login with the .netrc credentials" begin
    seen = Ref("")
    edl = HTTP.serve!("127.0.0.1", 0) do req
        seen[] = HTTP.header(req, "Authorization", "")
        seen[] == "Basic " * base64encode("me:s3cret") ?
            HTTP.Response(200, "{\"access_token\":\"tok-1\",\"token_type\":\"Bearer\"}") :
            HTTP.Response(401, "bad credentials")
    end
    try
        mktempdir() do dir
            good = joinpath(dir, "good")
            write(good, "machine urs.earthdata.nasa.gov login me password s3cret\n")
            url = "http://127.0.0.1:$(HTTP.port(edl))/api/users/find_or_create_token"
            @test earthdata_token(good; url) == "tok-1"
            bad = joinpath(dir, "bad")
            write(bad, "machine urs.earthdata.nasa.gov login me password wrong\n")
            @test_throws "refused the credentials" earthdata_token(bad; url)
        end
    finally
        close(edl)
    end
end
