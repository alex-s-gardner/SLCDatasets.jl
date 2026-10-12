# Reading a real granule over the network.
#
# Needs an Earthdata Login token, from `EARTHDATA_TOKEN` or the `~/.netrc` entry for
# `urs.earthdata.nasa.gov`, and transfers a few megabytes, so it runs only when
# `SAR_LIVE_TEST=1`. What it proves that the fixture tests cannot: that a multi-gigabyte product on a
# DAAC opens through the redirect chain, and that opening it and reading a window of it transfer a small
# fraction of it.
#
# The values asserted are the committed fixture's, which was harvested from this granule with `h5py`.
# So a pass here means the Julia reader and the reference agree on a real product, bitwise.

using SLCDatasets
using SLCDatasets: AbstractSLCSource
import ChunkManifests
using ChunkManifests: HTTPTransport, ByteRange
using Test

const LIVE_URL = get(ENV, "SAR_LIVE_URL",
    "https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_L1_RSLC_BETA_V1/" *
    "NISAR_L1_PR_RSLC_003_170_D_053_7700_SHNA_A_20251028T235201_20251028T235238_X05009_N_P_J_001/" *
    "NISAR_L1_PR_RSLC_003_170_D_053_7700_SHNA_A_20251028T235201_20251028T235238_X05009_N_P_J_001.h5")

# `RemoteHTTP` read through a transport that counts the bytes it fetches.
mutable struct Counted{T} <: ChunkManifests.AbstractTransport
    inner::T
    bytes::Threads.Atomic{Int}
end
function ChunkManifests.fetchrange(t::Counted, uri::AbstractString, r::ByteRange)
    b = ChunkManifests.fetchrange(t.inner, uri, r)
    Threads.atomic_add!(t.bytes, length(b))
    return b
end
for f in (:concurrency, :maxblock, :maxgap)
    @eval ChunkManifests.$f(t::Counted) = ChunkManifests.$f(t.inner)
end
ChunkManifests.objectsize(t::Counted, uri::AbstractString) = ChunkManifests.objectsize(t.inner, uri)

struct CountedHTTP <: AbstractSLCSource
    url::String
    transport::Counted
end
SLCDatasets.source_uri(s::CountedHTTP) = s.url
SLCDatasets.source_transport(s::CountedHTTP) = s.transport

let token = get(() -> earthdata_token(), ENV, "EARTHDATA_TOKEN")
    t = Counted(HTTPTransport(; headers = ["Authorization" => "Bearer $token"]), Threads.Atomic{Int}(0))
    s = SLCDatasets._open_nisar(CountedHTTP(LIVE_URL, t); frequency = nothing)
    o = orbit(s)
    total = ChunkManifests.objectsize(t.inner, LIVE_URL)

    @testset "opening transferred a small fraction of the granule" begin
        @test total > 1_000_000_000            # a real granule, not a stub
        @test t.bytes[] < total ÷ 100
        @info "live granule opened" total transferred = t.bytes[]
    end

    @testset "a window reads by range" begin
        before = t.bytes[]
        w = pixels(s)[20001:20512, 20001:20512]
        @test size(w) == (512, 512)
        @test any(!iszero, w)
        @test t.bytes[] - before < total ÷ 100
    end

    @testset "the live read matches the committed fixture bitwise" begin
        g = s.geometry
        @test g.starting_range === gx(FIXTURE.geometry.starting_range)
        @test g.far_range === gx(FIXTURE.geometry.far_range)
        @test g.range_pixel_spacing === gx(FIXTURE.geometry.range_pixel_spacing)
        @test g.wavelength === gx(FIXTURE.geometry.wavelength)
        @test g.prf === gx(FIXTURE.geometry.prf)
        @test g.sensing_start === gx(FIXTURE.geometry.sensing_start)
        @test g.sensing_stop === gx(FIXTURE.geometry.sensing_stop)
        @test g.nlines == FIXTURE.geometry.nlines
        @test g.nsamples == FIXTURE.geometry.nsamples

        @test length(o.time) == FIXTURE.orbit.n
        @test all(o.time[i] === gx(FIXTURE.orbit.time[i]) for i in eachindex(o.time))
        @test all(o.position[i][c] === gx(FIXTURE.orbit.position[i][c])
                  for i in eachindex(o.position), c in 1:3)
        @test all(o.velocity[i][c] === gx(FIXTURE.orbit.velocity[i][c])
                  for i in eachindex(o.velocity), c in 1:3)

        @test s.identification.mission == FIXTURE.identification.mission
        @test s.identification.absolute_orbit == FIXTURE.identification.absolute_orbit
        @test s.identification.start_time == FIXTURE.identification.start_time
    end
end
