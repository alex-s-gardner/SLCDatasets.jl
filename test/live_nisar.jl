# Reading a real granule over the network.
#
# Needs an Earthdata Login token, from `EARTHDATA_TOKEN` or the `~/.netrc` entry for
# `urs.earthdata.nasa.gov`, and transfers a few megabytes, so it runs only when
# `SAR_LIVE_TEST=1`. What it proves that the fixture tests cannot: that a multi-gigabyte product on a
# DAAC opens through the redirect chain, read as `open_slc` reads it — through `RemoteHTTP`'s transport,
# which follows the chain once — and that a window of it reads by range.
#
# The values asserted are the committed fixture's, which was harvested from this granule with `h5py`.
# So a pass here means the Julia reader and the reference agree on a real product, bitwise.

using SLCDatasets
using Test

const LIVE_URL = get(ENV, "SAR_LIVE_URL",
    "https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_L1_RSLC_BETA_V1/" *
    "NISAR_L1_PR_RSLC_003_170_D_053_7700_SHNA_A_20251028T235201_20251028T235238_X05009_N_P_J_001/" *
    "NISAR_L1_PR_RSLC_003_170_D_053_7700_SHNA_A_20251028T235201_20251028T235238_X05009_N_P_J_001.h5")

let opened = @elapsed(s = SLCDatasets._open_nisar(RemoteHTTP(LIVE_URL); frequency = nothing))
    o = orbit(s)

    @testset "a window reads by range" begin
        read = @elapsed(w = pixels(s)[20001:20512, 20001:20512])
        @test size(w) == (512, 512)
        @test any(!iszero, w)
        @info "live granule" opened_seconds = round(opened; digits = 1) window_seconds = round(read; digits = 1)
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
