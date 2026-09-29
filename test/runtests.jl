using SLCDatasets
using Test
import DiskArrays

include("fixture.jl")
include("sentinel1_fixture.jl")
include("tiff_fixture.jl")

@testset verbose = true "SLCDatasets.jl" begin
    @time @testset "measurement raster" begin include("tiff.jl") end
    @time @testset "NISAR reader" begin include("nisar.jl") end
    @time @testset "Sentinel-1 reader" begin include("sentinel1.jl") end
    @time @testset "burst grid" begin include("burstgrid.jl") end
    @time @testset "merged bursts" begin include("merge.jl") end
    @time @testset "ASF bursts" begin include("asf.jl") end
    @time @testset "access layer" begin include("remote.jl") end
    @time @testset "pairing" begin include("pairing.jl") end
    @time @testset "TOPS mode" begin include("tops.jl") end
    # Transfers a few megabytes from a DAAC and needs Earthdata credentials in `~/.netrc`.
    if get(ENV, "SAR_LIVE_TEST", "") == "1"
        @time @testset "live NISAR granule" begin include("live_nisar.jl") end
    else
        @info "skipping the live granule test; set SAR_LIVE_TEST=1 to run it"
    end
end
