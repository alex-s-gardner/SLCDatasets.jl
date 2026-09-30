# Real Sentinel-1 IW geometry, rounded to values easy to check by hand: a 2.33 m range pixel, IW1
# near range 800000 m, each subswath 20000 samples wide with a 200-sample overlap into the next.
function _iw_annotation(swath::Int, mission::String, starting_range::Float64)
    t = SLCDatasets.UtcTime(DateTime(2020, 1, 1), 0.0)
    poly = [RangePolynomial(t, starting_range, (0.0, 0.0, 0.0))]
    return SLCDatasets.SubswathAnnotation(
        swath, "vv", 12345, mission, "SLC", "ascending",
        starting_range, 2.33,
        0.05546576, 0.002, 1500, 20000,
        [t], [0], [1499], [0], [19999],
        0.0, poly, poly)
end

function _iw123(mission::String)
    # IW1: [800000, 800000 + 19999*2.33] = [800000, 846597.67]
    # IW2 starts 200 samples before IW1's far range: 846597.67 - 199*2.33 = 846134.0
    # IW3 starts 200 samples before IW2's far range
    iw1 = _iw_annotation(1, mission, 800_000.0)
    iw1_far = iw1.starting_range + 19999 * iw1.range_pixel_spacing
    iw2 = _iw_annotation(2, mission, iw1_far - 199 * 2.33)
    iw2_far = iw2.starting_range + 19999 * iw2.range_pixel_spacing
    iw3 = _iw_annotation(3, mission, iw2_far - 199 * 2.33)
    return [iw1, iw2, iw3]
end

_product(mission::String) =
    SLCDatasets.Sentinel1Product("fake.SAFE", "fake.EOF", "vv", [1, 2, 3], _iw123(mission))

@testset "subswath_borders" begin
    ref = _product("S1A")

    @testset "same platform: not applied, but still computable" begin
        b = subswath_borders(ref, _product("S1A"))
        @test b.same_platform
        @test b.reference_platform == "S1A"
    end

    @testset "cross-platform border arithmetic" begin
        b = subswath_borders(ref, _product("S1B"))
        @test !b.same_platform
        @test b.reference_platform == "S1A"

        iw1, iw2, iw3 = _iw123("S1A")
        r0, sp0 = iw1.starting_range, iw1.range_pixel_spacing
        col(r) = round(Int, (r - r0) / sp0)
        far(a) = a.starting_range + 19999 * a.range_pixel_spacing

        @test b.ncols == col(far(iw3))
        @test b.border12 == (col(iw2.starting_range) + col(far(iw1))) / 2
        @test b.border23 == (col(iw3.starting_range) + col(far(iw2))) / 2
        # Every column index is >= 0 and increasing subswath-to-subswath, matching a real mosaic —
        # not a check on the arithmetic itself, which the two asserts above already pin down exactly.
        @test 0 < b.border12 < b.border23 < b.ncols
    end

    @testset "requires all three subswaths" begin
        two_swath = SLCDatasets.Sentinel1Product("fake.SAFE", "fake.EOF", "vv", [1, 2], _iw123("S1A")[1:2])
        @test_throws "all three subswaths" subswath_borders(two_swath, _product("S1B"))
        @test_throws "all three subswaths" subswath_borders(ref, two_swath)
    end
end
