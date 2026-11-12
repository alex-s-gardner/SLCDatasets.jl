# Whether an acquisition's samples carry a TOPS azimuth ramp.
#
# This is the question a consumer of `pixels` has to ask before interpolating them, and it lives here
# because the answer is a property of how the product was collected — which this package reads and a
# geometry or correlation package cannot infer.
#
# What makes it worth testing rather than obvious: the answer has to be right for four different backend
# kinds, and one of them does not inherit it. `MergedBurstBackend` describes several bursts and subtypes
# `AbstractSLCBackend` rather than `AbstractBurstBackend`, so it would default to `false` without a method
# of its own. A merge silently reported as non-TOPS is exactly the case that would reach a resampler.

using Dates: DateTime
using SLCDatasets
using SLCDatasets: AbstractBurstBackend, AbstractSLCBackend, MergedBurstBackend, NisarBackend,
                   Sentinel1Backend, is_tops, deramp_parameters, burst_at, DerampParameters,
                   RangePolynomial, nearest_polynomial
using Test

@testset "a NISAR acquisition is not TOPS" begin
    # Stripmap at zero Doppler, so the samples are directly interpolable and a consumer needs no
    # deramping. Read from the committed fixture, so this needs no granule.
    mktempdir() do dir
        s = open_slc(write_fixture_product(joinpath(dir, "a.h5")))
        @test !is_tops(s)
        @test !is_tops(s.backend)
        @test s.backend isa NisarBackend
    end
end

@testset "every Sentinel-1 form is TOPS" begin
    # A mosaic, a single burst, and a merge of bursts. All three describe TOPS data, and the three reach
    # the answer by different routes — the first two through `AbstractBurstBackend`, the third through its
    # own method — so all three are asserted rather than one standing for the family.
    if !isempty(S1_PRODUCTS)
        for p in S1_PRODUCTS
            pol = String(p.gold.polarization)
            mosaic = open_slc(p.safe; orbit = p.eof, swath = 2, polarization = pol)
            @test is_tops(mosaic)
            @test mosaic.backend isa Sentinel1Backend
            # The type hierarchy is what carries the answer for a burst, so that a backend added later
            # inherits it rather than having to remember to declare it.
            @test mosaic.backend isa AbstractBurstBackend

            b = bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol)
            @test is_tops(b[1])
            @test is_tops(b[end])
            @test all(is_tops, b)

            # The case that does not inherit: a merge subtypes `AbstractSLCBackend` directly.
            m = merge_bursts(b)
            @test m.backend isa MergedBurstBackend
            @test !(m.backend isa AbstractBurstBackend)
            @test is_tops(m)
            # A merge of one burst too, since that is the degenerate placement and takes a different path
            # through the grid.
            @test is_tops(merge_bursts(b[1:1]))
        end
    else
        @info "skipping the Sentinel-1 TOPS cases; no committed products resolved"
    end
end

@testset "a backend that does not say is not TOPS" begin
    # The default is `false` rather than an error, so that a non-TOPS sensor added later needs no negative
    # declaration. Checked against a backend type declared here, which the package has never seen — if
    # this needed a method the default would not be doing its job.
    struct _PlainBackend <: AbstractSLCBackend end
    @test !is_tops(_PlainBackend())
end

@testset "a non-TOPS acquisition has no ramp to remove" begin
    # Stripmap, so there is nothing to deramp. That is a different situation from "this reader cannot tell
    # you", and the message says which one it is rather than raising a `MethodError` a caller would read as
    # a missing feature.
    mktempdir() do dir
        s = open_slc(write_fixture_product(joinpath(dir, "a.h5")))
        err = try
            deramp_parameters(s)
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing
        @test occursin("not a TOPS acquisition", err)
        @test occursin("is_tops", err)
    end
end

@testset "a range polynomial evaluates about its reference range" begin
    # `r0` is a range, not an index: the argument is `range - r0`, so the constant term is the value there.
    # Getting this wrong would evaluate every polynomial at a range hundreds of kilometers off and is not
    # visible in a round trip, since remove and reapply would share the error.
    t = NanoDate(DateTime(2020, 1, 1, 0, 0, 0))
    p = RangePolynomial(t, 800_000.0, (2.0, 3.0, 4.0))
    @test p(800_000.0) == 2.0
    @test p(800_001.0) == 2.0 + 3.0 + 4.0
    @test p(799_999.0) == 2.0 - 3.0 + 4.0
    # Integer and Float32 ranges reach the same answer, since a caller's spacing arithmetic may produce
    # either.
    @test p(800_001) == p(800_001.0)
end

@testset "a burst selects the polynomial nearest its own time" begin
    # The annotation estimates these on its own schedule rather than once per burst, so selection is by
    # time. An off-by-one here would deramp with a neighbour's coefficients: plausible, and wrong in a way
    # that grows with how fast the estimates vary.
    base = DateTime(2020, 1, 1, 0, 0, 0)
    at(s) = SLCDatasets._advance(NanoDate(base), s)
    polys = [RangePolynomial(at(0.0), 800_000.0, (1.0, 0.0, 0.0)),
             RangePolynomial(at(10.0), 800_000.0, (2.0, 0.0, 0.0)),
             RangePolynomial(at(20.0), 800_000.0, (3.0, 0.0, 0.0))]

    @test nearest_polynomial(polys, at(0.0)).coeffs[1] == 1.0
    @test nearest_polynomial(polys, at(9.0)).coeffs[1] == 2.0
    @test nearest_polynomial(polys, at(11.0)).coeffs[1] == 2.0
    @test nearest_polynomial(polys, at(20.0)).coeffs[1] == 3.0
    # Before the first and after the last clamp to the ends rather than extrapolating.
    @test nearest_polynomial(polys, at(-100.0)).coeffs[1] == 1.0
    @test nearest_polynomial(polys, at(100.0)).coeffs[1] == 3.0
    # An exact tie takes the earlier entry, which is what a scan keeping the strict minimum does. Asserted
    # because it is a choice, not a consequence.
    @test nearest_polynomial(polys, at(5.0)).coeffs[1] == 1.0
    @test nearest_polynomial(polys, at(15.0)).coeffs[1] == 2.0

    @test_throws "lists no polynomials" nearest_polynomial(RangePolynomial[], at(0.0))
end

@testset "a product without the deramp fields is still readable" begin
    # The committed golden inputs carry none of the three, since they were not among the values dumped from
    # those granules. So the geometry and the amplitudes must still read, and only `deramp_parameters` may
    # object — a product this reader handles today must not become unreadable for want of a field nothing
    # else uses.
    if !isempty(S1_PRODUCTS)
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        b = bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol)
        @test nlines(b[1]) > 0
        @test is_tops(b[1])

        err = try
            deramp_parameters(b[1])
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing
        # Which field is missing, and where it lives, so the gap is actionable.
        @test occursin("azimuthFmRateList", err)
        @test occursin("must not be interpolated", err)
    else
        @info "skipping the missing-field case; no committed products resolved"
    end
end

@testset "a merged line resolves to the burst whose ramp it carries" begin
    # The ramp is quadratic about each burst's own center, so a merged image carries one per burst. A line
    # attributed to the wrong burst would be deramped with a phase referenced up to half a burst away.
    if !isempty(S1_PRODUCTS)
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        b = bursts(p.safe; orbit = p.eof, swath = 2, polarization = pol)
        m = merge_bursts(b)
        grid = m.backend.grid

        # Every placement's own rows resolve to it, at both ends and in the middle.
        for (k, pl) in enumerate(grid.placements)
            for line in (first(pl.grid_rows), last(pl.grid_rows),
                         (first(pl.grid_rows) + last(pl.grid_rows)) ÷ 2)
                burst, burst_line = burst_at(m, line)
                @test burst == k
                # The returned row is the burst's own, and the mapping is the placement's.
                @test burst_line == pl.burst_rows[line - first(pl.grid_rows) + 1]
            end
        end

        # Consecutive placements are disjoint and adjacent, so the seam between two bursts attributes each
        # side to its own — the property that makes the lookup unambiguous.
        for k in 2:length(grid.placements)
            @test first(burst_at(m, last(grid.placements[k - 1].grid_rows))) == k - 1
            @test first(burst_at(m, first(grid.placements[k].grid_rows))) == k
        end

        @test_throws "not one of them" burst_at(m, 0)
        @test_throws "not one of them" burst_at(m, grid.nlines + 1)

        # A single burst is all its own lines, so asking is most likely a confusion with the merged image
        # and gets told so rather than a plausible `(1, line)`.
        @test_throws "not a merge of bursts" burst_at(b[1], 1)
    else
        @info "skipping the merged burst_at cases; no committed products resolved"
    end
end

@testset "the deramp parameters of a synthetic product" begin
    # The three fields written into the fixture annotation, so the parse and the per-burst selection are
    # exercised end to end without needing a granule that carries them.
    #
    # Both coefficient spellings appear: newer IPF versions write one whitespace-separated text node, older
    # ones separate sibling elements. A reader handling only the current form would fail on an archive
    # product, and the failure would look like bad data rather than an unread format.
    if !isempty(S1_PRODUCTS)
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        inputs = deepcopy(p.inputs)
        for swath in ("1", "2", "3")
            sw = inputs["swaths"][swath]
            times = sw["burstAzimuthTimes"]
            sw["azimuthSteeringRate"] = "1.590368784"
            # One estimate per burst here, tagged at that burst's own recorded time, so the nearest-in-time
            # selection has an unambiguous answer to check against.
            sw["azimuthFmRateList"] = [Dict("azimuthTime" => t, "t0" => "0.005",
                                            "coefficients" => ["-2300.0", "0.5", "-1.0e-5"])
                                       for t in times]
            sw["dcEstimateList"] = [Dict("azimuthTime" => t, "t0" => "0.005",
                                         "spelling" => "old",
                                         "coefficients" => ["-40.0", "0.25", "-2.0e-6"])
                                    for t in times]
        end

        mktempdir() do dir
            safe, eof = write_s1_fixture(dir, inputs)
            b = bursts(safe; orbit = eof, swath = 2, polarization = pol)
            d = deramp_parameters(b[1])

            @test d isa DerampParameters
            @test d.azimuth_steering_rate ≈ 1.590368784
            # `t0` is a two-way time, so the reference range is `t0 * c / 2`.
            @test d.azimuth_fm_rate.r0 ≈ 0.005 * SLCDatasets.SPEED_OF_LIGHT / 2
            @test d.azimuth_fm_rate.coeffs == (-2300.0, 0.5, -1.0e-5)
            # The old spelling reaches the same coefficients as the new one.
            @test d.doppler_centroid.coeffs == (-40.0, 0.25, -2.0e-6)

            # The scales come from the same annotation the geometry uses, so a deramp and a solve agree.
            @test d.lines_per_burst == nlines(b[1])
            @test d.azimuth_time_interval > 0
            @test d.starting_range > 0

            # The mid-time is half a burst after its start, which is what the ramp is referenced to and what
            # an orbit is interpolated at.
            g = b[1].geometry
            @test SLCDatasets.seconds_between(g.epoch, d.burst_mid) - g.sensing_start ≈
                  (d.lines_per_burst - 1) * d.azimuth_time_interval / 2

            # A later burst selects its own estimate, not the first one. With one estimate per burst that is
            # a different object, which is the whole point of selecting by time.
            if length(b) >= 2
                d2 = deramp_parameters(b[2])
                @test d2.azimuth_fm_rate.time != d.azimuth_fm_rate.time
                @test d2.burst_mid > d.burst_mid
            end

            # A merge takes a burst index over its own range, and burst k of the merge is burst k of the
            # annotation slice it merged.
            m = merge_bursts(b)
            @test deramp_parameters(m, 1).burst_mid == d.burst_mid
            @test_throws "not one of them" deramp_parameters(m, length(b) + 1)
            # And a single burst has exactly one.
            @test_throws "only burst index" deramp_parameters(b[1], 2)

            # The mosaic has no single ramp: its subswaths sit at different slant ranges.
            mosaic = open_slc(safe; orbit = eof, swath = 2, polarization = pol)
            @test_throws "mosaic across subswaths" deramp_parameters(mosaic)
        end
    else
        @info "skipping the synthetic deramp cases; no committed products resolved"
    end
end
