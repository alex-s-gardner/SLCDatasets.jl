# The placement of a subswath's bursts on one uniform azimuth grid, against the arithmetic
# hyp3-autorift's `merge_bursts_in_swath` performs.
#
# `reference/burst_grid.json` holds what that function computes for the granules the golden values come
# from, dumped by `reference/dump_burst_grid.py`. Two of its conventions differ from this reader's and
# the tests assert the difference rather than absorbing it: the reference writes bursts in turn and lets
# a later one overwrite the odd row a midpoint split leaves shared, so its *effective* extents are what
# a reader placing each burst once must match; and it slices range samples with an exclusive end against
# an inclusive index, dropping the last valid sample of every burst, which this reader keeps.

using SLCDatasets: read_annotation, annotation_xml, burst_grid, nbursts, valid_lines, valid_samples,
                   BurstGrid, MAX_BURST_GRID_RESIDUAL, placement_at
using JSON
using Test

const BURST_GRID_GOLD = let f = joinpath(@__DIR__, "reference", "burst_grid.json")
    isfile(f) ? JSON.parse(read(f, String)) : nothing
end

# The annotation of each committed case, however the product is reached.
function burst_grid_cases()
    out = NamedTuple[]
    BURST_GRID_GOLD === nothing && return out
    for p in S1_PRODUCTS
        safe_name = String(p.gold.safe)
        haskey(BURST_GRID_GOLD, Symbol(safe_name)) || continue
        pol = String(p.gold.polarization)
        for swath in 1:3
            a = read_annotation(annotation_xml(p.safe, swath, pol), swath)
            push!(out, (; safe_name, swath,
                        annotation = a,
                        gold = BURST_GRID_GOLD[Symbol(safe_name)][Symbol(string(swath))]))
        end
    end
    return out
end

const BURST_GRID_CASES = burst_grid_cases()

if isempty(BURST_GRID_CASES)
    @info "no burst-grid golden values found, so those tests are skipped."
else
    @testset "placement matches hyp3-autorift's merge" begin
        for case in BURST_GRID_CASES
            g = burst_grid(case.annotation)
            gold = case.gold
            @testset "$(case.safe_name) IW$(case.swath)" begin
                @test g.nlines == gold.num_az_lines
                @test g.nsamples == gold.samples_per_burst
                @test length(g.placements) == length(gold.rows)

                # A merged subswath is shorter than the stack of its bursts, because the overlap
                # appears once rather than twice. That difference is the whole point of the grid.
                @test g.nlines < gold.lines_per_burst * length(gold.rows)

                for (p, r) in zip(g.placements, gold.rows)
                    @test p.burst == r.burst
                    # The reference counts rows from zero and ends its slices exclusively.
                    @test first(p.grid_rows) == r.merge_start + 1
                    @test last(p.grid_rows) == r.eff_merge_end
                    @test first(p.burst_rows) == r.burst_start + 1
                    @test last(p.burst_rows) == r.eff_burst_end

                    @test first(p.cols) == r.first_valid_sample + 1
                    # Inclusive, so one sample wider than the reference's exclusive slice takes.
                    @test last(p.cols) == r.last_valid_sample + 1
                    @test length(p.cols) ==
                          (r.last_valid_sample - r.first_valid_sample) + 1
                end
            end
        end
    end

    # The properties that make the merged image one a `RadarGeometry` can describe: every row comes
    # from exactly one burst, and the rows a burst gives equal the rows it takes.
    @testset "the grid tiles the merged image" begin
        for case in BURST_GRID_CASES
            g = burst_grid(case.annotation)
            a = case.annotation

            for p in g.placements
                @test length(p.grid_rows) == length(p.burst_rows)
                @test first(p.burst_rows) >= a.first_valid_line[p.burst]
                @test last(p.burst_rows) <= a.last_valid_line[p.burst]
                @test p.cols == valid_samples(a, p.burst)
                @test 1 <= first(p.grid_rows) <= last(p.grid_rows) <= g.nlines
            end

            # Contiguous and disjoint: consecutive placements meet exactly.
            for k in 2:length(g.placements)
                @test first(g.placements[k].grid_rows) ==
                      last(g.placements[k - 1].grid_rows) + 1
            end

            # The residual is what makes placing a burst at a whole row a copy rather than a
            # resampling. Every granule measured is orders of magnitude inside the bound.
            @test g.residual < MAX_BURST_GRID_RESIDUAL / 10

            # Row lookup must agree with the placements it searches.
            @test placement_at(g, 1) === nothing ||
                  first(g.placements[1].grid_rows) == 1
            for p in g.placements
                @test placement_at(g, first(p.grid_rows)) === p
                @test placement_at(g, last(p.grid_rows)) === p
            end
            @test placement_at(g, 0) === nothing
            @test placement_at(g, g.nlines + 1) === nothing
            # The margin before the first burst's valid region belongs to no burst.
            if first(g.placements[1].grid_rows) > 1
                @test placement_at(g, first(g.placements[1].grid_rows) - 1) === nothing
            end
        end
    end

    # The reason for the uniform grid: a row's azimuth time is `sensing_start + row / prf` throughout,
    # which is exactly what `RadarGeometry` claims. A stack of whole bursts breaks it at every seam.
    @testset "a row's azimuth time is uniform across the merged image" begin
        for case in BURST_GRID_CASES
            a = case.annotation
            g = burst_grid(a)
            anchor = a.burst_start[1]
            dt = a.azimuth_time_interval
            for p in g.placements
                # The grid row this burst's own first line lands on, from the placement alone.
                row = first(p.grid_rows) - (first(p.burst_rows) - 1)
                # ...against the row its recorded time puts it on.
                from_time = seconds_between(anchor, a.burst_start[p.burst]) / dt
                @test abs((row - 1) - from_time) < MAX_BURST_GRID_RESIDUAL
            end
        end
    end

    @testset "a subset of a subswath's bursts merges" begin
        case = first(BURST_GRID_CASES)
        a = case.annotation
        n = nbursts(a)
        whole = burst_grid(a)

        part = burst_grid(a, 2:4)
        @test length(part.placements) == 3
        @test [p.burst for p in part.placements] == [2, 3, 4]
        @test part.nlines < whole.nlines
        @test part.nsamples == whole.nsamples

        # A subset is anchored on its own first burst, so it starts at row 1 rather than carrying the
        # offset it had in the whole subswath.
        @test first(part.placements[1].grid_rows) ==
              first(burst_grid(a, 2:2).placements[1].grid_rows)

        one = burst_grid(a, 5:5)
        @test length(one.placements) == 1
        @test one.nlines == a.lines_per_burst
        @test one.placements[1].burst_rows == valid_lines(a, 5)

        @test_throws "not a range of them" burst_grid(a, 0:3)
        @test_throws "not a range of them" burst_grid(a, 1:(n + 1))
        @test_throws "at least one burst" burst_grid(a, 3:2)
    end

    # A product whose bursts do not land on a shared axis would need interpolation, which this does
    # not do; placing them at whole rows anyway would shift their samples in azimuth silently.
    @testset "a burst off the uniform grid is refused" begin
        case = first(BURST_GRID_CASES)
        a = case.annotation
        @test_throws "from the uniform azimuth grid" burst_grid(a; tolerance = 1.0e-9)
    end
end
