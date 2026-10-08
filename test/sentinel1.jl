# The Sentinel-1 reader against golden values dumped from isce3 + s1reader.
#
# The scalars the two must agree on bitwise are compared with `===` against the reference's hex
# literals, so these assert agreement with what ISCE3 computed rather than agreement to a printed
# precision. Times are compared as absolute instants — `epoch` plus the offset — because this package's
# epoch is truncated to the second where the reference's keeps microseconds, and only the instant they
# denote is a claim about the product.
#
# The granules themselves are multi-gigabyte and not in the repository, so the products these run
# against are rebuilt from `reference/sentinel1_inputs*.json`, which carries the annotation fields and
# state vectors of each granule verbatim — the text as the annotation writes it, so the reader parses
# what it would have parsed from the granule. The comparison against ISCE3 is therefore hardcoded and
# runs anywhere.
#
# Setting `SLCDATASETS_S1_DIR` to a directory holding the granules named in the golden files reads those
# instead. Both modes assert the same values, which is what keeps the rebuilt products honest: if a
# reader change made them diverge, the granule run would fail where the fixture run passed.

using SLCDatasets
using SLCDatasets: LookRight, SPEED_OF_LIGHT, Sentinel1Backend, UtcTime, seconds_between,
           parse_utc, read_eof_state_vectors, safe_polarizations, default_polarization,
           is_safe_product, annotation_xml, epoch_of, S1_ORBIT_PADDING,
           read_annotation, valid_lines, valid_samples
using Dates
using StaticArrays: SVector
using JSON
using Test

norm3(v) = sqrt(v[1]^2 + v[2]^2 + v[3]^2)

# An instant named as an epoch plus an offset in seconds, which is how both readers report their times.
struct Instant
    epoch::UtcTime
    offset::Float64
end

# Ours, whose epoch is a whole second, and the reference's, whose epoch string carries microseconds.
ours(epoch::DateTime, offset::Float64) = Instant(UtcTime(epoch, 0.0), offset)
reference(epoch_string, offset::Float64) = Instant(parse_utc(epoch_string), offset)
utc(stamp) = Instant(parse_utc(stamp), 0.0)

# Seconds from `b` to `a`, taking the difference of the epochs before adding the offsets. Differencing
# the epochs first keeps every intermediate value the size of an acquisition rather than of the seconds
# since some distant origin, where a `Float64` no longer resolves a nanosecond.
seconds_apart(a::Instant, b::Instant) =
    seconds_between(b.epoch, a.epoch) + (a.offset - b.offset)

# The products the golden values are checked against. By default each is rebuilt from the committed
# inputs, which carry the annotation fields and state vectors of the granule verbatim. Setting
# `SLCDATASETS_S1_DIR` to a directory holding the granules themselves reads those instead, which is what
# confirms the rebuilt products still stand for them.
const S1_GRANULE_DIR = get(ENV, "SLCDATASETS_S1_DIR", "")
const S1_FROM_GRANULES = !isempty(S1_GRANULE_DIR)
const S1_FIXTURE_DIR = mktempdir(; cleanup = true)

function s1_products()
    out = NamedTuple[]
    for f in s1_fixtures()
        if S1_FROM_GRANULES
            safe = joinpath(S1_GRANULE_DIR, String(f.gold.safe))
            eof = joinpath(S1_GRANULE_DIR, String(f.gold.orbit_file))
            # A granule named by the golden file but absent from the directory is a gap in coverage,
            # not a pass; `@test` on the paths below makes it visible.
            push!(out, (; f.gold, f.inputs, safe, eof, f.name))
        else
            dir = mkpath(joinpath(S1_FIXTURE_DIR, splitext(f.name)[1]))
            safe, eof = write_s1_fixture(dir, f.inputs)
            push!(out, (; f.gold, f.inputs, safe, eof, f.name))
        end
    end
    return out
end

const S1_PRODUCTS = s1_products()

if isempty(S1_PRODUCTS)
    @info "no Sentinel-1 golden-value files found, so those tests are skipped."
else
    S1_FROM_GRANULES ||
        @info """Sentinel-1 golden values checked against products rebuilt from
                 test/reference/sentinel1_inputs*.json. Set SLCDATASETS_S1_DIR to a directory holding
                 the granules named in the golden files to check against those instead."""
    for p in S1_PRODUCTS
        gold = p.gold
        @testset "$(gold.safe)" begin
            # A granule named by a golden file but absent from `SLCDATASETS_S1_DIR` is a gap in
            # coverage rather than a pass, so the paths are asserted before anything reads them.
            @test ispath(p.safe)
            @test isfile(p.eof)
            (ispath(p.safe) && isfile(p.eof)) || continue

            @test is_safe_product(p.safe)
            pol = String(gold.polarization)
            @test pol in safe_polarizations(p.safe)
            @test default_polarization(p.safe) == pol

            @testset "subswath mosaic matches isce3" begin
                m = gold.mosaic
                s = open_slc(p.safe; orbit = p.eof, polarization = pol)
                g = s.geometry

                # Ranges, spacings and the wavelength are derived from the annotation by the same
                # arithmetic isce3 uses, so they agree to the bit.
                @test g.starting_range === gx(m.starting_range)
                @test g.far_range === gx(m.far_range)
                @test g.range_pixel_spacing === gx(m.range_pixel_spacing)
                @test g.wavelength === gx(m.wavelength)
                @test g.prf === gx(m.prf)
                @test g.nlines == m.nlines
                @test g.nsamples == m.nsamples
                @test g.look_side == LookRight
                @test nlines(s) == m.nlines
                @test nsamples(s) == m.nsamples

                # The mosaic spans the subswaths, so its first line precedes IW1's: the epoch anchors
                # to IW1's first burst but the sensing window opens earlier.
                @test seconds_apart(ours(g.epoch, g.sensing_start),
                                    utc(m.sensing_start_utc)) ≈ 0 atol = 1e-9
                # The stop time is reached by accumulating ~13000 steps of `1/prf`, so its tolerance is
                # that accumulation rather than a differently-read field.
                @test seconds_apart(ours(g.epoch, g.sensing_stop),
                                    utc(m.sensing_stop_utc)) ≈ 0 atol = 1e-6
                # The same instant against the reference's own epoch and offset, which is the pair a
                # consumer of the ISCE3 metadata actually sees.
                @test seconds_apart(ours(g.epoch, g.sensing_start),
                                    reference(gold.orbit.epoch, gx(m.sensing_start))) ≈ 0 atol = 1e-9
            end

            @testset "mosaic orbit matches isce3" begin
                s = open_slc(p.safe; orbit = p.eof, polarization = pol)
                o = orbit(s)
                ref = gold.orbit

                # The window is the anchor burst's, padded, so the count is a property of the
                # filtering rather than of the file: a wider or narrower window changes it.
                @test length(o.time) == ref.n
                @test length(o.position) == ref.n
                @test length(o.velocity) == ref.n
                @test o.epoch == s.geometry.epoch
                @test o.kind == "POEORB"

                for i in eachindex(o.time)
                    @test seconds_apart(ours(o.epoch, o.time[i]),
                                        reference(ref.epoch, gx(ref.time[i]))) ≈ 0 atol = 1e-9
                end
                # Positions and velocities are parsed decimals, not derived, so any difference here is
                # a misread column rather than arithmetic. A transposed read would swap x and z and
                # still pass a norm check, which is why every component is compared.
                for i in eachindex(o.position), c in 1:3
                    @test o.position[i][c] === gx(ref.position[i][c])
                    @test o.velocity[i][c] === gx(ref.velocity[i][c])
                end
            end

            @testset "every burst matches isce3" begin
                for (key, bursts) in pairs(gold.bursts)
                    swath = parse(Int, String(key))
                    @test nbursts(p.safe; swath, polarization = pol) == length(bursts)
                    for (i, ref) in enumerate(bursts)
                        s = open_slc(p.safe; orbit = p.eof, swath, burst = i, polarization = pol)
                        g = s.geometry
                        @test g.starting_range === gx(ref.starting_range)
                        @test g.far_range === gx(ref.far_range)
                        @test g.range_pixel_spacing === gx(ref.range_pixel_spacing)
                        @test g.wavelength === gx(ref.wavelength)
                        @test g.prf === gx(ref.prf)
                        @test g.nlines == ref.nlines
                        @test g.nsamples == ref.nsamples
                        @test g.look_side == LookRight
                        # A burst's epoch is its own sensing start less two days, so this is the test
                        # that a burst is not silently read with a neighbor's timing.
                        @test seconds_apart(ours(g.epoch, g.sensing_start),
                                            utc(ref.sensing_start_utc)) ≈ 0 atol = 1e-9
                        @test seconds_apart(ours(g.epoch, g.sensing_stop),
                                            utc(ref.sensing_start_utc)) ≈
                              (ref.nlines - 1) * gx(ref.azimuth_time_interval) atol = 1e-9
                        # And against the reference's own epoch/offset pair, which for a burst is its
                        # subswath's epoch rather than the mosaic's.
                        @test seconds_apart(ours(g.epoch, g.sensing_start),
                                            reference(ref.epoch,
                                                      gx(ref.sensing_start))) ≈ 0 atol = 1e-9
                    end
                end
            end

            @testset "the epoch is the anchor burst less two days" begin
                # The offset is not free to change: an ISCE3-derived Doppler LUT carries a matching
                # one, and disagreeing puts every azimuth time exactly two days out.
                iw1 = gold.bursts[Symbol(1)][1]
                s = open_slc(p.safe; orbit = p.eof, swath = 1, burst = 1, polarization = pol)
                anchor = parse_utc(iw1.sensing_start_utc)
                @test s.geometry.epoch == epoch_of(anchor)
                @test Dates.Day(anchor.datetime - s.geometry.epoch) == Dates.Day(2)
                # The mosaic reports against IW1's first burst, not its own earliest line.
                mosaic = open_slc(p.safe; orbit = p.eof, polarization = pol)
                @test mosaic.geometry.epoch == epoch_of(anchor)
            end

            @testset "identification" begin
                s = open_slc(p.safe; orbit = p.eof, polarization = pol)
                id = s.identification
                @test id.mission == uppercase(first(gold.safe, 3))
                @test id.product_type == "SLC"
                # The granule name carries the absolute orbit two fields before the product's unique
                # ID, zero-padded to six digits.
                @test id.absolute_orbit ==
                      parse(Int, split(splitext(String(gold.safe))[1], '_')[end - 2])
                @test id.look_direction == "Right"
                @test id.pass_direction in ("ascending", "descending")
                # The mosaic's bounds are its own, so they match the geometry rather than one
                # subswath's `productFirstLineUtcTime`.
                @test seconds_apart(utc(id.start_time),
                                    ours(s.geometry.epoch, s.geometry.sensing_start)) ≈ 0 atol = 1e-6
                @test start_datetime(id) <= stop_datetime(id)
            end

            @testset "invariants that need no reference" begin
                s = open_slc(p.safe; orbit = p.eof, polarization = pol)
                g = s.geometry
                # Range axis closes on itself: a misread spacing or width breaks this.
                @test g.far_range ≈ g.starting_range + (g.nsamples - 1) * g.range_pixel_spacing
                @test g.sensing_stop > g.sensing_start
                @test g.nlines ≈ (g.sensing_stop - g.sensing_start) * g.prf + 1 rtol = 1e-6
                @test g.wavelength ≈ 0.0554 atol = 1e-4   # C-band

                o = orbit(s)
                @test issorted(o.time)
                # The state vectors bracket the anchor burst, and the padding bounds how far past it
                # they may reach.
                @test first(o.time) <= g.sensing_start
                @test last(o.time) >= g.sensing_start
                @test all(isfinite, o.time)
                @test all(v -> 6.9e6 < norm3(v) < 7.2e6, o.position)   # LEO radius
                @test all(v -> 7.0e3 < norm3(v) < 7.9e3, o.velocity)   # LEO speed
            end

            @testset "the orbit is read once and held" begin
                s = open_slc(p.safe; orbit = p.eof, polarization = pol)
                @test s.orbit === nothing
                first_read = orbit(s)
                @test s.orbit !== nothing
                @test orbit(s) === first_read
            end

            @testset "subswath selection" begin
                # A single subswath is narrower than the mosaic and starts at its own range.
                mosaic = open_slc(p.safe; orbit = p.eof, polarization = pol)
                iw1 = open_slc(p.safe; orbit = p.eof, swath = 1, polarization = pol)
                iw3 = open_slc(p.safe; orbit = p.eof, swath = 3, polarization = pol)
                @test iw1.geometry.nsamples < mosaic.geometry.nsamples
                @test iw1.geometry.starting_range === mosaic.geometry.starting_range
                @test iw3.geometry.starting_range > iw1.geometry.starting_range
                @test open_slc(p.safe; orbit = p.eof, swaths = [1, 2, 3],
                               polarization = pol).geometry.nsamples == mosaic.geometry.nsamples
            end

            @testset "show" begin
                s = open_slc(p.safe; orbit = p.eof, polarization = pol)
                @test occursin("SLC", sprint(show, s))
                long = sprint(show, MIME"text/plain"(), s)
                @test occursin("right", long)
                @test occursin("not read", long)
            end
        end
    end
end

@testset "parse_utc" begin
    t = parse_utc("2015-11-20T08:02:04.790228")
    @test t.datetime == DateTime("2015-11-20T08:02:04")
    # The microseconds are kept out of the `DateTime`, which holds only milliseconds.
    @test t.seconds ≈ 0.790228 atol = 1e-12
    @test parse_utc("2015-11-20T08:02:04").seconds === 0.0
    @test parse_utc("2015-11-20T08:02:04Z").datetime == DateTime("2015-11-20T08:02:04")
    @test_throws "cannot read a UTC time" parse_utc("2015-11-20")
end

@testset "seconds_between is exact across a second boundary" begin
    # Subtracting two `DateTime`s would round the sub-second parts away; a time that straddles a
    # second must come out exact.
    a = parse_utc("2015-11-20T08:02:04.900000")
    b = parse_utc("2015-11-20T08:02:05.100000")
    @test seconds_between(a, b) ≈ 0.2 atol = 1e-12
    @test seconds_between(b, a) ≈ -0.2 atol = 1e-12
    @test seconds_between(a, a) === 0.0
end

if !isempty(S1_PRODUCTS)
    @testset "a Sentinel-1 product needs an orbit file" begin
        p = first(S1_PRODUCTS)
        # Silently returning an SLC whose state vectors are unavailable would defer the failure to
        # whatever consumes the orbit, so opening one is refused outright.
        @test_throws "state vectors live in a separate" open_slc(p.safe)
        @test_throws "is not a readable orbit file" open_slc(p.safe; orbit = "no_such.EOF")
        @test_throws "not both" open_slc(p.safe; orbit = p.eof, swath = 1, swaths = [1, 2])
        @test_throws "it needs `swath = 1`" open_slc(p.safe; orbit = p.eof, burst = 1)
        @test_throws "but 2 were named" open_slc(p.safe; orbit = p.eof, swaths = [1, 2], burst = 1)
        @test_throws "subswaths are 1, 2 and 3" open_slc(p.safe; orbit = p.eof, swath = 4)
        @test_throws "does not exist" open_slc(p.safe; orbit = p.eof, swath = 1, burst = 9999)
        @test_throws "repeats a subswath" open_slc(p.safe; orbit = p.eof, swaths = [1, 1])
        @test_throws "name at least one subswath" open_slc(p.safe; orbit = p.eof, swaths = Int[])
        wrong = String(first(S1_PRODUCTS).gold.polarization) == "hh" ? "vv" : "hh"
        @test_throws "not $(uppercase(wrong))" open_slc(p.safe; orbit = p.eof,
                                                        polarization = wrong)
    end

    @testset "an orbit file from the wrong granule is rejected" begin
        ps = S1_PRODUCTS
        if length(ps) >= 2
            # The two products are years apart, so neither orbit file covers the other's window.
            # Reading zero state vectors and reporting an empty orbit would be the silent failure.
            @test_throws "different granule" orbit(open_slc(ps[1].safe; orbit = ps[2].eof,
                                                           polarization = String(ps[1].gold.polarization)))
        end
    end
end

@testset "a malformed orbit file throws" begin
    mktempdir() do dir
        path = joinpath(dir, "bad.EOF")
        write(path, "<Earth_Explorer_File><Data_Block/></Earth_Explorer_File>")
        @test_throws "not a Sentinel-1 orbit file" read_eof_state_vectors(path)
    end
end

@testset "a non-SAFE input is not taken for one" begin
    mktempdir() do dir
        plain = joinpath(dir, "plain.txt")
        write(plain, "not a product")
        @test !is_safe_product(plain)
        @test !is_safe_product(dir)
        @test_throws "neither an HDF5 file nor a Sentinel-1 SAFE product" open_slc(plain)
        @test_throws "not a Sentinel-1 SAFE product" bursts(plain; orbit = plain)
    end
end

# The orbit window is half-open at its leading edge: a vector exactly `S1_ORBIT_PADDING` before the
# burst start is excluded. `s1reader` draws the boundary there, so admitting it would hand an
# interpolator one more vector than ISCE3 sees — a difference no golden-value field would reveal, since
# every kept vector's own values are unchanged.
@testset "the orbit padding window is half-open at its leading edge" begin
    start = parse_utc("2019-01-01T12:00:00.000000")
    stop = parse_utc("2019-01-01T12:00:03.000000")
    mktempdir() do dir
        path = joinpath(dir, "edge.EOF")
        # One vector exactly at `-padding`, one just inside it, and one exactly at `+padding`.
        stamps = ["2019-01-01T11:59:00.000000",   # start - 60 s exactly: excluded
                  "2019-01-01T11:59:00.500000",   # inside: kept
                  "2019-01-01T12:00:01.000000",   # inside: kept
                  "2019-01-01T12:01:03.000000"]   # stop + 60 s exactly: kept
        osvs = join(("<OSV><UTC>UTC=$s</UTC><X>1.0</X><Y>2.0</Y><Z>3.0</Z>" *
                     "<VX>4.0</VX><VY>5.0</VY><VZ>6.0</VZ></OSV>" for s in stamps))
        write(path, "<Earth_Explorer_File><Data_Block><List_of_OSVs count=\"4\">" *
                    osvs * "</List_of_OSVs></Data_Block></Earth_Explorer_File>")
        table = read_eof_state_vectors(path; from = start, to = stop,
                                       padding = S1_ORBIT_PADDING)
        @test length(table) == 3
        @test table.time[1] == parse_utc("2019-01-01T11:59:00.500000")
        # Unfiltered, every record is kept.
        @test length(read_eof_state_vectors(path)) == 4
    end
end

@testset "an OSV record missing a field throws rather than reading a zero vector" begin
    mktempdir() do dir
        path = joinpath(dir, "short.EOF")
        write(path, "<Earth_Explorer_File><Data_Block><List_of_OSVs count=\"1\">" *
                    "<OSV><UTC>UTC=2019-01-01T12:00:00.000000</UTC>" *
                    "<X>1.0</X><Y>2.0</Y><Z>3.0</Z></OSV>" *
                    "</List_of_OSVs></Data_Block></Earth_Explorer_File>")
        @test_throws "missing one of its X/Y/Z/VX/VY/VZ fields" read_eof_state_vectors(path)
    end
end

@testset "state vector records must be index-matched" begin
    @test_throws DimensionMismatch SLCDatasets.StateVectorTable(
        [parse_utc("2019-01-01T12:00:00")], SVector{3,Float64}[], SVector{3,Float64}[])
    @test_throws "index-matched" StateVectors(
        [0.0, 1.0], [SVector(1.0, 2.0, 3.0)], [SVector(1.0, 2.0, 3.0)],
        DateTime(2019), "Hermite", "POEORB")
end

# `bursts` shares one parse across a subswath, so every burst it yields must equal the one `open_slc`
# builds on its own; a shared-state bug would show as bursts differing from their independent reads.
if !isempty(S1_PRODUCTS)
    p = first(S1_PRODUCTS)
    pol = String(p.gold.polarization)
    @testset "bursts() agrees with open_slc per burst" begin
        for swath in 1:3
            series = bursts(p.safe; orbit = p.eof, swath, polarization = pol)
            @test series isa AbstractVector{<:SLC}
            @test length(series) == nbursts(p.safe; swath, polarization = pol)
            @test eachindex(series) == 1:length(series)
            @test_throws BoundsError series[length(series) + 1]
            for i in eachindex(series)
                direct = open_slc(p.safe; orbit = p.eof, swath, burst = i, polarization = pol)
                @test series[i].geometry == direct.geometry
                @test series[i].identification == direct.identification
                @test orbit(series[i]).time == orbit(direct).time
            end
        end
    end

    @testset "a Sentinel1Product refuses a subswath it was not opened for" begin
        prod = Sentinel1Product(p.safe; orbit = p.eof, polarization = pol, swaths = [1, 2])
        @test nbursts(prod, 2) == nbursts(p.safe; swath = 2, polarization = pol)
        @test_throws "was opened for subswaths IW1, IW2, not IW3" nbursts(prod, 3)
    end
end

# The valid region of each burst, against the reduction `s1reader`'s `burst_from_xml` performs. The
# golden values are the four scalars per burst held in the inputs, which is what the granules' per-line
# arrays reduce to; `read_valid_region` re-derives them from arrays the fixture writes back out, so
# this asserts the reduction rather than the round trip.
@testset "the valid region matches the s1reader reduction" begin
    for p in S1_PRODUCTS
        pol = String(p.gold.polarization)
        for swath in 1:3
            a = read_annotation(annotation_xml(p.safe, swath, pol), swath)
            regions = p.inputs["swaths"][string(swath)]["burstValidRegions"]
            @test length(regions) == nbursts(a)
            for (i, region) in enumerate(regions)
                # The inputs hold 1-based lines and the annotation's own 0-based samples.
                @test a.first_valid_line[i] == Int(region["firstValidLine"])
                @test a.last_valid_line[i] == Int(region["lastValidLine"])
                @test a.first_valid_sample[i] == Int(region["firstValidSample"]) + 1
                @test a.last_valid_sample[i] == Int(region["lastValidSample"]) + 1

                # The region must be inside the burst's own extent, or a merge clipping to it would
                # index outside the raster.
                @test 1 <= a.first_valid_line[i] <= a.last_valid_line[i] <= a.lines_per_burst
                @test 1 <= a.first_valid_sample[i] <= a.last_valid_sample[i] <= a.samples_per_burst

                # Inclusive bounds, matching `s1reader`'s `inlength`/`inwidth`.
                @test length(valid_lines(a, i)) ==
                      a.last_valid_line[i] - a.first_valid_line[i] + 1
                @test length(valid_samples(a, i)) ==
                      a.last_valid_sample[i] - a.first_valid_sample[i] + 1
            end
        end
    end
end

# Each of these is an annotation whose valid region is not the rectangle the reader reduces it to.
# Reducing one anyway yields a rectangle covering lines or samples that carry nothing, which a merge
# would then read as data.
@testset "a valid region that is not a rectangle is refused" begin
    lines = 6
    samples = 100
    burst(firsts, lasts) = """<product>
      <adsHeader><missionId>S1A</missionId><productType>SLC</productType>
        <polarisation>HH</polarisation><absoluteOrbitNumber>1</absoluteOrbitNumber></adsHeader>
      <generalAnnotation><productInformation><pass>Descending</pass>
        <rangeSamplingRate>6.4e7</rangeSamplingRate>
        <radarFrequency>5.4e9</radarFrequency></productInformation></generalAnnotation>
      <imageAnnotation><imageInformation><azimuthTimeInterval>0.002</azimuthTimeInterval>
        <slantRangeTime>0.005</slantRangeTime></imageInformation></imageAnnotation>
      <swathTiming><linesPerBurst>$lines</linesPerBurst>
        <samplesPerBurst>$samples</samplesPerBurst>
        <burstList count="1"><burst>
          <azimuthTime>2020-01-01T00:00:00.000000</azimuthTime>
          <firstValidSample>$firsts</firstValidSample>
          <lastValidSample>$lasts</lastValidSample>
        </burst></burstList></swathTiming>
    </product>"""

    valid = "-1 10 10 10 10 -1"
    ends = "-1 90 90 90 90 -1"

    # Nothing imaged at all.
    @test_throws "no valid lines" read_annotation(burst("-1 -1 -1 -1 -1 -1", ends), 2)

    # Valid lines in two runs: a gap in the middle is not a rectangle.
    @test_throws "more than one run" read_annotation(burst("-1 10 -1 10 10 -1", ends), 2)

    # One array shorter than the burst is tall.
    @test_throws "one per line" read_annotation(burst("-1 10 10 10 10", ends), 2)
    @test_throws "one of each" read_annotation(burst(valid, "-1 90 90 90 -1"), 2)

    # An interior line narrower than the two the reduction reads: taking the rectangle from the
    # boundary lines alone would claim samples this line does not carry.
    @test_throws "varies its valid sample range" read_annotation(
        burst(valid, "-1 90 50 90 90 -1"), 2)
    @test_throws "varies its valid sample range" read_annotation(
        burst("-1 10 40 10 10 -1", ends), 2)

    # A burst whose first and last valid lines share no sample.
    @test_throws "no sample valid on both" read_annotation(
        burst("-1 10 95 95 95 -1", "-1 20 99 99 99 -1"), 2)
end
