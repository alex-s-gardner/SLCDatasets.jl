# The NISAR reader against the committed fixture.
#
# Floats are compared with `===` against the fixture's hex literals, so these assert bit-exact
# agreement with what `h5py` read from the granule, not agreement to a printed precision.

using SLCDatasets
using SLCDatasets: LookSide, LookLeft, LookRight, SPEED_OF_LIGHT, parse_cf_epoch, NisarBackend, nisar_band,
           nisar_product_type, GEOCODED_TYPES, grid
using Dates
using HDF5
using StaticArrays: SVector
using Test

const FX = FIXTURE

mktempdir() do dir
    path = write_fixture_product(joinpath(dir, "fixture_rslc.h5"))
    s = open_slc(path)

    @testset "identification" begin
        id = s.identification
        @test id.mission == FX.identification.mission
        @test id.product_type == FX.identification.product_type
        @test id.absolute_orbit == FX.identification.absolute_orbit
        @test id.pass_direction == FX.identification.pass_direction
        @test id.look_direction == FX.identification.look_direction
        # Kept as strings: the product records nanoseconds, which a `DateTime` would silently drop.
        @test id.start_time == FX.identification.start_time
        @test id.stop_time == FX.identification.stop_time
        @test start_datetime(id) == DateTime("2025-10-28T23:52:01.000")
        @test stop_datetime(id) == DateTime("2025-10-28T23:52:38.999")
    end

    @testset "geometry is bitwise" begin
        g = s.geometry
        @test g.starting_range === gx(FX.geometry.starting_range)
        @test g.far_range === gx(FX.geometry.far_range)
        @test g.range_pixel_spacing === gx(FX.geometry.range_pixel_spacing)
        @test g.wavelength === gx(FX.geometry.wavelength)
        @test g.prf === gx(FX.geometry.prf)
        @test g.sensing_start === gx(FX.geometry.sensing_start)
        @test g.sensing_stop === gx(FX.geometry.sensing_stop)
        @test g.nlines == FX.geometry.nlines
        @test g.nsamples == FX.geometry.nsamples
        @test g.look_side == (lowercase(FX.identification.look_direction) == "left" ?
                              LookLeft : LookRight)
        @test nlines(s) == g.nlines
        @test nsamples(s) == g.nsamples
    end

    @testset "wavelength is c over the center frequency" begin
        # The reference derives the wavelength this way, so the division must not be reordered.
        fc = SPEED_OF_LIGHT / gx(FX.geometry.wavelength)
        @test s.geometry.wavelength === SPEED_OF_LIGHT / fc
    end

    @testset "orbit is bitwise" begin
        o = orbit(s)
        @test length(o.time) == FX.orbit.n
        @test o.interp_method == FX.orbit.interp_method
        @test o.kind == FX.orbit.kind
        @test all(o.time[i] === gx(FX.orbit.time[i]) for i in eachindex(o.time))
        # The file stores (N, 3); a transposed read would swap x and z here and pass a norm check.
        @test all(o.position[i][c] === gx(FX.orbit.position[i][c])
                  for i in eachindex(o.position), c in 1:3)
        @test all(o.velocity[i][c] === gx(FX.orbit.velocity[i][c])
                  for i in eachindex(o.velocity), c in 1:3)
    end

    @testset "the orbit is read once and held" begin
        s2 = open_slc(path)
        @test s2.orbit === nothing
        first_read = orbit(s2)
        @test s2.orbit !== nothing
        @test orbit(s2) === first_read
    end

    @testset "one epoch for both clocks" begin
        # The azimuth times and the state vector times share an epoch in this format. A consumer
        # converting between them relies on it, so a product where they diverge must be visible.
        @test s.geometry.epoch == orbit(s).epoch
        @test s.geometry.epoch == parse_cf_epoch(FX.geometry.epoch)
    end

    @testset "the orbit brackets the acquisition" begin
        o = orbit(s)
        @test first(o.time) <= s.geometry.sensing_start
        @test s.geometry.sensing_stop <= last(o.time)
    end

    @testset "backend paths" begin
        b = s.backend
        @test b isa NisarBackend
        @test b.band == FX.band
        @test b.product_type == FX.product_type
        @test b.frequency == FX.frequency
        @test nisar_band(path) == FX.band
        @test nisar_product_type(path, FX.band) == FX.product_type
    end

    @testset "show" begin
        @test occursin("RSLC", sprint(show, s))
        long = sprint(show, MIME"text/plain"(), s)
        @test occursin("left", long)
        @test occursin("35 state vectors", sprint(show, MIME"text/plain"(), s))
    end
end

@testset "parse_cf_epoch" begin
    @test parse_cf_epoch("seconds since 2025-10-28T00:00:00") == DateTime(2025, 10, 28)
    @test parse_cf_epoch("seconds since 2025-10-28 00:00:00") == DateTime(2025, 10, 28)
    @test parse_cf_epoch("seconds since 2025-10-28T00:00:00Z") == DateTime(2025, 10, 28)
    # Rescaling a different unit would corrupt every time by a constant factor, so it throws.
    @test_throws "only seconds are supported" parse_cf_epoch("days since 2025-10-28T00:00:00")
    @test_throws "cannot read a reference epoch" parse_cf_epoch("2025-10-28")
end

@testset "unreadable inputs throw" begin
    @test_throws "is not a readable file" open_slc(joinpath(@__DIR__, "no_such_file.h5"))
    mktempdir() do dir
        plain = joinpath(dir, "plain.txt")
        write(plain, "not hdf5")
        @test_throws "neither an HDF5 file nor a Sentinel-1 SAFE product" open_slc(plain)

        empty = joinpath(dir, "empty.h5")
        h5open(empty, "w") do h
            h["unrelated"] = 1
        end
        @test_throws "not a NISAR-format product" open_slc(empty)
    end
end

# The stored width of a field is a property of how the product was written, not of the format: counts
# appear as `Int32` or `Int64` and times as `Float32` or `Float64` across generations. The reader
# converts, so a narrower product reads to the same values rather than failing or reinterpreting bytes.
@testset "field widths the product chose do not change what is read" begin
    mktempdir() do dir
        wide = write_fixture_product(joinpath(dir, "wide.h5"))
        reference = open_slc(wide)
        ref_orbit = orbit(reference)

        narrow = joinpath(dir, "narrow.h5")
        write_fixture_product(narrow)
        # Rewrite the orbit and the azimuth axis at single precision, and the orbit count as `Int64`.
        h5open(narrow, "r+") do h
            p = "science/LSAR/$(FIXTURE.product_type)"
            for name in ("$p/metadata/orbit/time", "$p/metadata/orbit/position",
                         "$p/metadata/orbit/velocity", "$p/swaths/zeroDopplerTime")
                units = haskey(HDF5.attributes(h[name]), "units") ?
                        read_attribute(h[name], "units") : nothing
                data = Float32.(read(h[name]))
                delete_object(h, name)
                h[name] = data
                units === nothing || write_attribute(h[name], "units", units)
            end
            delete_object(h, "science/LSAR/identification/absoluteOrbitNumber")
            h["science/LSAR/identification/absoluteOrbitNumber"] =
                Int64(reference.identification.absolute_orbit)
        end

        s = open_slc(narrow)
        @test s.identification.absolute_orbit == reference.identification.absolute_orbit
        @test s.geometry.nlines == reference.geometry.nlines
        @test s.geometry.nsamples == reference.geometry.nsamples
        @test s.geometry.epoch == reference.geometry.epoch
        # Single precision is what was stored, so the values agree only to that.
        @test s.geometry.sensing_start ≈ reference.geometry.sensing_start rtol = 1e-6
        o = orbit(s)
        @test length(o.time) == length(ref_orbit.time)
        @test o.epoch == ref_orbit.epoch
        @test eltype(o.time) === Float64
        @test eltype(o.position) === SVector{3,Float64}
        for i in eachindex(o.position)
            @test o.position[i] ≈ ref_orbit.position[i] rtol = 1e-6
            @test o.velocity[i] ≈ ref_orbit.velocity[i] rtol = 1e-6
        end
    end
end

# The epoch lives in a `units` attribute rather than a dataset, so a product missing it would otherwise
# be read with whatever epoch the parse of an empty string produced.
@testset "a missing units attribute is reported" begin
    mktempdir() do dir
        path = write_fixture_product(joinpath(dir, "no_units.h5"))
        h5open(path, "r+") do h
            delete_attribute(h["science/LSAR/$(FIXTURE.product_type)/swaths/zeroDopplerTime"],
                             "units")
        end
        @test_throws "has no `units` attribute" open_slc(path)
    end
end

# A name the layout says is a dataset but the file made a group would otherwise reach `read` and fail
# inside HDF5; the reader names the path instead.
@testset "a group where a dataset belongs is reported" begin
    mktempdir() do dir
        path = write_fixture_product(joinpath(dir, "wrong_kind.h5"))
        h5open(path, "r+") do h
            name = "science/LSAR/identification/absoluteOrbitNumber"
            delete_object(h, name)
            create_group(h, name)
        end
        @test_throws "is not a NISAR-format product" open_slc(path)
    end
end

@testset "a geocoded product is rejected, not misread" begin
    # A GSLC carries map-projected grids and no slant-range axis. Reading one as if it had a swath
    # would fail on a missing group; the reader says why instead.
    mktempdir() do dir
        path = joinpath(dir, "gslc.h5")
        h5open(path, "w") do h
            h["science/LSAR/identification/productType"] = "GSLC"
            h["science/LSAR/identification/absoluteOrbitNumber"] = Int32(1)
            h["science/LSAR/identification/listOfFrequencies"] = ["A"]
            h["science/LSAR/GSLC/grids/frequencyA/x"] = [0.0, 1.0]
        end
        @test "GSLC" in GEOCODED_TYPES
        @test_throws "carries no slant-range/azimuth geometry" open_slc(path)
    end
end

@testset "samples" begin
    # `(line, sample)` on the way out, whatever the file's layout: the dataset is written transposed, so
    # a reader that forgot to transpose would return this array's transpose and every assertion here
    # would fail on a shape.
    A = ComplexF32[complex(i, j) for i in 1:7, j in 1:5]
    mktempdir() do dir
        path = write_fixture_product(joinpath(dir, "rslc_with_samples.h5"); samples = A)
        s = open_slc(path)
        p = pixels(s)

        @test p isa SLCDatasets.NisarRaster{ComplexF32}
        @test size(p) == (7, 5)
        @test p[:, :] == A
        @test p[3, 4] === A[3, 4]
        @test p[2:5, 2:3] == A[2:5, 2:3]
        @test p[2:5, 3] == A[2:5, 3]
        @test p[4, 2:4] == A[4, 2:4]
        # A window is the same values as the whole read, which is what says the hyperslab offsets are
        # applied to the axis they belong to rather than to the other one.
        @test all(p[i:i, j:j][1, 1] === A[i, j] for i in axes(A, 1), j in axes(A, 2))

        @test amplitude(s) == abs.(A)
        @test_throws BoundsError p[1:8, 1:5]
    end
end

@testset "a product carrying no samples says so" begin
    # The fixture without a sample array is the state a metadata-only read leaves: the polarization is
    # listed and the dataset is absent, which is a different failure from an unreadable file.
    mktempdir() do dir
        path = write_fixture_product(joinpath(dir, "rslc_no_samples.h5"))
        @test_throws "has no HH" pixels(open_slc(path))
    end
end

@testset "a geocoded product opens through its own door" begin
    A = ComplexF32[complex(i, -j) for i in 1:4, j in 1:6]
    mktempdir() do dir
        path = write_fixture_geocoded(joinpath(dir, "gslc.h5"); samples = A)
        g = open_geocoded(path)

        @test g isa GeocodedProduct
        @test g.identification.product_type == "GSLC"
        @test nlines(g) == 4
        @test nsamples(g) == 6
        # The product stores pixel *centers* and a grid names the outer corner, so the origin sits half a
        # pixel out along each axis — the convention a geotransform states and an intersection needs.
        @test grid(g).origin == (-340558.75 - 1.25, -2.1067225e6 + 2.5)
        @test grid(g).spacing == (2.5, -5.0)
        @test grid(g).size == (4, 6)
        @test grid(g).epsg == 3413
        @test pixels(g)[:, :] == A
        @test pixels(g)[2:3, 4:5] == A[2:3, 4:5]

        # The two doors stay separate: a geocoded product has no slant-range axis to describe, and a
        # radar-geometry one has no map grid.
        @test_throws "carries no slant-range/azimuth geometry" open_slc(path)
        rslc = write_fixture_product(joinpath(dir, "rslc.h5"))
        @test_throws "lies in radar geometry" open_geocoded(rslc)
    end
end
