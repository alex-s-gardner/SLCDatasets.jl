# Rebuilding a NISAR-layout product from the committed fixture.
#
# The fixture records what `h5py` read from a real granule. Writing it back out in the same group
# layout gives the reader a product to open, so the whole path — root-path probe, product-type
# dispatch, group paths, epoch parsing, the (3, N) orbit transpose — is exercised without a multi-
# gigabyte file or a network.

using HDF5
using JSON3
using SLCDatasets: SPEED_OF_LIGHT

const FIXTURE = JSON3.read(read(joinpath(@__DIR__, "reference", "nisar_metadata.json"), String))

# Every float in the fixture carries a hex literal beside its decimal, and the hex is what the tests
# compare against: a decimal round-trip through JSON is not guaranteed to preserve the last bit.
gx(v) = parse(Float64, v.hex)

# The inverse, for a test that needs a value the granule does not carry.
hx(x::Real) = (dec = Float64(x), hex = _hex_literal(Float64(x)))

function _hex_literal(v::Float64)
    v == 0.0 && return signbit(v) ? "-0x0.0p+0" : "0x0.0p+0"
    sign = signbit(v) ? "-" : ""
    m, e = frexp(abs(v))
    frac = m * 2 - 1.0
    digits = ""
    for _ in 1:13
        frac *= 16
        d = floor(Int, frac)
        digits *= string(d; base = 16)
        frac -= d
    end
    exponent = e - 1
    return string(sign, "0x1.", digits, "p", exponent >= 0 ? "+" : "-", abs(exponent))
end

"""
    override(fx, path => value, ...) -> NamedTuple

`fx` with the named fields replaced, for a test needing a product the granule is not.

`path` is a `Symbol` for a top-level field or a `Symbol` pair for a nested one, as in
`:geometry => :epoch`.
"""
function override(fx, pairs::Pair...)
    out = Dict{Symbol,Any}(k => v for (k, v) in Base.pairs(fx))
    for (key, value) in pairs
        if key isa Pair
            outer, inner = key
            nested = Dict{Symbol,Any}(k => v for (k, v) in Base.pairs(out[outer]))
            nested[inner] = value
            out[outer] = NamedTuple(nested)
        else
            out[key] = value
        end
    end
    return NamedTuple(out)
end

"""
    write_fixture_product(path, fx = FIXTURE) -> String

Write `fx` as a NISAR-layout HDF5 product at `path`.

Only the datasets the reader reads are written, with the group names and the `units` attributes a real
product carries.
"""
function write_fixture_product(path::AbstractString, fx = FIXTURE; samples = nothing)
    band, product, freq = fx.band, fx.product_type, fx.frequency
    id, geom, orb = fx.identification, fx.geometry, fx.orbit

    n = orb.n
    time = [gx(v) for v in orb.time]
    # The fixture stores state vectors as n rows of 3, which is the file's own (N, 3) layout; HDF5.jl
    # writes a (3, n) Julia array to that shape.
    pos = [gx(orb.position[i][c]) for c in 1:3, i in 1:n]
    vel = [gx(orb.velocity[i][c]) for c in 1:3, i in 1:n]

    h5open(path, "w") do h
        p = "science/$band/$product"
        ip = "science/$band/identification"

        h["$ip/missionId"] = id.mission
        h["$ip/productType"] = id.product_type
        h["$ip/absoluteOrbitNumber"] = Int32(id.absolute_orbit)
        h["$ip/orbitPassDirection"] = id.pass_direction
        h["$ip/lookDirection"] = id.look_direction
        h["$ip/zeroDopplerStartTime"] = id.start_time
        h["$ip/zeroDopplerEndTime"] = id.stop_time
        h["$ip/boundingPolygon"] = "POLYGON EMPTY"
        h["$ip/listOfFrequencies"] = [freq]

        # The azimuth axis is regenerated from its ends and the line count so the file carries a full
        # `zeroDopplerTime` vector of the product's real length, as `nlines` is read from its size.
        t0, t1, nl = gx(geom.sensing_start), gx(geom.sensing_stop), geom.nlines
        h["$p/swaths/zeroDopplerTime"] = collect(range(t0, t1; length = nl))
        h["$p/swaths/zeroDopplerTimeSpacing"] = 1 / gx(geom.prf)

        r0, r1, ns = gx(geom.starting_range), gx(geom.far_range), geom.nsamples
        h["$p/swaths/frequency$freq/slantRange"] = collect(range(r0, r1; length = ns))
        h["$p/swaths/frequency$freq/slantRangeSpacing"] = gx(geom.range_pixel_spacing)
        h["$p/swaths/frequency$freq/processedCenterFrequency"] =
            SPEED_OF_LIGHT / gx(geom.wavelength)
        h["$p/swaths/frequency$freq/listOfPolarizations"] = ["HH"]

        h["$p/metadata/orbit/time"] = time
        h["$p/metadata/orbit/position"] = pos
        h["$p/metadata/orbit/velocity"] = vel
        h["$p/metadata/orbit/interpMethod"] = orb.interp_method
        h["$p/metadata/orbit/orbitType"] = orb.kind

        write_attribute(h["$p/swaths/zeroDopplerTime"], "units", geom.epoch)
        write_attribute(h["$p/metadata/orbit/time"], "units", orb.epoch)

        # The samples, written in the file's own layout: azimuth slowest, so the dataset is the
        # transpose of the `(line, sample)` array a reader hands back. Small, because what is under test
        # is the read and the transpose rather than a product's real extent.
        samples === nothing || (h["$p/swaths/frequency$freq/HH"] = permutedims(samples))
    end
    return path
end

"""
    write_fixture_geocoded(path; samples, origin, spacing, epsg) -> String

A minimal geocoded NISAR product: a `grids` group with its coordinate axes and a sample array.

`origin` names the **first pixel's center**, as a product does; the grid a reader reports pulls that back
by half a pixel to the corner.
"""
function write_fixture_geocoded(path::AbstractString; samples::AbstractMatrix,
                                origin::Tuple{Real,Real} = (-340558.75, -2.1067225e6),
                                spacing::Tuple{Real,Real} = (2.5, -5.0), epsg::Integer = 3413)
    ny, nx = size(samples)
    h5open(path, "w") do h
        ip = "science/LSAR/identification"
        h["$ip/missionId"] = "NISAR"
        h["$ip/productType"] = "GSLC"
        h["$ip/absoluteOrbitNumber"] = Int32(1)
        h["$ip/orbitPassDirection"] = "Descending"
        h["$ip/lookDirection"] = "Right"
        h["$ip/zeroDopplerStartTime"] = "2025-10-28T23:52:01.000000000"
        h["$ip/zeroDopplerEndTime"] = "2025-10-28T23:52:38.000000000"
        h["$ip/boundingPolygon"] = "POLYGON EMPTY"
        h["$ip/listOfFrequencies"] = ["A"]

        g = "science/LSAR/GSLC/grids/frequencyA"
        h["$g/xCoordinates"] = collect(range(origin[1]; step = spacing[1], length = nx))
        h["$g/yCoordinates"] = collect(range(origin[2]; step = spacing[2], length = ny))
        h["$g/xCoordinateSpacing"] = Float64(spacing[1])
        h["$g/yCoordinateSpacing"] = Float64(spacing[2])
        h["$g/projection"] = UInt32(epsg)
        h["$g/listOfPolarizations"] = ["HH"]
        h["$g/HH"] = permutedims(samples)
    end
    return path
end
