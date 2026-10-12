# The NISAR backend: HDF5 products under `/science/{LSAR,SSAR}/<productType>`.
#
# Group paths are built the way isce3's `nisar.products.readers.Base` builds them — a root path found
# by probing for the sensor band, then the product type, then a fixed subgroup name. The product type
# is read from the identification group and selects the reader, which is also how isce3's
# `open_product` dispatches. A product whose type is `SLC` is an early-mission RSLC.
#
# A product is read through libhdf5 for its metadata — the few datasets below, by path — and through
# ChunkManifests.jl's typed arrays for its samples, whose chunks are located once from the file's chunk index
# and then read with range requests. A local file and a remote one take the same path.

const SCIENCE_ROOT = "science"
const SENSOR_BANDS = ("LSAR", "SSAR")

# Speed of light, for converting the processed center frequency to a wavelength. This is the value
# isce3 uses (`isce3.core.speed_of_light`), and the conversion must match it to the bit.
const SPEED_OF_LIGHT = 299792458.0

"""
    NisarBackend <: AbstractSLCBackend

A NISAR-format HDF5 product: where it is, the transport it is read through, the sensor band group inside
it, and the product type.

`frequency` selects the sub-band (`"A"` or `"B"`); a product lists the ones it carries at
`identification/listOfFrequencies`.

No file is held open: each read opens the product, reads what it needs and closes it, so a backend has no
lifetime to manage. The samples are read through the transport alone once their chunks are located.
"""
struct NisarBackend <: AbstractSLCBackend
    path::String
    transport::LeafTransport
    band::String
    product_type::String
    frequency::String
end

# `g(f)` with the product open as `f`, closed afterwards.
function _with_file(g, b::NisarBackend)
    f = h5open_typed(b.path; transport = b.transport)
    try
        return g(f)
    finally
        close(f)
    end
end

product_path(b::NisarBackend) = string(SCIENCE_ROOT, "/", b.band, "/", b.product_type)
metadata_path(b::NisarBackend) = string(product_path(b), "/metadata")
swath_path(b::NisarBackend) = string(product_path(b), "/swaths")
frequency_path(b::NisarBackend) = string(swath_path(b), "/frequency", b.frequency)
identification_path(b::NisarBackend) = string(SCIENCE_ROOT, "/", b.band, "/identification")
orbit_path(b::NisarBackend) = string(metadata_path(b), "/orbit")

# A geocoded product stores its samples on a map grid under `grids` rather than `swaths`, so it has no
# slant-range axis for `RadarGeometry` to describe and is out of this package's scope. Recognized by name
# in order to say that, rather than failing on a missing group.
const GEOCODED_TYPES = ("GSLC", "GCOV", "GUNW", "GOFF")

"""
    nisar_file(path; transport = nothing) -> ChunkManifests.H5File

The product at `path` — a local file, an `https://` URL or an `s3://` URI — open for reading. Only the
metadata a read touches is transferred; a remote file is read in place, through `transport` when one is
given. Close it with `close`.
"""
function nisar_file(path::AbstractString; transport = nothing)
    t = transport !== nothing ? transport :
        startswith(path, "s3://") ? S3Transport(first(split(path[6:end], '/'; limit = 2))) :
        (startswith(path, "http://") || startswith(path, "https://")) ? HTTPTransport() :
        LocalTransport()
    return h5open_typed(path; transport = t)
end

"""
    nisar_band(f) -> String

The sensor band group of an open NISAR product: `"LSAR"` or `"SSAR"`.
"""
function nisar_band(f::H5File)::String
    h5_exists(f, SCIENCE_ROOT) || throw(ArgumentError(
        "`$(f.uri)` has no `/$SCIENCE_ROOT` group, so it is not a NISAR-format product"))
    for band in SENSOR_BANDS
        h5_exists(f, string(SCIENCE_ROOT, "/", band)) && return band
    end
    throw(ArgumentError("`" * f.uri * "` has none of " * join(SENSOR_BANDS, ", ") * " under `/" * SCIENCE_ROOT *
                        "`, so it is not a NISAR-format product"))
end

nisar_band(path::AbstractString) = _with_path(nisar_band, path)

"""
    nisar_product_type(f, band) -> String

The product type recorded in an open NISAR product's identification group, as the group name that holds
it. Early-mission products name the group `SLC` where later ones name it `RSLC`; both report `RSLC`.
"""
function nisar_product_type(f::H5File, band::AbstractString)::String
    root = string(SCIENCE_ROOT, "/", band)
    h5_exists(f, root * "/identification") || throw(ArgumentError(
        "`$(f.uri)` has no `identification` group under `/$root`"))
    declared = _text(f, root * "/identification/productType")
    # The declared type names the group except for the early-mission spelling.
    declared == "SLC" && return h5_exists(f, root * "/RSLC") ? "RSLC" : "SLC"
    return declared
end

nisar_product_type(path::AbstractString, band::AbstractString) =
    _with_path(f -> nisar_product_type(f, band), path)

function _with_path(g, path::AbstractString)
    f = nisar_file(path)
    try
        return g(f)
    finally
        close(f)
    end
end

# The accessors below read a dataset by its full path and convert to the type the layout calls for, so a
# product storing a count as `Int32` or times as `Float32` is read correctly here; libhdf5 converts as it
# reads. A dataset missing from the layout is an error naming the path.

function _require(f::H5File, name::AbstractString)
    h5_exists(f, name) || throw(ArgumentError(
        "`/$name` is not in `$(f.uri)`, so this is not a NISAR-format product"))
    h5_isgroup(f, name) && throw(ArgumentError(
        "`/$name` in `$(f.uri)` is a group where a dataset belongs, so this is not a NISAR-format product"))
    return String(name)
end

# A dataset's value where only a scalar is wanted: a zero-dimensional or one-element dataset.
function _scalar(f::H5File, name::AbstractString)::Float64
    v = h5_read(f, _require(f, name), Float64)
    length(v) == 1 || throw(ArgumentError(
        "`/$name` holds $(length(v)) values where one is expected, so this is not a NISAR-format product"))
    return v[1]
end

# A whole-number field, stored variously as an integer or a float across product generations.
function _scalar_int(f::H5File, name::AbstractString)
    v = _scalar(f, name)
    isinteger(v) || throw(ArgumentError(
        "`/$name` holds $v where a whole number is expected, so this is not a NISAR-format product"))
    return Int(v)
end

# A string dataset holding one string.
function _text(f::H5File, name::AbstractString)::String
    v = h5_read_strings(f, _require(f, name))
    length(v) == 1 || throw(ArgumentError(
        "`/$name` holds $(length(v)) strings where one is expected, so this is not a NISAR-format product"))
    return v[1]
end

# The first and last element of a one-dimensional dataset, and its length. Reading the dataset whole
# would transfer the entire axis — 27,000 azimuth times for a real RSLC — where only the ends are used.
function _axis_bounds(f::H5File, name::AbstractString)
    p = _require(f, name)
    dims = h5_dims(f, p)
    length(dims) == 1 || throw(ArgumentError("`/$name` is not one-dimensional"))
    n = dims[1]
    n > 0 || throw(ArgumentError("`/$name` is empty, so it has no bounds"))
    ends = h5_read_elements(f, p, Float64, [0, n - 1])
    return ends[1], ends[2], n
end

function _units_epoch(f::H5File, name::AbstractString)
    u = h5_attr_string(f, _require(f, name), "units")
    u === nothing && throw(ArgumentError(
        "`/$name` has no `units` attribute, so the epoch its times are measured from is unknown"))
    return parse_cf_epoch(u)
end

read_identification(b::NisarBackend) = _with_file(f -> read_identification(b, f), b)

function read_identification(b::NisarBackend, f::H5File)
    g = identification_path(b)
    # An absent optional field reads as empty rather than failing: a product missing one is still
    # usable, where one missing `absoluteOrbitNumber` is not.
    get_str(name)::String = h5_exists(f, g * "/" * name) ? _text(f, g * "/" * name) : ""
    return Identification(
        get_str("missionId"),
        get_str("productType"),
        _scalar_int(f, g * "/absoluteOrbitNumber"),
        get_str("orbitPassDirection"),
        get_str("lookDirection"),
        get_str("zeroDopplerStartTime"),
        get_str("zeroDopplerEndTime"),
        get_str("boundingPolygon"),
    )
end

read_geometry(b::NisarBackend) = _with_file(f -> read_geometry(b, f), b)

function read_geometry(b::NisarBackend, f::H5File)
    b.product_type in GEOCODED_TYPES && throw(ArgumentError(
        "`$(b.path)` is a $(b.product_type) product, which stores its samples on a map grid and so " *
        "carries no slant-range/azimuth geometry. This package reads SLCs in radar geometry; open an " *
        "RSLC, or use a raster library for a geocoded product"))
    freq = frequency_path(b)
    swaths = swath_path(b)

    near_range, far_range, nsamples = _axis_bounds(f, freq * "/slantRange")
    sensing_start, sensing_stop, nlines = _axis_bounds(f, swaths * "/zeroDopplerTime")

    look = _text(f, identification_path(b) * "/lookDirection")
    side = lowercase(look) == "left" ? LookLeft :
           lowercase(look) == "right" ? LookRight :
           throw(ArgumentError(
               "`$(b.path)` records lookDirection \"$look\"; expected \"Left\" or \"Right\""))

    return RadarGeometry(
        near_range,
        far_range,
        _scalar(f, freq * "/slantRangeSpacing"),
        SPEED_OF_LIGHT / _scalar(f, freq * "/processedCenterFrequency"),
        1 / _scalar(f, swaths * "/zeroDopplerTimeSpacing"),
        sensing_start,
        sensing_stop,
        nlines,
        nsamples,
        side,
        # The azimuth times and the orbit's state vector times share one epoch in this format, so the
        # epoch is read once here and reported for both.
        _units_epoch(f, swaths * "/zeroDopplerTime"),
    )
end

read_orbit(b::NisarBackend) = _with_file(f -> read_orbit(b, f), b)

function read_orbit(b::NisarBackend, f::H5File)
    g = orbit_path(b)
    time = h5_read(f, _require(f, g * "/time"), Float64)
    pos = _state_vectors(f, g * "/position")
    vel = _state_vectors(f, g * "/velocity")
    length(pos) == length(vel) || throw(DimensionMismatch(
        "orbit in `$(b.path)` has $(length(pos)) positions but $(length(vel)) velocities"))
    length(pos) == length(time) || throw(DimensionMismatch(
        "orbit in `$(b.path)` has $(length(time)) times but $(length(pos)) state vectors"))
    get_str(name, default)::String = h5_exists(f, g * "/" * name) ? _text(f, g * "/" * name) : default
    return StateVectors(
        time,
        pos,
        vel,
        _units_epoch(f, g * "/time"),
        get_str("interpMethod", "Hermite"),
        get_str("orbitType", "Custom"),
    )
end

# An (N, 3) dataset of state vectors, as N 3-vectors. HDF5 stores it row-major, so each vector's three
# components are adjacent.
function _state_vectors(f::H5File, name::AbstractString)
    p = _require(f, name)
    dims = h5_dims(f, p)
    (length(dims) == 2 && dims[2] == 3) || throw(ArgumentError(
        "`/$name` has dimensions $(join(dims, " x ")), expected N x 3"))
    v = h5_read(f, p, Float64)
    return [SVector{3,Float64}(v[3i - 2], v[3i - 1], v[3i]) for i in 1:dims[1]]
end

# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------

"""
    NisarRaster{T}

A NISAR product's sample array: its dataset as a ChunkManifests.jl typed array, permuted to
`(line, sample)`. The samples are read as `ComplexF32`.

**The file's layout is the transpose of Julia's.** A NISAR image is stored row-major with azimuth (or
northing) slowest, so the typed array is `(range, azimuth)`; this reverses that, so indexing is
`(line, sample)` for a radar-geometry product and `(northing, easting)` for a geocoded one — the order
every other array in this package uses.

Nothing is held: an RSLC band is 50511 x 57760 `ComplexF32` and a GSLC band 161280 x 80640, which are
23 GB and 104 GB of samples against files of 11 GiB and 10 GiB. A `DiskArrays` array, so a window reads
the chunks it touches, `eachchunk` reports the dataset's own storage chunks, and concurrent reads from
several tasks are safe. A chunk an unwritten region never allocated reads as the dataset's fill value.
"""
const NisarRaster{T} = DiskArrays.PermutedDiskArray{T,2,(2, 1),(2, 1),TypedArray{T}}

"""
    nisar_samples_path(b::NisarBackend, f, polarization = nothing) -> String

Where `b`'s polarization's samples live: under `swaths` for a radar-geometry product, `grids` for a
geocoded one.

The polarization is the first this frequency lists, which is the only one a single-polarization product
carries and the co-polarized channel of a dual one.
"""
function nisar_samples_path(b::NisarBackend, f::H5File, polarization = nothing)
    group = b.product_type in GEOCODED_TYPES ? "grids" : "swaths"
    base = string(product_path(b), "/", group, "/frequency", b.frequency)
    h5_isgroup(f, base) || throw(ArgumentError(
        "`/$base` is not a group in `$(f.uri)`, so this is not a NISAR-format product"))
    listing = base * "/listOfPolarizations"
    listed() = h5_read_strings(f, listing)
    pol = if polarization === nothing
        h5_exists(f, listing) || throw(ArgumentError(
            "`$(b.path)` lists no polarizations under $base"))
        first(listed())
    else
        String(polarization)
    end
    h5_exists(f, base * "/" * pol) || throw(ArgumentError(
        "`$(b.path)` has no $pol under $base; it lists " * join(listed(), ", ")))
    return string(base, "/", pol)
end

NisarRaster(b::NisarBackend, polarization = nothing) = _with_file(f -> NisarRaster(b, f, polarization), b)

NisarRaster(b::NisarBackend, f::H5File, polarization = nothing) =
    permutedims(h5_open_array(f, nisar_samples_path(b, f, polarization), ComplexF32), (2, 1))

read_pixels(b::NisarBackend) = NisarRaster(b)

# ---------------------------------------------------------------------------
# Geocoded products
# ---------------------------------------------------------------------------

"""
    GeocodedGrid

The map grid a geocoded product's samples lie on: `origin` at the outer corner of the first pixel,
signed `spacing`, `size` as `(rows, columns)`, and the EPSG code of the projection.

**`origin` is the corner, not the first pixel's center.** The product stores `xCoordinates` and
`yCoordinates` as centers; a geotransform names the corner, and a consumer intersecting two grids or
handing one to GDAL wants the latter.
"""
struct GeocodedGrid
    origin::Tuple{Float64,Float64}
    spacing::Tuple{Float64,Float64}
    size::Tuple{Int,Int}
    epsg::Int
end

"""
    GeocodedProduct <: AbstractSLC

A geocoded NISAR product: its identification, the map grid its samples lie on, and the samples.

Separate from [`SLC`](@ref) because it has no slant-range geometry to describe — `read_geometry` refuses
one — while still carrying an acquisition's identification and samples. [`pixels`](@ref) reads it.
"""
struct GeocodedProduct{B<:AbstractSLCBackend} <: AbstractSLC
    backend::B
    identification::Identification
    grid::GeocodedGrid
end

"""
    geocoded_grid(b::NisarBackend) -> GeocodedGrid

The map grid of a geocoded product's frequency group.
"""
geocoded_grid(b::NisarBackend) = _with_file(f -> geocoded_grid(b, f), b)

function geocoded_grid(b::NisarBackend, f::H5File)
    b.product_type in GEOCODED_TYPES || throw(ArgumentError(
        "`$(b.path)` is a $(b.product_type) product, which lies in radar geometry and so has no map " *
        "grid; `read_geometry` describes it instead"))
    g = string(product_path(b), "/grids/frequency", b.frequency)
    x0, _, nx = _axis_bounds(f, g * "/xCoordinates")
    y0, _, ny = _axis_bounds(f, g * "/yCoordinates")
    dx = _scalar(f, g * "/xCoordinateSpacing")
    dy = _scalar(f, g * "/yCoordinateSpacing")
    return GeocodedGrid((x0 - dx / 2, y0 - dy / 2), (dx, dy), (ny, nx), _scalar_int(f, g * "/projection"))
end

"""
    open_geocoded(src; frequency = nothing) -> GeocodedProduct

Open a geocoded NISAR product — a GSLC, GCOV, GUNW or GOFF — from a path, URL or
[`AbstractSLCSource`](@ref).

[`open_slc`](@ref) refuses these: they carry no slant-range axis, so there is no `RadarGeometry` to read.
What they do carry is a map grid and samples, which is what this returns.
"""
open_geocoded(path::AbstractString; frequency = nothing) = open_geocoded(source_for(path); frequency)

function open_geocoded(src::AbstractSLCSource; frequency = nothing)
    return _with_nisar(src; frequency) do b, f
        b.product_type in GEOCODED_TYPES || throw(ArgumentError(
            "`$(b.path)` is a $(b.product_type) product, which lies in radar geometry; use `open_slc`"))
        GeocodedProduct(b, read_identification(b, f), geocoded_grid(b, f))
    end
end

grid(g::GeocodedProduct) = g.grid
nlines(g::GeocodedProduct) = g.grid.size[1]
nsamples(g::GeocodedProduct) = g.grid.size[2]