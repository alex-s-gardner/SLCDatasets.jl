# The NISAR backend: HDF5 products under `/science/{LSAR,SSAR}/<productType>`.
#
# Group paths are built the way isce3's `nisar.products.readers.Base` builds them — a root path found
# by probing for the sensor band, then the product type, then a fixed subgroup name. The product type
# is read from the identification group and selects the reader, which is also how isce3's
# `open_product` dispatches. A product whose type is `SLC` is an early-mission RSLC.

const SCIENCE_ROOT = "science"
const SENSOR_BANDS = ("LSAR", "SSAR")

# Speed of light, for converting the processed center frequency to a wavelength. This is the value
# isce3 uses (`isce3.core.speed_of_light`), and the conversion must match it to the bit.
const SPEED_OF_LIGHT = 299792458.0

"""
    NisarBackend <: AbstractSLCBackend

A NISAR-format HDF5 product: a file, the sensor band group inside it, and the product type.

`frequency` selects the sub-band (`"A"` or `"B"`); a product lists the ones it carries at
`identification/listOfFrequencies`.
"""
struct NisarBackend <: AbstractSLCBackend
    path::String
    band::String
    product_type::String
    frequency::String
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
    nisar_band(h) -> String

The sensor band group of an open NISAR product: `"LSAR"` or `"SSAR"`.
"""
function nisar_band(h)::String
    haskey(h, SCIENCE_ROOT) || throw(ArgumentError(
        "`$(HDF5.filename(h))` has no `/$SCIENCE_ROOT` group, so it is not a NISAR-format product"))
    science = _group(h, SCIENCE_ROOT)
    for band in SENSOR_BANDS
        haskey(science, band) && return band
    end
    throw(ArgumentError("`$(HDF5.filename(h))` has none of $(join(SENSOR_BANDS, ", ")) under " *
                        "`/$SCIENCE_ROOT`, so it is not a NISAR-format product"))
end

nisar_band(path::AbstractString) = h5open(nisar_band, path, "r")

"""
    nisar_product_type(h, band) -> String

The product type recorded in an open NISAR product's identification group, as the group name that
holds it. Early-mission products name the group `SLC` where later ones name it `RSLC`; both report
`RSLC`.
"""
function nisar_product_type(h, band::AbstractString)::String
    g = _group(h, string(SCIENCE_ROOT, "/", band))
    haskey(g, "identification") || throw(ArgumentError(
        "`$(HDF5.filename(h))` has no `identification` group under `/$SCIENCE_ROOT/$band`"))
    declared = _text(g, "identification/productType")
    # The declared type names the group except for the early-mission spelling.
    declared == "SLC" && return haskey(g, "RSLC") ? "RSLC" : "SLC"
    return declared
end

nisar_product_type(path::AbstractString, band::AbstractString) =
    h5open(h -> nisar_product_type(h, band), path, "r")

# Indexing an HDF5 file or group gives back `Union{Attribute,Dataset,Datatype,Group}`, so `read` on the
# result infers as `Any` and every value derived from it dispatches at runtime. These accessors assert
# what the layout says a name is, which both keeps inference concrete and turns a product whose
# structure differs from the NISAR layout into an error naming the path rather than a `MethodError`
# deeper in.

@noinline _not_a(kind, parent, name) = throw(ArgumentError(
    "`$name` under `$(HDF5.name(parent))` is a $(nameof(typeof(parent))) member that is not a " *
    "$kind, so this is not a NISAR-format product"))

function _dataset(parent, name::AbstractString)
    d = parent[name]
    d isa HDF5.Dataset || _not_a("dataset", parent, name)
    return d
end

function _group(parent, name::AbstractString)
    g = parent[name]
    g isa HDF5.Group || _not_a("group", parent, name)
    return g
end

# `read` on a dataset infers as `Any`, because a dataset's element type is only known once the file is
# open. So each accessor below converts to the type the layout calls for and annotates that, which keeps
# inference concrete through the arithmetic that follows. The conversion is explicit rather than
# `read(d, T)`, which reinterprets the stored bytes as `T` and errors on a width mismatch instead of
# converting — a product storing a count as `Int32` or times as `Float32` is read correctly here.

# A dataset's value where only a scalar is wanted.
_scalar(parent, name::AbstractString)::Float64 = _only_number(_dataset(parent, name))

# A whole-number field, stored variously as an integer or a float across product generations.
function _scalar_int(parent, name::AbstractString)
    d = _dataset(parent, name)
    v = _only_number(d)
    isinteger(v) || throw(ArgumentError(
        "`$(HDF5.name(d))` holds $v where a whole number is expected, so this is not a " *
        "NISAR-format product"))
    return Int(v)
end

# A scalar dataset holds a number, a zero-dimensional array, or a one-element one depending on how the
# product was written; `only` covers the arrays and the bare number falls through.
_only_number(d)::Float64 = _as_number(read(d))
_as_number(x::Number) = Float64(x)
_as_number(x::AbstractArray) = Float64(only(x))

# A fixed-length string dataset, surfaced as a `String` or as bytes; `_string` narrows either and strips
# the NUL padding.
_text(parent, name::AbstractString)::String = _string(read(_dataset(parent, name)))

# The first and last element of a one-dimensional dataset, and its length. Reading the dataset whole
# would transfer the entire axis — 27,000 azimuth times for a real RSLC — where only the ends are used.
function _axis_bounds(parent, name::AbstractString)
    d = _dataset(parent, name)
    n = length(d)
    n > 0 || throw(ArgumentError("`$(HDF5.name(d))` is empty, so it has no bounds"))
    # Only the endpoints are read: the axis itself is tens of thousands of samples.
    return _as_number(d[1]), _as_number(d[n]), n
end

# `read_attribute` has no type-asserting form, so its `Any` is narrowed here rather than at each call.
function _units_epoch(parent, name::AbstractString)
    d = _dataset(parent, name)
    haskey(HDF5.attributes(d), "units") || throw(ArgumentError(
        "`$(HDF5.name(d))` has no `units` attribute, so the epoch its times are measured from is " *
        "unknown"))
    units = read_attribute(d, "units")
    return parse_cf_epoch(_string(units))
end

# Every read in the eager path goes through one open file handle. Opening a NISAR product costs about
# as much as the metadata reads themselves, so a reader that reopened per query would spend most of its
# time in `h5open`.
function read_identification(b::NisarBackend)
    return h5open(h -> read_identification(b, h), b.path, "r")
end

function read_identification(b::NisarBackend, h)
    g = _group(h, identification_path(b))
    # An absent optional field reads as empty rather than failing: a product missing one is still
    # usable, where one missing `absoluteOrbitNumber` is not.
    get_str(name)::String = haskey(g, name) ? _text(g, name) : ""
    return Identification(
        get_str("missionId"),
        get_str("productType"),
        _scalar_int(g, "absoluteOrbitNumber"),
        get_str("orbitPassDirection"),
        get_str("lookDirection"),
        get_str("zeroDopplerStartTime"),
        get_str("zeroDopplerEndTime"),
        get_str("boundingPolygon"),
    )
end

read_geometry(b::NisarBackend) = h5open(h -> read_geometry(b, h), b.path, "r")

function read_geometry(b::NisarBackend, h)
    b.product_type in GEOCODED_TYPES && throw(ArgumentError(
        "`$(b.path)` is a $(b.product_type) product, which stores its samples on a map grid and so " *
        "carries no slant-range/azimuth geometry. This package reads SLCs in radar geometry; open an " *
        "RSLC, or use a raster library for a geocoded product"))
    freq = _group(h, frequency_path(b))
    swaths = _group(h, swath_path(b))
    ident = _group(h, identification_path(b))

    near_range, far_range, nsamples = _axis_bounds(freq, "slantRange")
    sensing_start, sensing_stop, nlines = _axis_bounds(swaths, "zeroDopplerTime")

    look = _text(ident, "lookDirection")
    side = lowercase(look) == "left" ? LookLeft :
           lowercase(look) == "right" ? LookRight :
           throw(ArgumentError(
               "`$(b.path)` records lookDirection \"$look\"; expected \"Left\" or \"Right\""))

    return RadarGeometry(
        near_range,
        far_range,
        _scalar(freq, "slantRangeSpacing"),
        SPEED_OF_LIGHT / _scalar(freq, "processedCenterFrequency"),
        1 / _scalar(swaths, "zeroDopplerTimeSpacing"),
        sensing_start,
        sensing_stop,
        nlines,
        nsamples,
        side,
        # The azimuth times and the orbit's state vector times share one epoch in this format, so the
        # epoch is read once here and reported for both.
        _units_epoch(swaths, "zeroDopplerTime"),
    )
end

read_orbit(b::NisarBackend) = h5open(h -> read_orbit(b, h), b.path, "r")

function read_orbit(b::NisarBackend, h)
    g = _group(h, orbit_path(b))
    time = _read_vector(g, "time")
    # Stored (N, 3) row-major, which HDF5.jl presents as (3, N).
    pos = _read_matrix(g, "position")
    vel = _read_matrix(g, "velocity")
    size(pos, 1) == 3 || throw(ArgumentError(
        "orbit position in `$(b.path)` has leading dimension $(size(pos, 1)), expected 3"))
    axes(pos) == axes(vel) || throw(DimensionMismatch(
        "orbit position and velocity in `$(b.path)` have axes $(axes(pos)) and $(axes(vel))"))
    size(pos, 2) == length(time) || throw(DimensionMismatch(
        "orbit in `$(b.path)` has $(length(time)) times but $(size(pos, 2)) state vectors"))
    get_str(name, default)::String = haskey(g, name) ? _text(g, name) : default
    return StateVectors(
        time,
        _columns_as_svectors(pos),
        _columns_as_svectors(vel),
        _units_epoch(g, "time"),
        get_str("interpMethod", "Hermite"),
        get_str("orbitType", "Custom"),
    )
end

# A product may store these as `Float32`, so the element type is converted rather than propagated into
# `StateVectors`, whose fields are `Float64`. `convert` is a no-op when the file already holds `Float64`.
_read_vector(parent, name::AbstractString)::Vector{Float64} =
    convert(Vector{Float64}, read(_dataset(parent, name))::AbstractVector)

_read_matrix(parent, name::AbstractString)::Matrix{Float64} =
    convert(Matrix{Float64}, read(_dataset(parent, name))::AbstractMatrix)

# The (3, N) array HDF5.jl hands back, as N 3-vectors.
function _columns_as_svectors(a::AbstractMatrix{Float64})
    Base.require_one_based_indexing(a)
    return [SVector{3,Float64}(a[1, i], a[2, i], a[3, i]) for i in axes(a, 2)]
end

# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------

"""
    NisarRaster{T} <: AbstractMatrix{T}

A NISAR product's sample array, read a window at a time out of its HDF5 dataset.

**The file's layout is the transpose of Julia's.** A NISAR image is stored row-major with azimuth (or
northing) slowest, so `HDF5.jl` reports it as `(range, azimuth)`; this reverses that, so indexing is
`(line, sample)` for a radar-geometry product and `(northing, easting)` for a geocoded one — the order
every other array in this package uses.

Nothing is held: an RSLC band is 50511 x 57760 `ComplexF32` and a GSLC band 161280 x 80640, which are
23 GB and 104 GB of samples against files of 11 GiB and 10 GiB. Index with ranges, so that a window is
one hyperslab read rather than one read per element.
"""
struct NisarRaster{T} <: AbstractMatrix{T}
    path::String
    dataset::String
    dims::Tuple{Int,Int}
end

# **HDF5 is not thread-safe unless the library was built for it, and the shipped one is not.** Two tasks
# reading one file segfault inside the C library rather than racing visibly. A read is milliseconds, so
# serializing them costs far less than the work a caller does per window.
const HDF5_IO = ReentrantLock()

"""
    nisar_samples_path(b::NisarBackend) -> String

Where `b`'s polarization's samples live: under `swaths` for a radar-geometry product, `grids` for a
geocoded one.

The polarization is the first this frequency lists, which is the only one a single-polarization product
carries and the co-polarized channel of a dual one.
"""
function nisar_samples_path(b::NisarBackend, polarization = nothing)
    group = b.product_type in GEOCODED_TYPES ? "grids" : "swaths"
    base = string(product_path(b), "/", group, "/frequency", b.frequency)
    pol = polarization === nothing ? nothing : String(polarization)
    return h5open(b.path, "r") do h
        g = _group(h, base)
        if pol === nothing
            haskey(g, "listOfPolarizations") || throw(ArgumentError(
                "`$(b.path)` lists no polarizations under $base"))
            pol = String(first(read(g["listOfPolarizations"])))
        end
        haskey(g, pol) || throw(ArgumentError(
            "`$(b.path)` has no $pol under $base; it lists " *
            join(String.(read(g["listOfPolarizations"])), ", ")))
        return string(base, "/", pol)
    end
end

function NisarRaster(b::NisarBackend, polarization = nothing)
    ds = nisar_samples_path(b, polarization)
    T, nx, ny = h5open(b.path, "r") do h
        d = h[ds]
        (eltype(d), size(d, 1), size(d, 2))
    end
    return NisarRaster{T}(b.path, ds, (ny, nx))
end

Base.size(r::NisarRaster) = r.dims

function Base.getindex(r::NisarRaster{T}, rows::AbstractUnitRange{<:Integer},
                       cols::AbstractUnitRange{<:Integer}) where {T}
    @boundscheck checkbounds(r, rows, cols)
    return @lock HDF5_IO h5open(h -> permutedims(h[r.dataset][cols, rows]), r.path, "r")
end

Base.getindex(r::NisarRaster, i::Int, j::Int) = r[i:i, j:j][1, 1]
Base.getindex(r::NisarRaster, rows::AbstractUnitRange{<:Integer}, j::Int) = r[rows, j:j][:, 1]
Base.getindex(r::NisarRaster, i::Int, cols::AbstractUnitRange{<:Integer}) = r[i:i, cols][1, :]

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
function geocoded_grid(b::NisarBackend)
    b.product_type in GEOCODED_TYPES || throw(ArgumentError(
        "`$(b.path)` is a $(b.product_type) product, which lies in radar geometry and so has no map " *
        "grid; `read_geometry` describes it instead"))
    base = string(product_path(b), "/grids/frequency", b.frequency)
    return h5open(b.path, "r") do h
        g = _group(h, base)
        x = read(_dataset(g, "xCoordinates"))
        y = read(_dataset(g, "yCoordinates"))
        dx = _scalar(g, "xCoordinateSpacing")
        dy = _scalar(g, "yCoordinateSpacing")
        return GeocodedGrid((first(x) - dx / 2, first(y) - dy / 2), (dx, dy),
                            (length(y), length(x)), _scalar_int(g, "projection"))
    end
end

"""
    open_geocoded(path; frequency = nothing) -> GeocodedProduct

Open a geocoded NISAR product — a GSLC, GCOV, GUNW or GOFF.

[`open_slc`](@ref) refuses these: they carry no slant-range axis, so there is no `RadarGeometry` to read.
What they do carry is a map grid and samples, which is what this returns.
"""
function open_geocoded(path::AbstractString; frequency = nothing)
    ispath(path) || throw(ArgumentError("`$path` is not a readable file"))
    ishdf5(path) || throw(ArgumentError("`$path` is not an HDF5 file"))
    return h5open(path, "r") do h
        band = nisar_band(h)
        product_type = nisar_product_type(h, band)
        product_type in GEOCODED_TYPES || throw(ArgumentError(
            "`$path` is a $product_type product, which lies in radar geometry; use `open_slc`"))
        freq = frequency === nothing ? default_frequency(h, band) : String(frequency)
        b = NisarBackend(path, band, product_type, freq)
        return GeocodedProduct(b, read_identification(b, h), geocoded_grid(b))
    end
end

grid(g::GeocodedProduct) = g.grid
nlines(g::GeocodedProduct) = g.grid.size[1]
nsamples(g::GeocodedProduct) = g.grid.size[2]
