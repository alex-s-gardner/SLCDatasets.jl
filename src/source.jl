# The single entry point that opens a product, and the dispatch that picks its reader.
#
# A source is separate from a backend: the backend knows a sensor's group and dataset names, the
# source knows how to get a readable file. Splitting them means a new sensor costs a backend and a new
# access route costs a source, rather than one costing both.

"""
    open_slc(src; frequency = nothing, orbit = nothing, swath = nothing, swaths = nothing,
             burst = nothing, polarization = nothing) -> SLC

Open a SAR acquisition.

`src` is an [`AbstractSLCSource`](@ref), or a path or URL, which is resolved to one. The sensor
backend is chosen by inspecting the product rather than by the caller naming it, so the same call
serves every supported format.

The identification and geometry records are read here; the state vectors are read on first access.

# NISAR

`frequency` selects the sub-band, defaulting to the first the product lists.

# Sentinel-1

`orbit` is the path to the POEORB or RESORB `.EOF` file holding the state vectors, and is required:
a Sentinel-1 product does not contain its own orbit.

A TOPS product is three subswaths of bursts rather than one image, so there are two geometries to
ask for. By default the result spans `swaths` — all three — the way hyp3-autorift's
`loadMetadataSlc` does. Passing `burst` instead describes that single burst of the single subswath
`swath`, matching `loadMetadata`. Use [`nbursts`](@ref) for the burst count.

`polarization` defaults to the first of `vv`, `vh`, `hh`, `hv` the product carries.

# Examples

```julia
s = open_slc("NISAR_L1_PR_RSLC_....h5")
s.geometry.prf

mosaic = open_slc("S1A_IW_SLC_....zip"; orbit = "S1A_OPER_AUX_POEORB_....EOF")
one    = open_slc("S1A_IW_SLC_....zip"; orbit = "...EOF", swath = 2, burst = 3)
```
"""
function open_slc(src::AbstractSLCSource; frequency = nothing, orbit = nothing,
                  swath = nothing, swaths = nothing, burst = nothing, polarization = nothing)
    path = localpath(src)
    ispath(path) || throw(ArgumentError("`$path` is not a readable file"))
    if is_safe_product(path)
        return _open_sentinel1(path; orbit, swath, swaths, burst, polarization)
    end
    isfile(path) || throw(ArgumentError("`$path` is not a readable file"))
    ishdf5(path) || throw(ArgumentError(
        "`$path` is neither an HDF5 file nor a Sentinel-1 SAFE product; those are the only " *
        "formats currently read"))
    return _open_nisar(src; frequency)
end

# A remote product is a NISAR one: a Sentinel-1 raster is only read from an unpacked `.SAFE`.
open_slc(src::Union{RemoteHTTP,RemoteS3}; frequency = nothing) = _open_nisar(src; frequency)

function _open_nisar(src::AbstractSLCSource; frequency)
    backend = _nisar_backend(src; frequency)
    return SLC(backend, read_identification(backend), read_geometry(backend))
end

# The one scan of a product that the band probe, the type dispatch, the frequency default and every
# metadata and sample read then share.
function _nisar_backend(src::AbstractSLCSource; frequency)
    path = source_uri(src)
    z = nisar_group(path; transport = source_transport(src))
    band = nisar_band(z)
    product_type = nisar_product_type(z, band)
    freq = frequency === nothing ? default_frequency(z, band) : String(frequency)
    return NisarBackend(path, band, product_type, freq, z)
end

function _open_sentinel1(path::AbstractString; orbit, swath, swaths, burst, polarization)
    # A single burst belongs to one subswath, so naming a burst without a subswath would silently pick
    # one; requiring `swath` alongside `burst` keeps the choice the caller's. Checked before the product
    # is parsed, since it is a fault in the call rather than in the product.
    if burst !== nothing
        swath === nothing && swaths === nothing && throw(ArgumentError(
            "`burst` selects a burst of one subswath, so it needs `swath = 1`, `2` or `3` as well"))
        named = swath !== nothing ? 1 : count(_ -> true, swaths)
        named == 1 || throw(ArgumentError(
            "a burst belongs to one subswath, but $named were named"))
    end
    p = _sentinel1_product(path; orbit, swath, swaths, polarization)
    backend = burst === nothing ? Sentinel1Backend(p) :
              Sentinel1Backend(p, only(p.swaths), Int(burst))
    return SLC(backend)
end

# A Sentinel-1 acquisition is refused without an orbit rather than returned half-read: the state vectors
# are in a separately distributed file, so there is nothing in the product to fall back on.
function _require_orbit(orbit, what::AbstractString)
    orbit === nothing && throw(ArgumentError(
        "$what is Sentinel-1, whose state vectors live in a separate POEORB/RESORB `.EOF` file " *
        "rather than in the product. Pass `orbit = \"<file>.EOF\"`"))
    orbit_path = String(orbit)
    isfile(orbit_path) || throw(ArgumentError("`$orbit_path` is not a readable orbit file"))
    return orbit_path
end

# The validation `open_slc` and `bursts` share, and the one parse of the container they both build on.
function _sentinel1_product(path::AbstractString; orbit, swath, swaths, polarization)
    orbit_path = _require_orbit(orbit, string("`", path, "`"))

    swath !== nothing && swaths !== nothing && throw(ArgumentError(
        "pass either `swath` for one subswath or `swaths` for several, not both"))

    # `Sentinel1Product` sorts, so the order here is the caller's and only the membership is checked.
    selected = swath !== nothing ? [Int(swath)] :
               swaths !== nothing ? collect(Int, swaths) : collect(S1_SUBSWATHS)
    isempty(selected) && throw(ArgumentError("`swaths` is empty; name at least one subswath"))
    all(in(S1_SUBSWATHS), selected) || throw(ArgumentError(
        "Sentinel-1 IW subswaths are 1, 2 and 3; got $(join(selected, ", "))"))
    allunique(selected) || throw(ArgumentError(
        "`swaths` repeats a subswath: $(join(selected, ", "))"))

    # The polarization is checked against what the product carries inside the constructor, which has the
    # container's listing already open.
    return Sentinel1Product(path; orbit = orbit_path, polarization, swaths = selected)
end

"""
    bursts(src; orbit, swath = 1, polarization = nothing) -> SLCSeries
    bursts(p::Sentinel1Product, swath = first(p.swaths)) -> SLCSeries

Every burst of one Sentinel-1 subswath, as a vector of [`SLC`](@ref)s.

The subswath's annotation is parsed once and shared, so indexing the result costs no further reading of
the product — unlike calling [`open_slc`](@ref) once per burst, which re-reads it each time.

# Examples

```julia
b = bursts("S1A_IW_SLC_....zip"; orbit = "...EOF", swath = 2)
length(b)
b[1].geometry.sensing_start
[nlines(s) for s in b]
```
"""
function bursts(src::AbstractSLCSource; orbit = nothing, swath::Integer = 1, polarization = nothing)
    path = localpath(src)
    is_safe_product(path) || throw(ArgumentError(
        "`$path` is not a Sentinel-1 SAFE product; bursts are a Sentinel-1 TOPS concept"))
    p = _sentinel1_product(path; orbit, swath, swaths = nothing, polarization)
    return bursts(p, swath)
end

bursts(path::AbstractString; kwargs...) = bursts(source_for(path); kwargs...)

function bursts(p::Sentinel1Product, swath::Integer = first(p.swaths))
    n = nbursts(p, swath)
    s = Int(swath)
    return SLCSeries(i -> SLC(Sentinel1Backend(p, s, i)), n, SLC{Sentinel1Backend})
end

open_slc(path::AbstractString; kwargs...) = open_slc(source_for(path); kwargs...)

"""
    source_for(spec) -> AbstractSLCSource

The source a path, URL or S3 URI denotes.
"""
function source_for(spec::AbstractString)
    startswith(spec, "s3://") && return RemoteS3(spec)
    (startswith(spec, "http://") || startswith(spec, "https://")) && return RemoteHTTP(spec)
    return LocalFile(spec)
end

"""
    default_frequency(z, band) -> String

The first sub-band a scanned NISAR product lists, as the one to read when the caller names none.
"""
function default_frequency(z::Zarr.ZGroup, band::AbstractString)::String
    name = string(SCIENCE_ROOT, "/", band, "/identification/listOfFrequencies")
    _haskey(z, name) || return "A"
    listed = _read(_dataset(z, name))
    isempty(listed) && throw(ArgumentError("the product lists no frequencies"))
    return _string(first(listed))
end
