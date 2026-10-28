"""
    SLCDatasets

Read single-look complex SAR products into one type, whatever the sensor.

[`open_slc`](@ref) returns an [`SLC`](@ref) carrying an [`Identification`](@ref) and a
[`RadarGeometry`](@ref) — the slant-range/azimuth geometry in the units a geometry kernel wants — plus
a [`StateVectors`](@ref) record, read on first access.

The scope is SLCs in radar geometry, not SAR products generally. An interferogram or a covariance
product carries a multilooked grid, and a geocoded SLC carries map coordinates rather than a slant-range
axis — neither is described by [`RadarGeometry`](@ref), so neither is read here.

NISAR-format HDF5 range-Doppler products (RSLC) are read through libhdf5 for their metadata and
ChunkManifests.jl's typed arrays for their samples. The metadata of one is a few tens of kilobytes near the start of a file tens of gigabytes
long, so a remote product is opened, and its samples read a chunk at a time, without transferring it.

Sentinel-1 IW SLCs are read from a `.SAFE` directory or the zip of one, taking geometry from the
annotation XML and state vectors from a POEORB or RESORB `.EOF` file, which a Sentinel-1 product does
not carry. A geometry costs the annotation alone, so a zip need not be unpacked to read one. A TOPS
product is three subswaths of bursts rather than one image; [`open_slc`](@ref) describes either the
mosaic across them or one individual burst.

An annotation document is a couple of megabytes holding a few dozen scalars, so a product's subswaths
are parsed once into a [`Sentinel1Product`](@ref) and every geometry, burst and orbit derived from it
reads nothing further. [`bursts`](@ref) returns a whole subswath's bursts as an
[`SLCSeries`](@ref) — an `AbstractVector` of [`SLC`](@ref)s over that one parse, rather than one
`open_slc` and one parse per burst.

[`merge_bursts`](@ref) folds those bursts into one acquisition, placing them on a single uniform azimuth
grid so that a line's time is `sensing_start + line / prf` throughout. [`pixels`](@ref) and
[`amplitude`](@ref) read its samples as they are indexed, and [`validmask`](@ref) says which of them
were imaged. A zip stores each raster deflated, so a line is not addressable without inflating everything
before it: a zipped product's rasters are inflated whole into a directory when its samples are first read.
"""
module SLCDatasets

import Dates
import DiskArrays
using ChunkManifests: H5File, HTTPTransport, LeafTransport, LocalTransport, S3Transport, TypedArray,
                      h5_attr_string, h5_dims, h5_exists, h5_isgroup, h5_open_array, h5_read, h5_read_elements, h5_read_strings, h5open_typed,
                      http_request, http_resolve, fetchall, fetchranges, objectsize, ByteRange,
                      open_tiff, ZipMember, zip_members, read_member, extract_member, inflate_ranges!
using Dates: DateTime
import EzXML
using EzXML: parsexml, readxml, root, findfirst, findall, nodecontent, nodename, eachelement
using StaticArrays: SVector
import Zarr

export open_slc, bursts, orbit, nlines, nsamples, start_datetime, stop_datetime
export repeat_interval, epoch_offset, nbursts
export merge_bursts, pixels, amplitude, validmask
# Whether the samples carry a TOPS azimuth ramp, which decides whether they may be interpolated. Exported
# because a consumer of `pixels` has to ask before resampling them; `deramp_parameters` supplies what
# removing the ramp needs, and `burst_at` says which burst's ramp a line of a merged image carries.
export is_tops, deramp_parameters, burst_at
export DerampParameters, RangePolynomial, nearest_polynomial
export SLC, SLCSeries, Identification, RadarGeometry, StateVectors
export Sentinel1Product, subswath_borders, prefetch_rasters
export LocalFile, RemoteHTTP, RemoteS3, earthdata_token
export AsfBurst, asf_bursts, asf_safe_name, is_asf_burst_url
export ZipBurst, zip_bursts, zip_burst_spec, is_zip_burst_spec, burst_footprints
# `grid` stays unexported for the reason `LookSide` does: a raster or geometry package consuming this one
# defines its own, and exporting both makes the name ambiguous at every call site. `GeocodedProduct`
# carries it as a field, and `geocoded_grid` reads one from a backend.
export open_geocoded, GeocodedProduct, GeocodedGrid, geocoded_grid

# `LookSide`, `LookLeft` and `LookRight` are deliberately not exported: a geometry package consuming
# this one defines its own, and exporting both makes the name ambiguous at every call site.

include("util.jl")
include("time.jl")
include("types.jl")
include("measurement.jl")
include("zipped.jl")
include("nisar.jl")
# `extractor.jl` and `zipburst.jl` before `sentinel1.jl`: a product assembled from bursts holds them.
include("extractor.jl")
include("zipburst.jl")
include("sentinel1.jl")
include("burstgrid.jl")
include("concatenated.jl")
include("pixels.jl")
# `asf.jl` before `merge.jl`: a merge dispatches on how its bursts were delivered.
include("asf.jl")
include("merge.jl")
# `remote.jl` before `source.jl`: the remote sources are what `open_slc` dispatches on.
include("remote.jl")
include("source.jl")
include("pairing.jl")

end
