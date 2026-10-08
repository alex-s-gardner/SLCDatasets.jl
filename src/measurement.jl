# A Sentinel-1 measurement raster: one subswath's bursts stacked in an uncompressed TIFF of complex 16-bit
# samples, one strip per line. ChunkManifests.jl scans it into a Zarr array whose chunks are bands of whole
# lines, located from the strip offsets without reading a sample.

"""
    open_measurement(path) -> AbstractMatrix{Complex{Int16}}

The measurement raster at `path`, as a lazy `DiskArrays` array indexed `(line, sample)`.

The file stores lines one after another, which Zarr.jl reports as `(sample, line)`; the result is that
array permuted. A window reads the bands of lines it touches, and concurrent reads from several tasks are
safe.
"""
function open_measurement(path::AbstractString)
    a = scan(path; driver = GeoTIFFDriver(), level = 0)["0"]["data"]
    ndims(a) == 2 || throw(ArgumentError(
        "`$path` holds $(size(a, 1)) samples per pixel; a Sentinel-1 measurement raster holds one " *
        "sample per pixel"))
    eltype(a) === Complex{Int16} || throw(ArgumentError(
        "`$path` holds $(eltype(a)) samples; a Sentinel-1 measurement raster holds complex 16-bit " *
        "integers"))
    return permutedims(a, (2, 1))
end
