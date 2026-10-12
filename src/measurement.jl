# A Sentinel-1 measurement raster: one subswath's bursts stacked in an uncompressed TIFF of complex 16-bit
# samples, one strip per line. ChunkManifests.jl opens it as a typed array whose chunks are its strips,
# located from the strip offsets without reading a sample.

"""
    open_measurement(path) -> AbstractMatrix{Complex{Int16}}

The measurement raster at `path`, as a lazy `DiskArrays` array indexed `(line, sample)`.

The file stores lines one after another, which the typed array reports as `(sample, line)`; the result
is that array permuted. A window reads the lines it touches, and concurrent reads from several tasks are
safe. A file with more than one band, or samples other than complex 16-bit integers, is refused.
"""
open_measurement(path::AbstractString) =
    permutedims(open_tiff(path, Complex{Int16}; transport = LocalTransport()), (2, 1))::MeasurementRaster

"""
    MeasurementRaster

The type [`open_measurement`](@ref) returns.
"""
const MeasurementRaster = DiskArrays.PermutedDiskArray{Complex{Int16},2,(2, 1),(2, 1),TypedArray{Complex{Int16}}}
