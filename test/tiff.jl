using SLCDatasets: open_tiff, StripedTiff
import DiskArrays

@testset "reads a striped raster" begin
    mktempdir() do dir
        pixels = tiff_pattern(9, 13)
        t = open_tiff(write_tiff(joinpath(dir, "p.tiff"), pixels))

        @test size(t) == (9, 13)
        @test eltype(t) == Complex{Int16}
        @test t == pixels

        # Scalar and window indexing must agree, since they take different paths to the same bytes.
        @test t[3, 4] == pixels[3, 4]
        @test t[2:5, 3:9] == pixels[2:5, 3:9]
        @test t[1:1, 1:1] == pixels[1:1, 1:1]
        @test t[9:9, 13:13] == pixels[9:9, 13:13]
        @test collect(t) == pixels
    end
end

# A real writer puts the raster wherever the tags happen to land, which is not a multiple of the sample
# size. A reader loading a sample as a 32-bit word rather than assembling its bytes fails only here.
@testset "reads a raster at an unaligned offset" begin
    mktempdir() do dir
        pixels = tiff_pattern(5, 9)
        for pad in 1:3
            t = open_tiff(write_tiff(joinpath(dir, "pad$pad.tiff"), pixels; pad))
            @test t == pixels
            @test t[2:4, 3:8] == pixels[2:4, 3:8]
            @test t[3, 5] == pixels[3, 5]
        end
    end
end

@testset "reads a big-endian raster" begin
    mktempdir() do dir
        pixels = tiff_pattern(6, 7)
        t = open_tiff(write_tiff(joinpath(dir, "be.tiff"), pixels; bigendian = true))
        @test size(t) == (6, 7)
        @test t == pixels
        @test t[2:4, 2:6] == pixels[2:4, 2:6]
    end
end

@testset "indexing outside the raster is caught" begin
    mktempdir() do dir
        t = open_tiff(write_tiff(joinpath(dir, "p.tiff"), tiff_pattern(4, 5)))
        @test_throws BoundsError t[5, 1]
        @test_throws BoundsError t[1, 6]
        @test_throws BoundsError t[0, 1]
        @test_throws BoundsError t[1:5, 1:5]
        @test_throws BoundsError t[1:4, 1:6]
    end
end

# Each of these is a layout whose per-line offsets either do not exist or do not mean what the reader
# would take them to mean. Reading one anyway returns plausible garbage, so each must be refused, and
# the message must name the property rather than the byte that failed.
@testset "a raster laid out otherwise is refused" begin
    mktempdir() do dir
        pixels = tiff_pattern(4, 5)

        @test_throws "is compressed" open_tiff(
            write_tiff(joinpath(dir, "lzw.tiff"), pixels; compression = 5))

        @test_throws "one line per strip" open_tiff(
            write_tiff(joinpath(dir, "rows.tiff"), pixels; rows_per_strip = 2))

        @test_throws "tiles rather than strips" open_tiff(
            write_tiff(joinpath(dir, "tiled.tiff"), pixels; tiled = true))

        @test_throws "BigTIFF" open_tiff(
            write_tiff(joinpath(dir, "big.tiff"), pixels; magic = 43))

        @test_throws "one sample per pixel" open_tiff(
            write_tiff(joinpath(dir, "bands.tiff"), pixels; samples_per_pixel = 3))

        @test_throws "complex-integer samples" open_tiff(
            write_tiff(joinpath(dir, "float.tiff"), pixels; sample_format = 3))

        @test_throws "32-bit samples" open_tiff(
            write_tiff(joinpath(dir, "bits.tiff"), pixels; bits_per_sample = 16))

        # A strip table shorter than the raster is tall leaves lines unaddressed.
        @test_throws "one offset per line" open_tiff(
            write_tiff(joinpath(dir, "short.tiff"), pixels; nstrips = 2))

        # A byte count that is not `4 * nsamples` means a line is not the run the tags imply, so
        # indexing would read across line boundaries.
        @test_throws "occupy" open_tiff(
            write_tiff(joinpath(dir, "bytes.tiff"), pixels; strip_bytes = 12))
    end
end

@testset "a file that is not a TIFF is refused" begin
    mktempdir() do dir
        p = joinpath(dir, "not.tiff")
        write(p, "this is not a TIFF file at all, just some bytes")
        @test_throws "byte-order mark" open_tiff(p)

        write(joinpath(dir, "tiny.tiff"), UInt8[0x49, 0x49, 0x2a])
        @test_throws "too short" open_tiff(joinpath(dir, "tiny.tiff"))

        @test_throws "not a file" open_tiff(joinpath(dir, "absent.tiff"))
    end
end

# The one property the whole design rests on: a real measurement raster is uncompressed with one strip
# per line, so its lines are addressable without decoding. Runs only where a granule is available.
if haskey(ENV, "SLCDATASETS_S1_TIFF")
    @testset "a real measurement raster has the layout the reader assumes" begin
        t = open_tiff(ENV["SLCDATASETS_S1_TIFF"])
        @test size(t, 1) > 0
        @test size(t, 2) > 0
        # Window and scalar reads agree deep inside a multi-gigabyte file.
        i, j = size(t, 1) ÷ 2, size(t, 2) ÷ 2
        @test t[i:(i + 3), j:(j + 4)] == [t[a, b] for a in i:(i + 3), b in j:(j + 4)]
    end
end

@testset "the raster is a DiskArrays array" begin
    pixels = tiff_pattern(37, 23)
    mktempdir() do dir
        t = open_tiff(write_tiff(joinpath(dir, "da.tiff"), pixels))
        @test t isa DiskArrays.AbstractDiskArray

        # The strips are the file's granularity and are reported as such: one line each.
        @test DiskArrays.haschunks(t) isa DiskArrays.Chunked
        @test map(length, first(DiskArrays.eachchunk(t))) == (1, 23)
        @test size(DiskArrays.eachchunk(t)) == (37, 1)

        # Everything `DiskArrays` builds on `readblock!`, none of it written here.
        @test t[:, :] == pixels
        @test t[5:9, 3:7] == pixels[5:9, 3:7]
        @test t[1:2:20, 2:3:20] == pixels[1:2:20, 2:3:20]
        @test t[CartesianIndex(6, 4)] === pixels[6, 4]
        @test view(t, 3:6, 2:4)[:, :] == pixels[3:6, 2:4]
        @test t[11, :] == pixels[11, :]
        @test t[:, 12] == pixels[:, 12]
        @test_throws BoundsError t[1:38, 1:23]

        # The scalar method is kept rather than left to `DiskArrays`, which would route one element through
        # `readblock!`: on a real subswath that is the difference between 4.9 ns an element and a window read.
        @test t[6, 4] === pixels[6, 4]
        @test @inferred(t[6, 4]) isa Complex{Int16}

        # `readblock!` into a view, which is what `DiskArrays` hands it when batching.
        dest = fill(Complex{Int16}(0, 0), 8, 8)
        DiskArrays.readblock!(t, view(dest, 2:5, 3:5), 10:13, 7:9)
        @test dest[2:5, 3:5] == pixels[10:13, 7:9]
        @test all(iszero, dest[1, :])
    end
end
