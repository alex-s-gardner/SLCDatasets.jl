using SLCDatasets: open_measurement
import DiskArrays

@testset "reads a striped raster" begin
    mktempdir() do dir
        pixels = tiff_pattern(9, 13)
        t = open_measurement(write_tiff(joinpath(dir, "p.tiff"), pixels))

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
            t = open_measurement(write_tiff(joinpath(dir, "pad$pad.tiff"), pixels; pad))
            @test t == pixels
            @test t[2:4, 3:8] == pixels[2:4, 3:8]
            @test t[3, 5] == pixels[3, 5]
        end
    end
end

@testset "indexing outside the raster is caught" begin
    mktempdir() do dir
        t = open_measurement(write_tiff(joinpath(dir, "p.tiff"), tiff_pattern(4, 5)))
        @test_throws BoundsError t[5, 1]
        @test_throws BoundsError t[1, 6]
        @test_throws BoundsError t[0, 1]
        @test_throws BoundsError t[1:5, 1:5]
        @test_throws BoundsError t[1:4, 1:6]
    end
end

# Each of these is a layout whose samples are not what a measurement raster holds, or whose bytes the
# reader cannot serve as they lie. Reading one anyway returns plausible garbage, so each is refused by name.
@testset "a raster laid out otherwise is refused" begin
    mktempdir() do dir
        pixels = tiff_pattern(4, 5)

        @test_throws "LZW" open_measurement(
            write_tiff(joinpath(dir, "lzw.tiff"), pixels; compression = 5))

        @test_throws "opposite byte order" open_measurement(
            write_tiff(joinpath(dir, "be.tiff"), pixels; bigendian = true))

        @test_throws "holds Float32 samples" open_measurement(
            write_tiff(joinpath(dir, "float.tiff"), pixels; sample_format = 3))

        @test_throws "holds Complex{Int8} samples" open_measurement(
            write_tiff(joinpath(dir, "bits.tiff"), pixels; bits_per_sample = 16))

        @test_throws "one sample per pixel" open_measurement(
            write_tiff(joinpath(dir, "bands.tiff"), pixels; samples_per_pixel = 3))

        # A strip table shorter than the raster is tall leaves lines unaddressed.
        @test_throws "strip offsets" open_measurement(
            write_tiff(joinpath(dir, "short.tiff"), pixels; nstrips = 2))
    end
end

@testset "a file that is not a TIFF is refused" begin
    mktempdir() do dir
        p = joinpath(dir, "not.tiff")
        write(p, "this is not a TIFF file at all, just some bytes")
        @test_throws Exception open_measurement(p)
        @test_throws "no such file" open_measurement(joinpath(dir, "absent.tiff"))
    end
end

# A real measurement raster: uncompressed, one strip per line. Runs only where a granule is available.
if haskey(ENV, "SLCDATASETS_S1_TIFF")
    @testset "a real measurement raster reads" begin
        t = open_measurement(ENV["SLCDATASETS_S1_TIFF"])
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
        t = open_measurement(write_tiff(joinpath(dir, "da.tiff"), pixels))
        @test t isa DiskArrays.AbstractDiskArray

        # Chunks are bands of whole lines.
        @test DiskArrays.haschunks(t) isa DiskArrays.Chunked
        @test map(length, first(DiskArrays.eachchunk(t)))[2] == 23

        # Everything `DiskArrays` builds on `readblock!`, none of it written here.
        @test t[:, :] == pixels
        @test t[5:9, 3:7] == pixels[5:9, 3:7]
        @test t[1:2:20, 2:3:20] == pixels[1:2:20, 2:3:20]
        @test t[CartesianIndex(6, 4)] === pixels[6, 4]
        @test view(t, 3:6, 2:4)[:, :] == pixels[3:6, 2:4]
        @test t[11, :] == pixels[11, :]
        @test t[:, 12] == pixels[:, 12]
        @test_throws BoundsError t[1:38, 1:23]
        @test t[6, 4] === pixels[6, 4]

        # `readblock!` into a view, which is what `DiskArrays` hands it when batching.
        dest = fill(Complex{Int16}(0, 0), 8, 8)
        DiskArrays.readblock!(t, view(dest, 2:5, 3:5), 10:13, 7:9)
        @test dest[2:5, 3:5] == pixels[10:13, 7:9]
        @test all(iszero, dest[1, :])

        # Concurrent windows, as a blocked reader takes them.
        wins = [(i:(i + 9), 2:20) for i in 1:8]
        @test fetch.([Threads.@spawn t[r, c] for (r, c) in wins]) == [pixels[r, c] for (r, c) in wins]
    end
end
