# Bursts read out of an SLC zip.
#
# Offline, a zip is written here holding the committed fixture's annotations, cut to a few short bursts,
# beside measurement rasters laid out as the processor writes them — a copy of the directory at byte 8
# with its strip table zeroed, the samples, then the directory the header points to — and each burst read
# out of it must be the rows of its subswath's raster. Comparing with ASF's extractor needs a zip and the
# extractor's rasters of the same bursts, so it runs only where they are named: set
# `SLCDATASETS_ZIP_BURST_ZIP` to the zip (a path or URL), `SLCDATASETS_ZIP_BURST_URLS` to some of its
# bursts' extractor URLs (comma-separated) and `SLCDATASETS_ZIP_BURST_ORBIT` to its orbit file.

using SLCDatasets
using SLCDatasets: SubswathAnnotation, annotation_xml, open_measurement, zip_rasters, burst_files
import Zarr, ZarrCore
using Test

# A zip of `members` (name => bytes), each deflated, written field by field.
function deflated_zip(path::AbstractString, members::Vector{Pair{String,Vector{UInt8}}})
    io = IOBuffer()
    central = IOBuffer()
    le16(x) = htol(UInt16(x))
    le32(x) = htol(UInt32(x))
    for (name, bytes) in members
        # Zarr's zlib stream less its two-byte header and four-byte checksum is the raw deflate a zip holds.
        data = ZarrCore.zcompress(bytes, Zarr.ZlibCompressor())[3:(end - 4)]
        offset = position(io)
        write(io, le32(0x04034b50), le16(20), le16(0), le16(8), le16(0), le16(0), le32(0),
              le32(length(data)), le32(length(bytes)), le16(ncodeunits(name)), le16(0), codeunits(name), data)
        write(central, le32(0x02014b50), le16(20), le16(20), le16(0), le16(8), le16(0), le16(0), le32(0),
              le32(length(data)), le32(length(bytes)), le16(ncodeunits(name)), le16(0), le16(0), le16(0),
              le16(0), le32(0), le32(offset), codeunits(name))
    end
    cd = take!(central)
    cd_offset = position(io)
    write(io, cd, le32(0x06054b50), le16(0), le16(0), le16(length(members)), le16(length(members)),
          le32(length(cd)), le32(cd_offset), le16(0))
    write(path, take!(io))
    return path
end

# A measurement raster as the processor writes it: a copy of the directory at byte 8 and the values it
# points to, its strip table zeroed when `zeroed`, then the samples, `pad` bytes, and the directory the
# header points to.
function processor_tiff(pixels::Matrix{Complex{Int16}}; zeroed::Bool = true, pad::Int = 0)
    nl, ns = size(pixels)
    entries(offsets_at, counts_at) = [(256, 3, 1, ns), (257, 3, 1, nl), (258, 3, 1, 32), (259, 3, 1, 1),
                                      (262, 3, 1, 1), (273, 4, nl, offsets_at), (277, 3, 1, 1),
                                      (278, 3, 1, 1), (279, 4, nl, counts_at), (284, 3, 1, 1), (339, 3, 1, 5)]
    dirbytes = 2 + 12 * 11 + 4
    copy_offsets = 8 + dirbytes
    copy_counts = copy_offsets + 4nl
    data = copy_counts + 4nl
    tail = data + 4 * nl * ns + pad
    tail_offsets = tail + dirbytes
    tail_counts = tail_offsets + 4nl
    io = IOBuffer()
    w16(x) = write(io, htol(UInt16(x)))
    w32(x) = write(io, htol(UInt32(x)))
    writedir(es) = (w16(length(es)); for (t, ty, c, v) in es; w16(t); w16(ty); w32(c); w32(v); end; w32(0))
    write(io, UInt8('I'), UInt8('I')); w16(42); w32(tail)
    writedir(entries(copy_offsets, copy_counts))
    for i in 1:nl; w32(zeroed ? 0 : data + 4ns * (i - 1)); end
    for i in 1:nl; w32(zeroed ? 0 : 4ns); end
    for i in 1:nl, j in 1:ns
        write(io, htol(real(pixels[i, j])), htol(imag(pixels[i, j])))
    end
    write(io, zeros(UInt8, pad))
    writedir(entries(tail_offsets, tail_counts))
    for i in 1:nl; w32(data + 4ns * (i - 1)); end
    for i in 1:nl; w32(4ns); end
    return take!(io)
end

# An annotation document cut to `n` bursts of `S` samples, every burst's valid region the whole burst but
# its first and last lines and samples. The line count is the document's own: with the line interval it
# sets how far apart consecutive bursts may start.
function short_annotation(doc::String, n::Int, S::Int)
    L = parse(Int, match(r"<linesPerBurst>(\d+)</linesPerBurst>", doc)[1])
    m = match(r"<burstList count=\"\d+\">(.*)</burstList>"s, doc)
    items = [x.match for x in eachmatch(r"<burst>.*?</burst>"s, m[1])][1:n]
    valid(lo, hi) = join((2 <= i <= L - 1 ? lo : -1 for i in 1:L), " ")
    items = [replace(x, r"<firstValidSample>.*?</firstValidSample>"s => "<firstValidSample>$(valid(1, S - 2))</firstValidSample>",
                     r"<lastValidSample>.*?</lastValidSample>"s => "<lastValidSample>$(valid(S - 2, S - 2))</lastValidSample>")
             for x in items]
    doc = replace(doc, m.match => "<burstList count=\"$n\">" * join(items) * "</burstList>")
    doc = replace(doc, r"<samplesPerBurst>\d+</samplesPerBurst>" => "<samplesPerBurst>$S</samplesPerBurst>")
    return doc
end

if !isempty(S1_PRODUCTS)
    @testset "bursts read out of an SLC zip" begin
        p = first(S1_PRODUCTS)
        pol = String(p.gold.polarization)
        mission = Sentinel1Product(p.safe; orbit = p.eof, polarization = pol).annotations[1].mission
        name = "$(mission)_IW_SLC__1SDV_20000101T000000_20000101T000030_000001_000001_AAAA"
        n, S = 4, 9
        mktempdir() do dir
            rasters = Dict{Int,Matrix{Complex{Int16}}}()
            members = Pair{String,Vector{UInt8}}[]
            L = Dict{Int,Int}()
            for sw in 1:3
                stem = "$(lowercase(mission))-iw$sw-slc-$pol-20000101t000000-20000101t000030-000001-000001-00$sw"
                doc = short_annotation(annotation_xml(p.safe, sw, pol), n, S + sw)
                L[sw] = parse(Int, match(r"<linesPerBurst>(\d+)</linesPerBurst>", doc)[1])
                rasters[sw] = Complex{Int16}[complex(Int16((10sw + i + j) % 30000), Int16(-i - j))
                                             for i in 1:(n * L[sw]), j in 1:(S + sw)]
                push!(members, "$name.SAFE/annotation/$stem.xml" => Vector{UInt8}(doc))
                # The zeroed strip table, the set one, and a byte of padding before the directory.
                tiff = processor_tiff(rasters[sw]; zeroed = sw != 2, pad = sw == 3 ? 1 : 0)
                push!(members, "$name.SAFE/measurement/$stem.tiff" => tiff)
            end
            zp = deflated_zip(joinpath(dir, "$name.zip"), members)
            out = mkpath(joinpath(dir, "bursts"))

            bs = zip_bursts(zp; polarization = pol, dir = out)
            @test [(b.swath, b.burst) for b in bs] == [(sw, k) for sw in 1:3 for k in 1:n]
            @test all(b -> b.slc == name, bs)
            picked = [ZipBurst(zp, 2, pol, 3; dir = out), ZipBurst(zp, 1, pol, 2; dir = out),
                      ZipBurst(zp, 2, pol, 2; dir = out), ZipBurst(zp, 3, pol, 1; dir = out)]
            files = burst_files(picked)
            for (b, f) in zip(picked, files)
                l = L[b.swath]
                @test open_measurement(f)[:, :] == rasters[b.swath][((b.burst - 1) * l + 1):(b.burst * l), :]
            end
            # Read again, the rasters are not inflated a second time.
            t = mtime(first(files))
            @test zip_rasters(picked[1:1]) == files[1:1] && mtime(first(files)) == t

            q = Sentinel1Product(picked; orbit = p.eof)
            @test q.swaths == [1, 2, 3]
            @test [[b.burst for b in g] for g in q.bursts] == [[2], [2, 3], [1]]
            whole = Sentinel1Product(zp; orbit = p.eof, polarization = pol)
            a, w = SLCDatasets.annotation(q, 2), SLCDatasets.annotation(whole, 2)
            per_burst = (:burst_start, :first_valid_line, :last_valid_line, :first_valid_sample, :last_valid_sample)
            for f in fieldnames(SubswathAnnotation)
                f in (:azimuth_fm_rate, :doppler_centroid) && continue
                @test isequal(getfield(a, f), f in per_burst ? getfield(w, f)[2:3] : getfield(w, f))
            end
            r, = SLCDatasets._subswath_raster(q, 2)
            @test r[:, :] == rasters[2][(L[2] + 1):(3L[2]), :]
            @test q.path == asf_safe_name(picked)

            # The fixture's annotations carry neither burst times from the ascending node nor a grid.
            @test_throws "no `azimuthAnxTime` element" burst_footprints(zp, 2, pol)

            spec = zip_burst_spec(picked[1])
            @test spec == zp * "#IW2/$(uppercase(pol))/2" && is_zip_burst_spec(spec)
            @test ZipBurst(spec; dir = out) == picked[1]
            @test !is_zip_burst_spec(zp) && !is_zip_burst_spec(zp * "#IW4/HH/0")
            @test_throws "does not name a burst of an SLC zip" ZipBurst(zp * "#IW1/HH")
            @test_throws "not a Sentinel-1 IW SLC zip" ZipBurst(joinpath(dir, "x.zip"), 1, pol, 1)
            @test_throws "has $n bursts, so burst 5" burst_files([ZipBurst(zp, 1, pol, 5; dir = out)])
            @test_throws "not consecutive" Sentinel1Product([ZipBurst(zp, 1, pol, 1), ZipBurst(zp, 1, pol, 3)];
                                                            orbit = p.eof)
            @test_throws "named twice" Sentinel1Product([ZipBurst(zp, 1, pol, 1), ZipBurst(zp, 1, pol, 1)];
                                                        orbit = p.eof)
        end
    end
end

@testset "a measurement raster's samples are located from its head" begin
    px = Complex{Int16}[complex(Int16(i), Int16(j)) for i in 1:6, j in 1:5]
    for (zeroed, pad) in ((true, 0), (true, 1), (false, 0), (false, 1))
        t = processor_tiff(px; zeroed, pad)
        D = SLCDatasets._measurement_data_offset(t, 3, 5, 2, "fixture")
        @test reinterpret(Complex{Int16}, t[(D + 1):(D + 4 * 30)]) == vec(permutedims(px))
    end
    @test_throws "not laid out as a Sentinel-1 measurement raster" SLCDatasets._measurement_data_offset(
        processor_tiff(px), 2, 5, 2, "fixture")
    @test_throws "not a little-endian classic TIFF" SLCDatasets._measurement_data_offset(zeros(UInt8, 64), 3, 5, 2, "x")
end

const ZIP_BURST_ZIP = get(ENV, "SLCDATASETS_ZIP_BURST_ZIP", "")
const ZIP_BURST_URLS = get(ENV, "SLCDATASETS_ZIP_BURST_URLS", "")
const ZIP_BURST_ORBIT = get(ENV, "SLCDATASETS_ZIP_BURST_ORBIT", "")

if isempty(ZIP_BURST_ZIP) || isempty(ZIP_BURST_URLS) || isempty(ZIP_BURST_ORBIT)
    @info """no SLC zip named, so the comparison of its bursts with the extractor's is skipped. Set
             SLCDATASETS_ZIP_BURST_ZIP, SLCDATASETS_ZIP_BURST_URLS and SLCDATASETS_ZIP_BURST_ORBIT to run it."""
else
    @testset "a zip's bursts equal the extractor's" begin
        mktempdir() do dir
            extracted = [AsfBurst(strip(u); dir = mkpath(joinpath(dir, "asf"))) for u in split(ZIP_BURST_URLS, ',')]
            zipped = [ZipBurst(ZIP_BURST_ZIP, b.swath, b.polarization, b.burst; dir = mkpath(joinpath(dir, "zip")))
                      for b in extracted]
            @test all(((a, z),) -> a.slc == z.slc, zip(extracted, zipped))
            for (fa, fz) in zip(burst_files(extracted), burst_files(zipped))
                @test open_measurement(fa)[:, :] == open_measurement(fz)[:, :]
            end
            pa = Sentinel1Product(extracted; orbit = ZIP_BURST_ORBIT)
            pz = Sentinel1Product(zipped; orbit = ZIP_BURST_ORBIT)
            @test pa.path == pz.path
            for (x, y) in zip(pa.annotations, pz.annotations), f in fieldnames(SubswathAnnotation)
                @test isequal(getfield(x, f), getfield(y, f))
            end
        end
    end
end
