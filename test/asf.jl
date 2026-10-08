# Bursts as ASF's burst extractor delivers them.
#
# The URL check needs no network. Everything else needs burst files, so it runs only where a directory
# holding them is named: set `SLCDATASETS_ASF_DIR` to one, together with the granule the files came from
# and an orbit file for it. `asf_bursts` fetches into such a directory and does not fetch again, so a
# directory filled by an earlier run is what these read.

using SLCDatasets
using SLCDatasets: AsfBurstBackend, asf_burst_url, asf_annotation, ASF_BURST_HOST,
                   read_annotation, annotation_xml, measurement_path, open_measurement,
                   Sentinel1Product, grid
using Test

@testset "a burst names its files" begin
    b = AsfBurst("S1A_IW_SLC__1SSH_20151120T080202", 2, "hh", 4)
    # The extractor numbers bursts from zero; a burst has the same number here it has everywhere else,
    # so burst 4 is index 3 in the URL.
    @test asf_burst_url(b, "tiff") ==
          "https://$ASF_BURST_HOST/S1A_IW_SLC__1SSH_20151120T080202/IW2/HH/3.tiff"
    @test asf_burst_url(b, "xml") ==
          "https://$ASF_BURST_HOST/S1A_IW_SLC__1SSH_20151120T080202/IW2/HH/3.xml"
    @test AsfBurst("g", 1, "vv", 1).polarization == "VV"

    @test_throws "subswaths are 1, 2 and 3" AsfBurst("g", 4, "hh", 1)
    @test_throws "subswaths are 1, 2 and 3" AsfBurst("g", 0, "hh", 1)
    @test_throws "numbered from 1" AsfBurst("g", 1, "hh", 0)
    @test_throws "HH, HV, VH or VV" AsfBurst("g", 1, "xx", 1)
    @test_throws "separate POEORB/RESORB" open_slc(AsfBurst("g", 1, "hh", 1))
    @test_throws "not a readable orbit file" open_slc(AsfBurst("g", 1, "hh", 1); orbit = "absent.EOF")
end

const ASF_DIR = get(ENV, "SLCDATASETS_ASF_DIR", "")
const ASF_GRANULE = get(ENV, "SLCDATASETS_ASF_GRANULE", "")
const ASF_ORBIT = get(ENV, "SLCDATASETS_ASF_ORBIT", "")
const ASF_SAFE = get(ENV, "SLCDATASETS_ASF_SAFE", "")

if isempty(ASF_DIR) || isempty(ASF_GRANULE) || isempty(ASF_ORBIT)
    @info """no ASF burst directory found, so those tests are skipped. Set SLCDATASETS_ASF_DIR,
             SLCDATASETS_ASF_GRANULE and SLCDATASETS_ASF_ORBIT to read burst files, and
             SLCDATASETS_ASF_SAFE to compare them against the same product's `.SAFE`."""
else
    const ASF_SWATH = parse(Int, get(ENV, "SLCDATASETS_ASF_SWATH", "2"))
    const ASF_POL = uppercase(get(ENV, "SLCDATASETS_ASF_POL", "HH"))
    # Which bursts the directory holds, as 1-based indices.
    const ASF_BURSTS = [parse(Int, s) for s in
                        split(get(ENV, "SLCDATASETS_ASF_BURSTS", ""), ','; keepempty = false)]

    @testset "a burst's metadata carries its whole subswath" begin
        b = AsfBurst(ASF_GRANULE, ASF_SWATH, ASF_POL, first(ASF_BURSTS); dir = ASF_DIR)
        a = read_annotation(asf_annotation(b), ASF_SWATH)
        # The file wraps the full annotation, so it lists every burst of the parent granule rather
        # than the one fetched.
        @test nbursts(a) > 1
        @test a.swath == ASF_SWATH
        @test lowercase(a.polarization) == lowercase(ASF_POL)
        @test a.lines_per_burst > 0
    end

    # A burst is the same acquisition however it was delivered, which is the claim that lets one merge
    # serve both. Checked against the `.SAFE` of the same granule where one is available.
    if !isempty(ASF_SAFE)
        @testset "an ASF burst equals the same burst of the .SAFE" begin
            raster = open_measurement(measurement_path(
                Sentinel1Product(ASF_SAFE; orbit = ASF_ORBIT, polarization = lowercase(ASF_POL)),
                ASF_SWATH))
            a = read_annotation(annotation_xml(ASF_SAFE, ASF_SWATH, lowercase(ASF_POL)), ASF_SWATH)
            L = a.lines_per_burst

            for i in ASF_BURSTS
                from_asf = open_slc(AsfBurst(ASF_GRANULE, ASF_SWATH, ASF_POL, i; dir = ASF_DIR);
                                    orbit = ASF_ORBIT)
                from_safe = open_slc(ASF_SAFE; orbit = ASF_ORBIT, swath = ASF_SWATH, burst = i,
                                     polarization = lowercase(ASF_POL))
                @test from_asf.geometry == from_safe.geometry
                @test from_asf.identification == from_safe.identification

                # The ASF file holds this burst alone, so its rows are the burst's own where the
                # `.SAFE` raster has them at that burst's offset into the stack.
                px = pixels(from_asf)
                @test size(px) == (L, a.samples_per_burst)
                for r in (1, a.first_valid_line[i], L), c in (1, a.first_valid_sample[i],
                                                              a.samples_per_burst)
                    @test px[r, c] == raster[(i - 1) * L + r, c]
                end
            end
        end

        # The point of the whole arrangement: two delivery formats, one merged array.
        if length(ASF_BURSTS) > 1 && ASF_BURSTS == first(ASF_BURSTS):last(ASF_BURSTS)
            @testset "a merge of ASF bursts equals a merge of the same .SAFE bursts" begin
                range = first(ASF_BURSTS):last(ASF_BURSTS)
                from_asf = merge_bursts(asf_bursts(ASF_GRANULE, ASF_SWATH, ASF_POL, range;
                                                   orbit = ASF_ORBIT, dir = ASF_DIR))
                from_safe = merge_bursts(bursts(ASF_SAFE; orbit = ASF_ORBIT, swath = ASF_SWATH,
                                                polarization = lowercase(ASF_POL))[range])

                @test from_asf.geometry == from_safe.geometry
                @test from_asf.identification == from_safe.identification
                @test validmask(from_asf) == validmask(from_safe)

                pa, ps = pixels(from_asf), pixels(from_safe)
                @test size(pa) == size(ps)
                # Every sample of a band spanning the first seam, and of the top margin.
                seam = first(grid(from_asf).placements[2].grid_rows)
                @test pa[(seam - 4):(seam + 4), :] == ps[(seam - 4):(seam + 4), :]
                @test pa[1:16, :] == ps[1:16, :]
                @test amplitude(from_asf)[1:16, 1:64] == amplitude(from_safe)[1:16, 1:64]
            end
        end
    end

    @testset "bursts reached different ways do not merge together" begin
        i = first(ASF_BURSTS)
        from_asf = open_slc(AsfBurst(ASF_GRANULE, ASF_SWATH, ASF_POL, i; dir = ASF_DIR);
                            orbit = ASF_ORBIT)
        if !isempty(ASF_SAFE)
            from_safe = open_slc(ASF_SAFE; orbit = ASF_ORBIT, swath = ASF_SWATH, burst = i + 1,
                                 polarization = lowercase(ASF_POL))
            @test_throws "reached different ways" merge_bursts([from_asf, from_safe])
        end
    end
end
