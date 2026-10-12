# Bursts as ASF's burst extractor delivers them.
#
# The URL rules and the product assembled from bursts are checked offline: the metadata cache is seeded
# with extractor-style metadata built from the committed fixture's annotations, split between two parent
# products the way consecutive slices of a datatake split them, and the stacked raster is read from TIFFs
# written here. Comparing an assembled product with the SAFE `burst2safe` builds from the same bursts
# needs both, so it runs only where they are named: set `SLCDATASETS_ASF_URLS` to the bursts' extractor
# URLs (comma-separated), `SLCDATASETS_ASF_SAFE` to the SAFE and `SLCDATASETS_ASF_ORBIT` to its orbit file.

using SLCDatasets
using SLCDatasets: asf_burst_url, ASF_BURST_HOST, BurstStack, MeasurementRaster, RangePolynomial,
                   SubswathAnnotation, open_measurement, parse_utc, annotation_xml
import DiskArrays
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

@testset "a burst is read from its extractor URL" begin
    slc = "S1C_IW_SLC__1SDV_20250416T010149_20250416T010216_001911_003B60_6252"
    url = "https://$ASF_BURST_HOST/$slc/IW1/VV/3.tiff"
    b = AsfBurst(url; dir = "/scratch")
    @test (b.slc, b.swath, b.polarization, b.burst, b.dir) == (slc, 1, "VV", 4, "/scratch")
    @test asf_burst_url(b, "tiff") == url
    @test is_asf_burst_url(url) && !is_asf_burst_url("https://example.com/$slc/IW1/VV/3.tiff")
    @test AsfBurst(replace(url, ".tiff" => ".xml")) == AsfBurst(slc, 1, "VV", 4)
    @test_throws "not a burst extractor URL" AsfBurst("https://example.com/$slc/IW1/VV/3.tiff")
    @test_throws "not a burst extractor URL" AsfBurst("https://$ASF_BURST_HOST/$slc/IW4/VV/3.tiff")
    @test_throws "not a burst extractor URL" AsfBurst("https://$ASF_BURST_HOST/$slc/IW1/VV/3.jpg")
end

# The extractor's metadata file for one product: every subswath's annotation, each the children of a
# `.SAFE` annotation's root under `<content>`.
function extractor_metadata(annotations::Vector{String}, pol::String)
    products = map(enumerate(annotations)) do (swath, doc)
        body = replace(doc, r"^\s*<\?xml[^>]*>\s*<product>"s => "", r"</product>\s*$"s => "")
        "<product><swath>IW$swath</swath><polarisation>$pol</polarisation><content>$body</content></product>"
    end
    return "<burst><metadata>" * join(products) * "</metadata></burst>"
end

# An annotation document with only bursts `keep` of its burst list.
function keep_bursts(doc::String, keep)
    m = match(r"<burstList count=\"\d+\">(.*)</burstList>"s, doc)
    items = [x.match for x in eachmatch(r"<burst>.*?</burst>"s, m[1])]
    return replace(doc, m.match => "<burstList count=\"$(length(keep))\">" * join(items[keep]) * "</burstList>")
end

if !isempty(S1_PRODUCTS)
    @testset "bursts assemble into the product burst2safe builds" begin
        p = first(S1_PRODUCTS)
        pol = uppercase(String(p.gold.polarization))
        safe = Sentinel1Product(p.safe; orbit = p.eof, polarization = lowercase(pol))
        docs = [annotation_xml(p.safe, s, lowercase(pol)) for s in 1:3]
        n = nbursts(safe, 2)
        @test n >= 5
        # Two slices of one datatake: the first lists bursts 1-3 of every subswath, the second the rest.
        mission = safe.annotations[1].mission
        first_slc = "$(mission)_IW_SLC__1SDV_20000101T000000_20000101T000030_000001_000001_AAAA"
        second_slc = "$(mission)_IW_SLC__1SDV_20000101T000025_20000101T000055_000001_000001_BBBB"
        SLCDatasets._ASF_METADATA[(first_slc, pol)] = extractor_metadata([keep_bursts(d, 1:3) for d in docs], pol)
        SLCDatasets._ASF_METADATA[(second_slc, pol)] = extractor_metadata([keep_bursts(d, 4:n) for d in docs], pol)
        picked = [AsfBurst(first_slc, 2, pol, 2), AsfBurst(second_slc, 2, pol, 2), AsfBurst(first_slc, 2, pol, 3),
                  AsfBurst(second_slc, 2, pol, 1), AsfBurst(first_slc, 1, pol, 3)]

        q = Sentinel1Product(picked; orbit = p.eof)
        @test q.swaths == [1, 2]
        @test [length(b) for b in q.bursts] == [1, 4]
        # Ordered in azimuth, whatever order they were named in.
        @test [(b.slc, b.burst) for b in q.bursts[2]] ==
              [(first_slc, 2), (first_slc, 3), (second_slc, 1), (second_slc, 2)]
        a, whole = SLCDatasets.annotation(q, 2), SLCDatasets.annotation(safe, 2)
        per_burst = (:burst_start, :first_valid_line, :last_valid_line, :first_valid_sample, :last_valid_sample)
        # `isequal`: the fixture carries no steering rate, which reads as `NaN`.
        for f in fieldnames(SubswathAnnotation)
            @test isequal(getfield(a, f), f in per_burst ? getfield(whole, f)[2:5] : getfield(whole, f))
        end
        @test only(SLCDatasets.annotation(q, 1).burst_start) == SLCDatasets.annotation(safe, 1).burst_start[3]

        t = SLCDatasets._name_time
        @test q.path == asf_safe_name(picked)
        @test q.path == "$(mission)_IW_SLC__1S" * (pol == "HH" ? "SH" : "SV") * "_" *
                        t(min(whole.burst_start[2], SLCDatasets.annotation(safe, 1).burst_start[3])) * "_" *
                        t(whole.burst_start[5]) * "_000001_000001_0000"
        @test Sentinel1Product(picked; orbit = p.eof, swaths = [2]).swaths == [2]

        @test_throws "no burst of" Sentinel1Product(picked; orbit = p.eof, swaths = [3])
        @test_throws "not consecutive" Sentinel1Product(
            [AsfBurst(first_slc, 2, pol, 1), AsfBurst(first_slc, 2, pol, 3)]; orbit = p.eof)
        @test_throws "named twice" Sentinel1Product(
            [AsfBurst(first_slc, 2, pol, 1), AsfBurst(first_slc, 2, pol, 1)]; orbit = p.eof)
        @test_throws "of its product has 3 bursts" asf_safe_name([AsfBurst(first_slc, 2, pol, 4)])
        @test_throws "not one datatake" asf_safe_name(
            [AsfBurst(first_slc, 2, pol, 1),
             AsfBurst(replace(second_slc, "000001_000001" => "000002_000001"), 2, pol, 1)])
        other = pol == "HH" ? "HV" : "VH"
        @test_throws "share a polarization" asf_safe_name(
            [AsfBurst(first_slc, 2, pol, 1), AsfBurst(first_slc, 2, other, 2)])
    end
end

@testset "estimate lists merge as burst2safe merges them" begin
    at(s) = parse_utc("2024-06-18T02:55:" * s)
    poly(s) = RangePolynomial(at(s), 800_000.0, (1.0, 0.0, 0.0))
    first_list = [poly("10.0"), poly("12.0"), poly("14.0")]
    second_list = [poly("13.0"), poly("14.0"), poly("16.0")]
    merged = SLCDatasets._merged([second_list, first_list])
    @test [p.time for p in merged] == at.(["10.0", "12.0", "14.0", "16.0"])
    @test SLCDatasets._merged([RangePolynomial[], first_list]) == first_list
    @test isempty(SLCDatasets._merged([RangePolynomial[], RangePolynomial[]]))
    # Strictly inside the buffer on either side.
    @test [p.time for p in SLCDatasets._within(merged, at("13.0"), at("14.0"), 1.0)] == at.(["14.0"])
    @test [p.time for p in SLCDatasets._within(merged, at("13.0"), at("14.5"), 1.5)] == at.(["12.0", "14.0"])
end

@testset "a subswath's bursts read as one stacked raster" begin
    mktempdir() do dir
        a = Complex{Int16}[complex(Int16(10i + j), Int16(-i)) for i in 1:4, j in 1:6]
        b = Complex{Int16}[complex(Int16(100i + j), Int16(i)) for i in 1:4, j in 1:5]
        ra = open_measurement(write_tiff(joinpath(dir, "a.tiff"), a))
        rb = open_measurement(write_tiff(joinpath(dir, "b.tiff"), b))
        @test ra isa MeasurementRaster
        m = BurstStack([ra, rb], 4, 6)
        @test size(m) == (8, 6)
        # The narrower burst is zero past its own width.
        @test m[:, :] == [a; [b zeros(Complex{Int16}, 4)]]
        @test m[3:6, 5:6] == [a[3:4, 5:6]; [b[1:2, 5] zeros(Complex{Int16}, 2)]]
        @test m[5, 2] == b[1, 2]
        @test length(DiskArrays.eachchunk(m)) == 2
    end
end

const ASF_URLS = get(ENV, "SLCDATASETS_ASF_URLS", "")
const ASF_SAFE = get(ENV, "SLCDATASETS_ASF_SAFE", "")
const ASF_ORBIT = get(ENV, "SLCDATASETS_ASF_ORBIT", "")

if isempty(ASF_URLS) || isempty(ASF_SAFE) || isempty(ASF_ORBIT)
    @info """no extracted bursts named, so the comparison with burst2safe's SAFE is skipped. Set
             SLCDATASETS_ASF_URLS, SLCDATASETS_ASF_SAFE and SLCDATASETS_ASF_ORBIT to run it."""
else
    @testset "extracted bursts equal the SAFE burst2safe builds from them" begin
        mktempdir() do dir
            bs = [AsfBurst(strip(u); dir) for u in split(ASF_URLS, ',')]
            p = Sentinel1Product(bs; orbit = ASF_ORBIT)
            q = Sentinel1Product(ASF_SAFE; orbit = ASF_ORBIT, polarization = p.polarization, swaths = p.swaths)
            for (x, y) in zip(p.annotations, q.annotations), f in fieldnames(SubswathAnnotation)
                @test isequal(getfield(x, f), getfield(y, f))
            end
            # The SAFE's name ends with a checksum of its manifest where this one ends with `0000`.
            @test chop(p.path; tail = 4) == chop(basename(rstrip(ASF_SAFE, '/')); tail = 9)
            for sw in p.swaths
                rp, = SLCDatasets._subswath_raster(p, sw)
                rq, = SLCDatasets._subswath_raster(q, sw)
                @test size(rp) == size(rq)
                L = SLCDatasets.annotation(p, sw).lines_per_burst
                for r0 in (1, L - 3, size(rp, 1) - 50)
                    @test rp[r0:(r0 + 50), :] == rq[r0:(r0 + 50), :]
                end
            end
        end
    end
end
