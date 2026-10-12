# Rebuilding a Sentinel-1 SAFE product and its orbit file from the committed inputs.
#
# `reference/sentinel1_inputs*.json` holds the annotation fields and state vectors the reader reads,
# taken verbatim from the granules the golden values were dumped from — the text as the annotation
# writes it, not a re-rounded float. Writing them back into the group and element layout a real product
# uses gives the reader a product to open, so the whole Sentinel-1 path is exercised against ISCE3's
# values without the 4.8 GB of granules those values came from.
#
# `reference/extract_sentinel1_inputs.jl` is what produced these files.

using JSON

# The golden-value files, each naming the granule its values were dumped from.
const S1_GOLDEN = ["sentinel1_metadata.json", "sentinel1_metadata_s1b.json"]

# The annotation elements `read_annotation` looks for, in the paths it looks for them. Fields a real
# annotation carries but this reader ignores are omitted rather than filled with stand-ins, so what the
# reader depends on is visible here.
# The per-line `firstValidSample`/`lastValidSample` arrays of one burst, written back out from the
# rectangle the inputs record: the burst's invalid lines carry -1, its valid lines the sample bounds.
function s1_valid_sample_arrays(region::AbstractDict, lines_per_burst::Integer)
    lo = Int(region["firstValidLine"])
    hi = Int(region["lastValidLine"])
    first_sample = Int(region["firstValidSample"])
    last_sample = Int(region["lastValidSample"])
    valid(i) = lo <= i <= hi
    firsts = join((valid(i) ? first_sample : -1 for i in 1:lines_per_burst), " ")
    lasts = join((valid(i) ? last_sample : -1 for i in 1:lines_per_burst), " ")
    return firsts, lasts
end

function s1_annotation_document(sw::AbstractDict)
    lines_per_burst = parse(Int, String(sw["linesPerBurst"]))
    regions = sw["burstValidRegions"]
    bursts = join(begin
                      firsts, lasts = s1_valid_sample_arrays(regions[i], lines_per_burst)
                      "<burst><azimuthTime>$t</azimuthTime>" *
                      "<firstValidSample>$firsts</firstValidSample>" *
                      "<lastValidSample>$lasts</lastValidSample></burst>"
                  end
                  for (i, t) in enumerate(sw["burstAzimuthTimes"]))
    steer = haskey(sw, "azimuthSteeringRate") ?
            "<azimuthSteeringRate>$(sw["azimuthSteeringRate"])</azimuthSteeringRate>" : ""
    fm_rate = s1_polynomial_list(get(sw, "azimuthFmRateList", nothing), "azimuthFmRateList",
                                 "azimuthFmRatePolynomial")
    # `dcEstimateList` hangs off `dopplerCentroid`, not `generalAnnotation`.
    dc_list = s1_polynomial_list(get(sw, "dcEstimateList", nothing), "dcEstimateList",
                                 "dataDcPolynomial")
    doppler = isempty(dc_list) ? "" : "<dopplerCentroid>$dc_list</dopplerCentroid>"
    return """<?xml version="1.0" encoding="UTF-8"?>
    <product>
      <adsHeader>
        <missionId>$(sw["missionId"])</missionId>
        <productType>$(sw["productType"])</productType>
        <polarisation>$(sw["polarisation"])</polarisation>
        <absoluteOrbitNumber>$(sw["absoluteOrbitNumber"])</absoluteOrbitNumber>
      </adsHeader>
      <generalAnnotation>
        <productInformation>
          <pass>$(sw["pass"])</pass>
          <rangeSamplingRate>$(sw["rangeSamplingRate"])</rangeSamplingRate>
          <radarFrequency>$(sw["radarFrequency"])</radarFrequency>
          $steer
        </productInformation>
        $fm_rate
      </generalAnnotation>
      $doppler
      <imageAnnotation>
        <imageInformation>
          <azimuthTimeInterval>$(sw["azimuthTimeInterval"])</azimuthTimeInterval>
          <slantRangeTime>$(sw["slantRangeTime"])</slantRangeTime>
        </imageInformation>
      </imageAnnotation>
      <swathTiming>
        <linesPerBurst>$(sw["linesPerBurst"])</linesPerBurst>
        <samplesPerBurst>$(sw["samplesPerBurst"])</samplesPerBurst>
        <burstList count="$(length(sw["burstAzimuthTimes"]))">$bursts</burstList>
      </swathTiming>
    </product>
    """
end

# A list of polynomial elements, or nothing at all where the input carries none — which is the state of the
# committed golden inputs, since the deramp fields were not among the values dumped from those granules.
# So the reader's optional-field path is what those exercise, and a test wanting the fields supplies them.
#
# An entry may set `"spelling" => "old"` to write the coefficients as separate sibling elements rather than
# as one text node, which is how products processed before IPF ~2.36-2.82 carry them. Both spellings have
# to be read, so both have to be writable here.
function s1_polynomial_list(entries, list_name::AbstractString, poly_name::AbstractString)
    entries === nothing && return ""
    items = map(entries) do e
        coeffs = e["coefficients"]
        body = if get(e, "spelling", "new") == "old"
            join("<c$i>$c</c$i>" for (i, c) in enumerate(coeffs))
        else
            "<$poly_name>$(join(coeffs, " "))</$poly_name>"
        end
        "<item><azimuthTime>$(e["azimuthTime"])</azimuthTime><t0>$(e["t0"])</t0>$body</item>"
    end
    return "<$list_name count=\"$(length(items))\">$(join(items))</$list_name>"
end

function s1_eof_document(vectors)
    osvs = join(("<OSV>" * join(("<$f>$(v[f])</$f>" for f in ("UTC", "X", "Y", "Z", "VX", "VY", "VZ"))) *
                 "</OSV>" for v in vectors))
    return """<?xml version="1.0"?>
    <Earth_Explorer_File><Data_Block type="xml"><List_of_OSVs count="$(length(vectors))">$osvs</List_of_OSVs></Data_Block></Earth_Explorer_File>
    """
end

"""
    write_s1_fixture(dir, inputs) -> (safe, eof)

Write `inputs` as a `.SAFE` directory and a POEORB `.EOF` under `dir`.

The orbit file keeps the name from the granule it was extracted from, because the reader reports
POEORB or RESORB by matching that name.
"""
function write_s1_fixture(dir::AbstractString, inputs)
    pol = String(inputs["polarization"])
    mission = lowercase(String(inputs["swaths"]["1"]["missionId"]))
    safe = joinpath(dir, "fixture.SAFE")
    mkpath(joinpath(safe, "annotation", "calibration"))
    write(joinpath(safe, "manifest.safe"), "<xfdu:XFDU/>")
    for swath in 1:3
        sw = inputs["swaths"][string(swath)]
        name = "$mission-iw$swath-slc-$pol-20000101t000000-20000101t000030-000001-000001-00$swath.xml"
        write(joinpath(safe, "annotation", name), s1_annotation_document(sw))
        # The parallel `calibration` tree holds same-named files the reader must not pick up.
        write(joinpath(safe, "annotation", "calibration", name), "<product/>")
    end
    eof = joinpath(dir, String(inputs["orbit_file"]))
    write(eof, s1_eof_document(inputs["state_vectors"]))
    return safe, eof
end

"""
    s1_fixtures() -> Vector{NamedTuple}

The committed Sentinel-1 cases: each golden-value file paired with the inputs to rebuild its product.
"""
function s1_fixtures()
    out = NamedTuple[]
    for name in S1_GOLDEN
        gold_file = joinpath(@__DIR__, "reference", name)
        input_file = joinpath(@__DIR__, "reference", replace(name, "_metadata" => "_inputs"))
        (isfile(gold_file) && isfile(input_file)) || continue
        push!(out, (; name,
                    gold = JSON.parse(read(gold_file, String)),
                    inputs = JSON.parse(read(input_file, String); dicttype = Dict{String,Any})))
    end
    return out
end
