# Bursts as ASF's burst extractor delivers them (`extractor.jl`), one file per burst: each read as an
# acquisition of its own, merged, or assembled into the product `burst2safe` builds from them.
#
# The raster has the same layout a `.SAFE` measurement file does — uncompressed, one strip per line,
# complex 16-bit samples — so `open_measurement` reads it unchanged. Unlike a `.SAFE`, each is its own file
# holding one burst, which is why a merge places them at row offset zero rather than at a burst's offset
# into a stacked raster.
#
# The metadata wraps the *whole* product's annotation for every subswath:
# `<burst><metadata><product><swath>IW2</swath><content>` holds a complete product annotation listing
# every burst of the parent SLC, not just the one fetched. So `read_annotation` parses the `<content>`
# subtree as it stands, and which burst the file holds is known from its index rather than from the XML.

"""
    asf_annotation(b::AsfBurst) -> Node

The product annotation of the subswath a burst belongs to, from the burst's metadata file.

The file wraps a full annotation per subswath and polarization of the parent granule, each under a
`<content>` element whose children are the sections a `.SAFE` annotation has at its root. So that
element *is* the annotation, and what is returned lists every burst of the subswath rather than only
this one.
"""
function asf_annotation(b::AsfBurst)
    # The document is held by the returned node, so it must stay reachable: `parsexml`'s result is kept
    # alive by the node's own reference to its document.
    r = root(parsexml(only(asf_metadata([b]))))
    want = string("IW", b.swath)
    for node in findall("metadata/product", r)
        _findtext(node, "swath") == want || continue
        uppercase(_findtext(node, "polarisation")) == b.polarization || continue
        content = findfirst("content", node)
        content === nothing && break
        return content
    end
    throw(ArgumentError(
        "the metadata of `$(b.slc)` carries no product annotation for subswath IW$(b.swath) \
         polarization $(b.polarization)"))
end

# Each product's subswath annotations as parsed, so a product's metadata is parsed once per subswath
# however many of its bursts are read.
const _ASF_ANNOTATIONS = Dict{Tuple{String,String,Int},SubswathAnnotation}()

"""
    asf_subswath(b::AsfBurst) -> SubswathAnnotation

The annotation of the subswath `b` lies in, as its parent product lists it — every burst of the
subswath, not only `b`.
"""
function asf_subswath(b::AsfBurst)
    key = (b.slc, b.polarization, b.swath)
    @lock _ASF_METADATA_LOCK begin
        a = get(_ASF_ANNOTATIONS, key, nothing)
        a === nothing || return a
    end
    a = read_annotation(asf_annotation(b), b.swath)
    @lock _ASF_METADATA_LOCK _ASF_ANNOTATIONS[key] = a
    return a
end

"""
    AsfBurstBackend <: AbstractSLCBackend

One ASF-delivered burst: its own raster, and the annotation of the subswath it belongs to.
"""
struct AsfBurstBackend <: AbstractBurstBackend
    source::AsfBurst
    annotation::SubswathAnnotation
    orbit::String
end

_path(b::AsfBurstBackend) = asf_burst_url(b.source, "tiff")
_leading_annotation(b::AsfBurstBackend) = b.annotation

burst_index(b::AsfBurstBackend) = b.source.burst
burst_swath(b::AsfBurstBackend) = b.source.swath
burst_polarization(b::AsfBurstBackend) = lowercase(b.source.polarization)
burst_source(b::AsfBurstBackend) = b.source.slc
orbit_path(b::AsfBurstBackend) = b.orbit

function _asf_check_burst(b::AsfBurstBackend)
    n = nbursts(b.annotation)
    i = b.source.burst
    1 <= i <= n || throw(ArgumentError(
        "subswath IW$(b.source.swath) of `$(b.source.slc)` has $n bursts, so burst $i does not exist"))
    return i
end

_anchor(b::AsfBurstBackend) = b.annotation.burst_start[_asf_check_burst(b)]

# A burst's geometry and identification are its own, exactly as for a burst read from a `.SAFE`: the
# annotation an ASF burst carries is the whole subswath's, so the same fields describe it and the same
# helpers build the records.
read_geometry(b::AsfBurstBackend) =
    _s1_geometry(b.annotation, _anchor(b), b.annotation.lines_per_burst,
                 b.annotation.samples_per_burst)

read_identification(b::AsfBurstBackend) =
    _s1_identification(b.annotation, _anchor(b), _burst_stop(b.annotation, _anchor(b)))

function read_orbit(b::AsfBurstBackend)
    start = _anchor(b)
    return _s1_orbit(b.orbit, start, _burst_stop(b.annotation, start),
                     "burst $(b.source.burst) of subswath IW$(b.source.swath)")
end

# The file holds this burst alone, so its rows are the burst's own.
function read_pixels(b::AsfBurstBackend)
    a = b.annotation
    file = localpath(b.source)
    raster = open_measurement(file)
    size(raster, 1) == a.lines_per_burst || throw(ArgumentError(
        "`$file` has $(size(raster, 1)) lines but subswath IW$(b.source.swath) records " *
        "$(a.lines_per_burst) per burst"))
    size(raster, 2) == a.samples_per_burst || throw(ArgumentError(
        "`$file` is $(size(raster, 2)) samples wide but subswath IW$(b.source.swath) " *
        "records $(a.samples_per_burst)"))
    return raster
end

"""
    open_slc(b::AsfBurst; orbit) -> SLC

One ASF-delivered burst as an acquisition.

`orbit` is a POEORB or RESORB `.EOF` file, required for the same reason it is for a `.SAFE`: a
Sentinel-1 product carries no state vectors.
"""
function open_slc(b::AsfBurst; orbit = nothing)
    orbit_path = _require_orbit(orbit, string("burst ", b.burst, " of `", b.slc, "`"))
    return SLC(AsfBurstBackend(b, asf_subswath(b), orbit_path))
end

"""
    asf_bursts(slc, swath, polarization, bursts; orbit, dir = nothing) -> Vector{SLC}

Several bursts of one subswath from ASF's burst extractor.

`bursts` is a range or vector of 1-based burst indices. The metadata file is fetched and parsed once
rather than per burst, since each of them carries the whole subswath's annotation; the rasters are
fetched as they are read.

The result merges with [`merge_bursts`](@ref).

# Examples

```julia
b = asf_bursts("S1A_IW_SLC__1SSH_20151120T080202_...", 2, "HH", 3:7; orbit = eof)
ref = merge_bursts(b)
```
"""
function asf_bursts(slc::AbstractString, swath::Integer, polarization::AbstractString,
                    bursts; orbit = nothing, dir = nothing)
    idx = collect(Int, bursts)
    isempty(idx) && throw(ArgumentError("no bursts were named"))
    orbit_path = _require_orbit(orbit, string("bursts of `", slc, "`"))

    # One fetch and one parse: every burst's file holds the same subswath annotation.
    a = asf_subswath(AsfBurst(slc, swath, polarization, first(idx); dir))
    return [SLC(AsfBurstBackend(AsfBurst(slc, swath, polarization, i; dir), a, orbit_path))
            for i in idx]
end

# How a merge of these reads its samples is in `merge.jl`, where `MergedBurstBackend` is defined.

# ---------------------------------------------------------------------------
# The product `burst2safe` assembles
# ---------------------------------------------------------------------------
#
# `burst2safe` builds one SAFE from a set of extracted bursts: per subswath, the bursts stacked in azimuth
# order in one measurement raster, under an annotation it builds from the bursts' parent products'
# (`burst2safe/product.py`). What that annotation holds of what this package reads:
#
# - the scalar fields — slant-range time, spacings, wavelength, steering rate, the header — are those of
#   the parent of the subswath's earliest burst;
# - `linesPerBurst` is the bursts' common line count, and `samplesPerBurst` the widest burst's width; a
#   narrower burst's samples start at the raster's first column, zero past its own width;
# - the burst list is the bursts given, each as its parent lists it;
# - the FM-rate and Doppler-centroid lists are the parents' lists merged — the first parent's whole, each
#   later one's from past the last time already kept — and cut to the estimates within 500 s (FM rate) or
#   3 s (Doppler) of the bursts' span, from the earliest burst's start to the latest burst's end.

"""
    Sentinel1Product(bursts::AbstractVector{AsfBurst}; orbit, swaths = nothing)

The product `burst2safe` assembles from extracted `bursts`, opened for `swaths` — by default every
subswath a burst lies in. Each subswath's bursts must be consecutive, of one polarization, with no burst
named twice. The metadata of the bursts' parent products is fetched here; every raster is downloaded,
into each burst's `dir`, when the first subswath's raster is read.

The product's `path` is the name `burst2safe` gives the SAFE (see [`asf_safe_name`](@ref)), since there
is no file: it is what the product is reported as.
"""
function Sentinel1Product(bursts::AbstractVector{AsfBurst}; orbit::AbstractString, swaths = nothing)
    name = asf_safe_name(bursts)
    have = sort!(unique([b.swath for b in bursts]))
    selected = swaths === nothing ? have : sort!(collect(Int, swaths))
    for s in selected
        s in have || throw(ArgumentError("no burst of `$name` lies in subswath IW$s"))
    end
    groups = Vector{AsfBurst}[]
    annotations = SubswathAnnotation[]
    for s in selected
        g, a = _assembled_swath(AsfBurst[b for b in bursts if b.swath == s], s)
        push!(groups, g)
        push!(annotations, a)
    end
    return Sentinel1Product(name, String(orbit), lowercase(first(bursts).polarization), selected, annotations,
                            groups, nothing)
end

"""
    asf_safe_name(bursts) -> String

The name `burst2safe` gives the SAFE it assembles from `bursts`, without `.SAFE`:
`<platform>_IW_SLC__1S<pol>_<first>_<last>_<orbit>_<datatake>_0000`, where `<first>` and `<last>` are
the earliest and latest burst start, to the second, and the platform, absolute orbit and datatake are
the parent products'. `burst2safe` ends the name with a checksum of the SAFE's manifest, which is
written with the time it was made; with no manifest, this ends it with the `0000` `burst2safe` names a
SAFE with before it has one.
"""
function asf_safe_name(bursts::AbstractVector{AsfBurst})
    isempty(bursts) && throw(ArgumentError("no bursts were named"))
    pol = first(bursts).polarization
    f = split(first(bursts).slc, '_')
    for b in bursts
        b.polarization == pol || throw(ArgumentError(
            "the bursts of one product share a polarization, but $(b.polarization) and $pol were named"))
        g = split(b.slc, '_')
        length(g) == 10 || throw(ArgumentError("`$(b.slc)` is not a Sentinel-1 product name"))
        (g[1] == f[1] && g[8] == f[8] && g[9] == f[9]) || throw(ArgumentError(
            "`$(b.slc)` and `$(first(bursts).slc)` are not one datatake, so their bursts make no product"))
    end
    asf_metadata(bursts)
    lo = hi = nothing
    for b in bursts
        a = asf_subswath(b)
        1 <= b.burst <= nbursts(a) || throw(ArgumentError(
            "`$(asf_burst_url(b, "tiff"))` names no burst: subswath IW$(b.swath) of its product has \
             $(nbursts(a)) bursts"))
        t = a.burst_start[b.burst]
        (lo === nothing || t < lo) && (lo = t)
        (hi === nothing || hi < t) && (hi = t)
    end
    # `burst2safe`'s code for a single-polarization product.
    code = pol == "HH" ? "SH" : pol == "VV" ? "SV" : pol
    # Built a few pieces at a time: a concatenation of many arguments is compiled for an argument tuple
    # of unknown length, which a trimmed executable cannot dispatch.
    head = f[1] * "_" * f[2] * "_" * f[3] * "__1S" * code
    times = _name_time(something(lo)) * "_" * _name_time(something(hi))
    return head * "_" * times * "_" * f[8] * "_" * f[9] * "_0000"
end

# `t` to the whole second as `yyyymmddTHHMMSS`, the form product and burst names carry.
function _name_time(t::UtcTime)
    d = t.datetime + Dates.Second(floor(Int, t.seconds))
    return _zpad(Dates.year(d), 4) * _zpad(Dates.month(d), 2) * _zpad(Dates.day(d), 2) * "T" *
           _zpad(Dates.hour(d), 2) * _zpad(Dates.minute(d), 2) * _zpad(Dates.second(d), 2)
end

# One subswath's bursts in azimuth order, and the annotation `burst2safe` writes for them.
function _assembled_swath(bs::Vector{AsfBurst}, swath::Int)
    parents = Dict{String,SubswathAnnotation}()
    for b in bs
        haskey(parents, b.slc) || (parents[b.slc] = asf_subswath(b))
    end
    start(b::AsfBurst) = parents[b.slc].burst_start[b.burst]
    sorted = sort(bs; by = start)
    for k in 2:length(sorted)
        prev, b = sorted[k - 1], sorted[k]
        gap = seconds_between(start(prev), start(b))
        gap > 0 || throw(ArgumentError(
            "`$(asf_burst_url(b, "tiff"))` is named twice, or as the same burst of another product"))
        # Consecutive bursts start a burst cycle apart, a little less than a burst's length.
        a = parents[prev.slc]
        gap < 1.5 * a.lines_per_burst * a.azimuth_time_interval || throw(ArgumentError(
            "the IW$swath bursts are not consecutive: `$(asf_burst_url(b, "tiff"))` starts \
             $(round(gap; digits = 2)) s after the burst before it"))
    end
    order = unique([b.slc for b in sorted])
    head = parents[first(order)]
    lpb = head.lines_per_burst
    for b in sorted
        parents[b.slc].lines_per_burst == lpb || throw(ArgumentError(
            "the IW$swath bursts differ in length: `$(asf_burst_url(b, "tiff"))` has \
             $(parents[b.slc].lines_per_burst) lines, the first $lpb"))
    end
    pick(f) = Int[f(parents[b.slc])[b.burst] for b in sorted]
    starts = UtcTime[start(b) for b in sorted]
    span_end = maximum(k -> _advance(starts[k], (lpb - 1) * parents[sorted[k].slc].azimuth_time_interval),
                       eachindex(sorted))
    lists(f) = [f(parents[p]) for p in order]
    annotation = SubswathAnnotation(
        swath, head.polarization, head.absolute_orbit, head.mission, head.product_type, head.pass_direction,
        head.starting_range, head.range_pixel_spacing, head.wavelength, head.azimuth_time_interval,
        lpb, maximum(b -> parents[b.slc].samples_per_burst, sorted),
        starts, pick(a -> a.first_valid_line), pick(a -> a.last_valid_line), pick(a -> a.first_valid_sample),
        pick(a -> a.last_valid_sample), head.azimuth_steering_rate,
        _within(_merged(lists(a -> a.azimuth_fm_rate)), first(starts), span_end, 500.0),
        _within(_merged(lists(a -> a.doppler_centroid)), first(starts), span_end, 3.0))
    return sorted, annotation
end

# Several products' estimate lists as one, in time order: the earliest list whole, then from each later
# one the estimates past the latest time already kept.
function _merged(lists::Vector{Vector{RangePolynomial}})
    nonempty = sort!(filter(!isempty, lists); by = l -> first(l).time)
    isempty(nonempty) && return RangePolynomial[]
    out = copy(first(nonempty))
    latest = maximum(p -> p.time, out)
    for l in nonempty[2:end]
        keep = filter(p -> latest < p.time, l)
        isempty(keep) && continue
        append!(out, keep)
        latest = maximum(p -> p.time, keep)
    end
    return out
end

# The estimates strictly within `buffer` seconds of `[lo, hi]`.
_within(list::Vector{RangePolynomial}, lo::UtcTime, hi::UtcTime, buffer::Float64) =
    filter(p -> seconds_between(lo, p.time) > -buffer && seconds_between(p.time, hi) > -buffer, list)

"""
    BurstStack <: DiskArrays.AbstractDiskArray{Complex{Int16},2}

A subswath's extracted bursts stacked as a measurement raster stacks them, indexed `(line, sample)`:
burst `k` in lines `(k - 1) L + 1` to `k L`, its samples from the first column and zero past its own
width. Each burst is read from its downloaded file a window at a time; a chunk is one burst.
"""
struct BurstStack <: DiskArrays.AbstractDiskArray{Complex{Int16},2}
    rasters::Vector{MeasurementRaster}
    lines_per_burst::Int
    samples::Int
end

Base.size(m::BurstStack) = (length(m.rasters) * m.lines_per_burst, m.samples)

DiskArrays.haschunks(::BurstStack) = DiskArrays.Chunked()
DiskArrays.eachchunk(m::BurstStack) = DiskArrays.GridChunks(size(m), (m.lines_per_burst, m.samples))

function DiskArrays.readblock!(m::BurstStack, out::AbstractMatrix, rows::AbstractUnitRange, cols::AbstractUnitRange)
    fill!(out, zero(Complex{Int16}))
    L = m.lines_per_burst
    for k in cld(first(rows), L):cld(last(rows), L)
        r = m.rasters[k]
        br = intersect(rows, ((k - 1) * L + 1):(k * L))
        bc = intersect(cols, 1:size(r, 2))
        isempty(bc) && continue
        out[br .- (first(rows) - 1), bc .- (first(cols) - 1)] = r[br .- (k - 1) * L, bc]
    end
    return out
end

# An assembled product's measurement raster for one subswath: its bursts, downloaded, stacked.
function _stacked_measurement(p::Sentinel1Product, swath::Integer)
    a = annotation(p, swath)
    bs = p.bursts[findfirst(==(Int(swath)), p.swaths)]
    rasters = MeasurementRaster[]
    # The first read downloads every subswath's rasters, as a zip's first read inflates every subswath:
    # requested together, they are generated together and download back to back.
    every = reduce(vcat, p.bursts)
    files = asf_rasters(every)
    for b in bs
        file = files[something(findfirst(==(b), every))]
        r = open_measurement(file)
        (size(r, 1) == a.lines_per_burst && size(r, 2) <= a.samples_per_burst) || throw(ArgumentError(
            "`$file` is not the shape of a subswath IW$swath burst, $(a.lines_per_burst) lines of at \
             most $(a.samples_per_burst) samples"))
        push!(rasters, r)
    end
    return (p.path * " IW" * string(swath), BurstStack(rasters, a.lines_per_burst, a.samples_per_burst))
end
