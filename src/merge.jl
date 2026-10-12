# Folding a subswath's bursts into one acquisition.
#
# The bursts of a subswath are already `SLC`s, so the merge takes them rather than a product and a range:
# whatever vector a caller has — a whole subswath, a slice of one, the result of a `filter` — is the
# selection. What comes back is an `SLC` like any other, so everything downstream of this package reads
# it without knowing it was merged.
#
# The work that is not layout is reconciling the clocks. Each burst reports its times against an epoch
# derived from its own start, so a subswath of nine bursts carries nine different epochs and their
# `sensing_start` values are not on one scale — differencing them across a subswath is out by the whole
# spread of the epochs, some tens of seconds. The merged acquisition adopts the first burst's epoch and
# re-references every other burst onto it, which is the same reconciliation `repeat_interval` performs
# for a pair of products.

"""
    MergedBurstBackend <: AbstractSLCBackend

Several bursts of one Sentinel-1 subswath, described as a single acquisition.

Built by [`merge_bursts`](@ref). Holds the parsed product, the subswath, and the layout the bursts were
placed with, so the geometry and the samples both derive from one grid.
"""
struct MergedBurstBackend{B<:AbstractSLCBackend} <: AbstractSLCBackend
    # The backends of the bursts merged, in azimuth order and parallel to `grid.placements`. Holding
    # them rather than a product is what lets bursts of a `.SAFE` and bursts delivered one file each
    # merge by the same path: each knows how to reach its own samples.
    sources::Vector{B}
    annotation::SubswathAnnotation
    swath::Int
    bursts::UnitRange{Int}
    grid::BurstGrid
    orbit_path::String
end

# A merge of bursts is as TOPS as the bursts it merged. `MergedBurstBackend` subtypes
# `AbstractSLCBackend` rather than `AbstractBurstBackend` — it describes several bursts, not one — so it
# does not inherit the burst answer and states it here, deferring to its sources rather than asserting it:
# were a non-TOPS sensor ever merged this way, the answer would follow the data.
is_tops(b::MergedBurstBackend) = any(is_tops, b.sources)

_path(b::MergedBurstBackend) = _path(first(b.sources))
_leading_annotation(b::MergedBurstBackend) = b.annotation

# The epoch the merged times are reported against: the first burst placed, as for any single burst. So a
# merge of one burst reports exactly what that burst reports.
_anchor(b::MergedBurstBackend) = b.annotation.burst_start[first(b.bursts)]

# `burst` counts this acquisition's bursts from one, so it indexes `b.bursts` rather than the annotation's
# whole list: a merge of bursts 5:8 has bursts 1 through 4.
function deramp_parameters(b::MergedBurstBackend, burst::Integer)
    n = length(b.bursts)
    (1 <= burst <= n) || throw(ArgumentError(
        "this merge spans $n bursts, so burst $burst is not one of them"))
    return _deramp_parameters(b.annotation, b.bursts[burst])
end

"""
    burst_at(s::SLC, line::Integer) -> (burst, burst_line)

Which burst a line of a merged image belongs to, and which line of that burst it is.

Both count from one: `burst` indexes the acquisition's own bursts, as [`deramp_parameters`](@ref) takes it,
and `burst_line` is the row within that burst.

This is what a consumer removing the TOPS azimuth ramp needs. The ramp is quadratic about each burst's
*own* center, so a merged image carries one ramp per burst rather than a single one, and a phase evaluated
against merged-grid rows alone would be wrong by up to half a burst. Merging places bursts at whole-row
offsets and divides their overlap at the midpoint, so every row belongs to exactly one burst.

Throws for a line outside the image, and for an acquisition that is not a merge — a single burst is all its
own lines, and a mosaic has no single burst structure.
"""
burst_at(s::AbstractSLC, line::Integer) = burst_at(s.backend, line)

function burst_at(b::MergedBurstBackend, line::Integer)
    (1 <= line <= b.grid.nlines) || throw(ArgumentError(
        "this merged image has $(b.grid.nlines) lines, so line $line is not one of them"))
    for (k, p) in enumerate(b.grid.placements)
        line in p.grid_rows && return (k, p.burst_rows[line - first(p.grid_rows) + 1])
    end
    # Placements cover the valid region of each burst, and merging trims the grid to those; a row between
    # two of them would mean the grid outran the data it was built from.
    throw(ArgumentError(
        "line $line of this merged image falls in no burst's valid region, which a merged grid should " *
        "not contain. This is a bug in the merge rather than a bad argument."))
end

"""
    merge_bursts(bursts::AbstractVector{<:SLC}; tolerance = MAX_BURST_GRID_RESIDUAL) -> SLC

Place consecutive bursts of one subswath on a single uniform azimuth grid.

The bursts must come from one product, subswath and polarization, and be consecutive and in order —
which is what [`bursts`](@ref) returns, so a whole subswath or any slice of one merges directly. The
result is an [`SLC`](@ref) whose `nlines` spans from the first burst's first line to the last burst's
last, with the overlap between neighbours divided at its midpoint rather than repeated.

A line's azimuth time is `sensing_start + line / prf` throughout the result, which is what
[`RadarGeometry`](@ref) means by `prf` and what a stack of whole bursts would break at every seam. The
samples are placed rather than resampled: each burst's start lands on the shared grid within a fraction
of a line, and `tolerance` is how far it may be off before the merge is refused.

Use [`pixels`](@ref) or [`amplitude`](@ref) for the samples, which are read as they are indexed, and
[`validmask`](@ref) for which of them were imaged.

# Examples

```julia
b = bursts(safe; orbit = eof, swath = 2)
ref = merge_bursts(b)          # the whole subswath
ref = merge_bursts(b[3:7])     # five of its bursts
nlines(ref)
amplitude(ref)[1:512, 1:512]
```
"""
function merge_bursts(slcs::AbstractVector{<:AbstractSLC};
                      tolerance::Real = MAX_BURST_GRID_RESIDUAL)
    isempty(slcs) && throw(ArgumentError(
        "a merged acquisition needs at least one burst, but no bursts were given"))

    backends = [_merge_backend(s, i) for (i, s) in enumerate(slcs)]
    swath, range = _merge_extent(backends)
    a = _leading_annotation(first(backends))
    _check_merge_agreement(slcs, a, swath)

    grid = burst_grid(a, range; tolerance)
    return SLC(MergedBurstBackend(backends, a, swath, range, grid,
                                  _merge_orbit_path(backends)))
end

# The backend of one burst, refused unless it is one. A mosaic is a `Sentinel1Backend` without a subswath
# to place bursts in, so it is turned away by name; anything that is not a burst backend at all cannot
# answer the questions a merge asks of one.
_merge_backend(s::AbstractSLC, i::Integer) = _as_burst(s.backend, i)

function _as_burst(b::Sentinel1Backend, i::Integer)
    b.swath === nothing && throw(ArgumentError(
        "acquisition $i describes the mosaic across `$(_path(b))`'s subswaths rather than a burst of " *
        "one. Merge the bursts of a single subswath, which `bursts(path; swath)` returns"))
    return b
end

_as_burst(b::AbstractBurstBackend, ::Integer) = b

_as_burst(b::AbstractSLCBackend, i::Integer) = throw(ArgumentError(
    "acquisition $i is a $(nameof(typeof(b))) rather than a Sentinel-1 burst, and only Sentinel-1 " *
    "bursts merge"))

# The one subswath and consecutive burst range the given bursts amount to.
function _merge_extent(backends::AbstractVector)
    swath = burst_swath(first(backends))
    same = typeof(first(backends))
    source = burst_source(first(backends))
    polarization = burst_polarization(first(backends))
    for (i, b) in enumerate(backends)
        b isa same || throw(ArgumentError(
            "burst $i is a $(nameof(typeof(b))) but the first a $(nameof(same)); bursts reached " *
            "different ways are not known to come from one product"))
        # Consecutive burst numbers of two products are not consecutive bursts of anything. Two slices
        # of one datatake share an absolute orbit and a range origin, so nothing further downstream
        # would catch it: the grid would place one product's bursts by the other's annotation.
        burst_source(b) == source || throw(ArgumentError(
            "burst $i comes from `$(burst_source(b))` but the first from `$source`; bursts of " *
            "different products describe different acquisitions and do not merge"))
        burst_swath(b) == swath || throw(ArgumentError(
            "burst $i is in subswath IW$(burst_swath(b)) but the first in IW$swath. Subswaths " *
            "lie at different slant ranges, so merging them would need one range origin for two, " *
            "and their samples are not on one grid"))
        # Two channels of one subswath share their mission, orbit and range geometry, so this is the
        # only check that separates them; merging them would interleave two channels' samples.
        burst_polarization(b) == polarization || throw(ArgumentError(
            "burst $i is polarization $(uppercase(burst_polarization(b))) but the first " *
            "$(uppercase(polarization)); a merged subswath carries one channel"))
    end

    idx = [burst_index(b) for b in backends]
    expected = first(idx):(first(idx) + length(idx) - 1)
    idx == collect(expected) || throw(ArgumentError(
        "the bursts given are " * join(map(string, idx), ", ") * ", which is not a consecutive ascending " *
        "run. A merged subswath " *
        "spans one unbroken range of bursts, since a gap would leave rows no burst covers and a " *
        "repeat would place one twice"))
    return swath, expected
end

function _merge_orbit_path(backends::AbstractVector)
    path = orbit_path(first(backends))
    for (i, b) in enumerate(backends)
        orbit_path(b) == path || throw(ArgumentError(
            "burst $i takes its state vectors from `$(orbit_path(b))` but the first from " *
            "`$path`; one acquisition has one orbit"))
    end
    return path
end

# What the bursts report about themselves, against the annotation the grid is built from. `_merge_extent`
# has already established they name one product, subswath and channel; this catches an acquisition whose
# own geometry does not match that annotation, which a caller assembling backends by hand could produce.
function _check_merge_agreement(slcs, a::SubswathAnnotation, swath::Integer)
    ref = first(slcs).identification
    for (i, s) in enumerate(slcs)
        id = s.identification
        id.mission == ref.mission || throw(ArgumentError(
            "burst $i is from $(id.mission) but the first from $(ref.mission)"))
        id.absolute_orbit == ref.absolute_orbit || throw(ArgumentError(
            "burst $i is on absolute orbit $(id.absolute_orbit) but the first on " *
            "$(ref.absolute_orbit); bursts of different orbits are different acquisitions"))
        # Every burst of a subswath shares its range geometry, so a disagreement here means the
        # annotation they were built from is not the one being placed.
        s.geometry.starting_range == a.starting_range || throw(ArgumentError(
            "burst $i starts at range $(s.geometry.starting_range) m but subswath IW$swath at " *
            "$(a.starting_range) m"))
    end
    return nothing
end

# The merged window: the whole grid, from the anchor burst's first line, one line per `1/prf` throughout.
# That last part is what makes the window and the row count describe each other, and it is the property a
# stack of whole bursts would not have.
function _merged_window(b::MergedBurstBackend)
    start = _anchor(b)
    return start, _advance(start,
                           (b.grid.nlines - 1) * _leading_annotation(b).azimuth_time_interval)
end

read_geometry(b::MergedBurstBackend) =
    _s1_geometry(_leading_annotation(b), _anchor(b), b.grid.nlines, b.grid.nsamples)

read_identification(b::MergedBurstBackend) =
    _s1_identification(_leading_annotation(b), _merged_window(b)...)

# The orbit must cover the merged window rather than one burst's, since a solve anywhere in the image has
# to interpolate rather than extrapolate.
function read_orbit(b::MergedBurstBackend)
    start, stop = _merged_window(b)
    return _s1_orbit(b.orbit_path, start, stop,
                     "the merged window of `$(_path(b))` subswath IW$(b.swath)")
end

# Every burst backend reaches its own samples, so a merge asks each of its sources for a
# [`BurstRaster`](@ref) and places what comes back. A `.SAFE` answers with the subswath's raster and the
# burst's offset into it; an ASF burst with its own file at offset zero. Neither shape is special here.
#
# The sources are asked through `burst_rasters` rather than one at a time, because bursts that share a
# file share the work of opening it: asking each of a subswath's nine bursts on its own would map and
# validate the same raster nine times.
read_pixels(b::MergedBurstBackend) = ConcatenatedBursts(burst_rasters(b.sources), b.grid)

"""
    grid(s::SLC) -> BurstGrid

Where each burst sits in a merged acquisition.

Throws for an acquisition that is not a merge of bursts.
"""
grid(s::AbstractSLC) = grid(s.backend)
grid(b::MergedBurstBackend) = b.grid
grid(b::AbstractSLCBackend) = throw(ArgumentError(
    "a $(nameof(typeof(b))) is not a merge of bursts, so it has no burst grid"))
