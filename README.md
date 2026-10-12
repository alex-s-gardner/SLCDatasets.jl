# SLCDatasets

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alex-s-gardner.github.io/SLCDatasets.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alex-s-gardner.github.io/SLCDatasets.jl/dev/)
[![Build Status](https://github.com/alex-s-gardner/SLCDatasets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/SLCDatasets.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/alex-s-gardner/SLCDatasets.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/alex-s-gardner/SLCDatasets.jl)

Read single-look complex SAR products into one type, whatever the sensor.

```julia
using SLCDatasets

s = open_slc("NISAR_L1_PR_RSLC_....h5")
s.geometry.starting_range   # 895255.2277025
s.geometry.prf              # 1520.0
s.geometry.look_side        # LookLeft
orbit(s)                    # state vectors, read on first access
```

`open_slc` returns a `SLC` carrying an `Identification` — what the acquisition is — and a
`RadarGeometry`: the slant-range/azimuth geometry in the units a geometry kernel wants. Both are parsed
when the product is opened, since a product whose metadata cannot be read is not usable. The state
vectors are deferred, and the sensor is chosen by inspecting the product rather than by the caller
naming it.

The field set is sensor-neutral by construction. `look_side` and `prf` are fields rather than constants
because the two reference loaders they follow disagree on both: NISAR looks left with the pulse
repetition frequency from the zero-Doppler time spacing, Sentinel-1 looks right with it from the azimuth
time interval.

## Reading a remote product without transferring it

A NISAR RSLC is tens of gigabytes, nearly all of it image samples. Its metadata is tens of kilobytes
lying within the first few megabytes.

A product is read through [ChunkManifests.jl](https://github.com/alex-s-gardner/ChunkManifests.jl),
which scans its metadata into a Zarr group with byte-range requests and then reads the samples a chunk
at a time as they are indexed. A remote product is read in place: opening it transfers the metadata,
and a window of its samples the chunks that window covers.

```julia
s = open_slc(RemoteHTTP(url))                  # Earthdata Login from ~/.netrc
s = open_slc(RemoteS3("s3://bucket/key.h5"))   # needs `using AWSS3`
pixels(s)[20001:22048, 30001:32048]            # reads only the chunks it needs
```

Earthdata URLs authenticate with an Earthdata Login bearer token: `RemoteHTTP(url; token)`,
`EARTHDATA_TOKEN`, or — with neither — one obtained from the `machine urs.earthdata.nasa.gov` entry of
`~/.netrc`, so that entry is the only setup. The host answers with a redirect to a pre-signed URL, and the
token is not sent on to it.

## Feeding a geometry kernel

This package reads products and says nothing about geometry. Converting an acquisition into a geometry
package's types belongs to that package, since they are its types:
[ImagePairGeometry.jl](https://github.com/alex-s-gardner/ImagePairGeometry.jl) extends its own
constructors over a `SLC`, so loading both is all that is needed.

```julia
using ImagePairGeometry, SLCDatasets

pair = CoregisteredPair(open_slc(url1), open_slc(url2))
pair.coordinate    # a RadarCoordinate, incidence angle included
pair.dt / 86400    # the repeat interval in days
```

What does live here is what a pair of *products* means: `repeat_interval` spans both epochs, and
`epoch_offset` gives the seconds between a product's epoch and midnight — the constant a kernel indexing
azimuth lines against midnight needs.

Verified against the reference implementations on real granules. For NISAR, against `h5py` and `isce3`:
all 26 geometry and identification values agree bitwise, and the scene-center range, azimuth time,
position and velocity agree to 0 ULP. For Sentinel-1, against `isce3` + `s1reader` on two granules from
different missions — every range, spacing, wavelength, PRF and image dimension of the mosaic and of all
27 bursts agrees bitwise, state vectors agree bitwise, and azimuth times agree to a nanosecond.

Those comparisons are part of the test suite and need no granule: the annotation fields and state
vectors each golden value was computed from are committed alongside it, and the products are rebuilt
from them. `test/reference/dump_sentinel1.py` regenerates the golden values from a granule and
`test/reference/extract_sentinel1_inputs.jl` the inputs; setting `SLCDATASETS_S1_DIR` to a directory
holding the granules checks against those directly instead.

## Sentinel-1

An IW SLC is read from a `.SAFE` directory or the zip of one. The state vectors come from a POEORB or
RESORB `.EOF` file, which the product does not contain, so `orbit` is required — an acquisition whose
orbit is unavailable is refused rather than returned half-read.

```julia
safe = "S1A_IW_SLC__1SSH_20151120T080202_....zip"
eof  = "S1A_OPER_AUX_POEORB_OPOD_....EOF"

mosaic = open_slc(safe; orbit = eof)                       # across subswaths 1-3
burst  = open_slc(safe; orbit = eof, swath = 2, burst = 3)  # one burst
nbursts(safe; swath = 2)                                    # 9
```

A TOPS product is three subswaths of bursts at different slant ranges rather than one image, so there
are two defensible geometries and both are available: the mosaic spanning the subswaths, and a single
burst. Reading a geometry touches the annotation XML alone, so the measurement rasters cost nothing and
a zip need not be unpacked; the samples are read only when asked for.

To work through a subswath's bursts, `bursts` parses the annotation once and returns them as an
`AbstractVector` of `SLC`s, where a call to `open_slc` per burst would re-read the product each time.

```julia
b = bursts(safe; orbit = eof, swath = 2)   # an AbstractVector{SLC}, one parse
length(b)                                  # 9
b[3].geometry.starting_range               # built on indexing
[nlines(s) for s in b]
```

### One image from a subswath's bursts

`merge_bursts` folds a vector of bursts into a single `SLC`. Whatever vector you have is the selection,
so a whole subswath or a slice of one merges directly:

```julia
b = bursts(safe; orbit = eof, swath = 2)
ref = merge_bursts(b)                       # across the subswath
ref = merge_bursts(b[3:7])                  # five of its bursts

amplitude(ref)[1:512, 1:512]                # this window, and no more of the product
validmask(ref)                              # which samples were imaged
```

The bursts are placed on one uniform azimuth time axis, with the overlap between neighbours divided at
its midpoint rather than repeated. A line's azimuth time is then `sensing_start + line / prf`
throughout, which is what `RadarGeometry` means by `prf`. Stacking whole bursts instead — what a
`.SAFE` measurement raster holds, and what
[burst2safe](https://github.com/ASFHyP3/burst2safe) writes — breaks that relation at every seam, by
about 166 lines each and some 1330 lines by the last. The merged image is correspondingly shorter:
12244 lines against 13572 stacked, for the granule the tests run on.

Placing a burst is a copy rather than a resampling. A burst's start lands on the shared axis within
2.2e-4 of a line on every granule measured, so the offset is a whole number of rows; the residual is
checked rather than assumed. The azimuth arithmetic reproduces hyp3-autorift's `merge_bursts_in_swath`
placement for placement, with one deliberate difference: the reference slices range samples with an
exclusive end against an inclusive index and so drops the last valid sample of every burst, which this
keeps.

Only each burst's valid region is placed. The margin around it is not zeroed in the file, so a reader
that copied whole bursts would admit samples the processor never imaged; outside every region the image
reads zero and `validmask` is `false`. A sample that *was* imaged can be zero too, which is why the mask
is what a correlator should be given rather than a test against zero.

Samples are read as they are indexed, so a subswath costs the windows taken from it and not its size.
Index with ranges where the shape of the read allows — a window is one read per burst it spans, about
four times faster per sample than resolving a row at a time.

A zip stores each raster deflated, so a line is not addressable without inflating everything before it.
A zipped product — local, over HTTP(S) or in S3 — therefore has its rasters inflated whole into a directory
(`Sentinel1Product(zip; orbit, dir)`, a temporary one by default) the first time its samples are read:
every subswath it was opened for at once, so a remote archive is read through several streams and each
member's bytes cross the network once. Reading a zipped product's metadata reads only the annotations.

### Bursts from ASF

ASF's burst extractor serves an SLC as its individual bursts, so a few bursts of a subswath can be had
without the several gigabytes of the whole product:

```julia
b = asf_bursts("S1A_IW_SLC__1SSH_20151120T080202_...", 2, "HH", 3:5; orbit = eof, dir = "bursts")
ref = merge_bursts(b)
```

Bursts are numbered from 1 here as everywhere else in this package, and the conversion to the
extractor's 0-based URLs happens inside; `AsfBurst(url)` reads a burst from the URL ASF's search lists
for it. Authentication is from `~/.netrc`. The extractor builds a burst on first request and answers
`202` until it is ready, which is retried. A burst's metadata — its parent product's whole annotation —
is held in memory, and its raster downloaded into `dir`, where one already there is not fetched again.

A merge of ASF bursts and a merge of the same bursts from a `.SAFE` give the same array — verified
sample for sample, including across seams — so which delivery format you have does not change what you
read. What differs is that each ASF file holds one burst, where a `.SAFE` raster holds a subswath's
bursts stacked.

A set of bursts, from one or more subswaths and from consecutive products of one datatake, also opens as
the product [burst2safe](https://github.com/ASFHyP3/burst2safe) assembles from them:

```julia
bs = [AsfBurst(url; dir = "bursts") for url in urls]
p = Sentinel1Product(bs; orbit = eof)      # the subswaths the bursts lie in
p.path                                     # the SAFE's name, e.g. S1A_IW_SLC__1SSV_..._069D69_0000
```

Its annotations and stacked rasters equal those of the SAFE `burst2safe` writes from the same bursts.

### A `Sentinel1Product`

The parsed product itself, when several subswaths are wanted from one read:

```julia
p = Sentinel1Product(safe; orbit = eof)    # subswaths 1-3, parsed once
nbursts(p, 2)
bursts(p, 3)
```

Times are reported against an `epoch` two days before the anchor burst's sensing start, matching the
offset `s1reader` and ISCE3's Doppler LUTs use; `epoch + sensing_start` is the instant the product
records.

## Scope

SLCs in radar geometry: NISAR-format HDF5 (RSLC) and Sentinel-1 IW SLC. Metadata for both; samples for a
merged Sentinel-1 subswath, from an unpacked `.SAFE`.

Mosaicking across subswaths is metadata only. IW1–3 lie at different slant ranges, so one
`starting_range` and one `prf` cannot describe their samples without resampling them onto a common grid,
and a merged image is per-subswath for that reason.

Other SAR products are out of scope rather than unimplemented. An interferogram (RIFG, RUNW) or a
covariance product (GCOV) carries a multilooked grid, and a geocoded SLC (GSLC) carries map coordinates
rather than a slant-range axis. None is described by `RadarGeometry`, which is the type this package
exists to produce, so a geocoded product is refused by name with that reason.
