# Reaching an acquisition's samples, and the amplitude of them.
#
# Reading samples is not something every product this package opens can do: a Sentinel-1 product read
# from a zip cannot reach its raster at all — the entry is deflated, so a line is not addressable without
# inflating everything before it — and a product rebuilt from metadata alone has no sample array to read.
# So `pixels` is a backend method a backend may not have, and a backend without one says why rather than
# raising a `MethodError` naming an internal type.
#
# Amplitude is a view rather than a read. The reference pipeline correlates `abs` of the samples and
# discards the phase, but `AutoRIFT.jl` can also correlate the complex samples directly, so the phase is
# kept and dropped where it is asked for instead of at the read.

"""
    pixels(s::SLC) -> AbstractMatrix

The acquisition's samples, as complex values.

The result is read from the product as it is indexed rather than held, so a subswath costs the windows
taken from it. Index it with ranges where the shape of the read allows: `p[rows, cols]` is one read per
burst spanned, while scalar indexing resolves a row per sample.

Not every product can reach its samples. A Sentinel-1 product read from a zip cannot — the raster is
deflated inside it — and neither can a remote product whose image was never transferred; both say so
rather than returning something partial. A NISAR product reads its samples out of its HDF5 dataset a
window at a time; see [`SLCDatasets.NisarRaster`](@ref). [`amplitude`](@ref) for the real magnitudes, and
[`validmask`](@ref) for which samples were imaged.
"""
pixels(s::AbstractSLC) = read_pixels(s.backend)

"""
    validmask(s::SLC) -> AbstractMatrix{Bool}

Which of the acquisition's samples were imaged.

`false` where no data was recorded — the corners a merged subswath leaves, and the margin around each
burst's valid region. A sample that was imaged may still be zero, so this is the honest test for
absence rather than comparing the samples against zero.

Pass it alongside the samples to a consumer that takes a validity mask, rather than letting it derive
one: it is computed from the layout and reads nothing.
"""
validmask(s::AbstractSLC) = validmask(pixels(s))

@noinline _no_pixels(b::AbstractSLCBackend) = throw(ArgumentError(
    "a $(nameof(typeof(b))) reads metadata only, so this acquisition's samples are not available"))

"""
    read_pixels(b::AbstractSLCBackend) -> AbstractMatrix

The samples of the product `b` reads, for a backend that can reach them.

The optional half of the backend protocol: a backend that reads metadata alone leaves this
unimplemented, and [`pixels`](@ref) then reports that rather than failing on dispatch.
"""
read_pixels(b::AbstractSLCBackend) = _no_pixels(b)

"""
    burst_raster(b::AbstractBurstBackend) -> BurstRaster

One burst's samples and where in its raster the burst begins.

What a merge asks of each of its sources, so that a burst stacked with its neighbours in a subswath's
raster and a burst delivered as its own file are placed the same way. A backend that cannot reach its
samples reports that, as [`read_pixels`](@ref) does.
"""
burst_raster(b::AbstractBurstBackend) = BurstRaster(read_pixels(b))

"""
    burst_rasters(sources::AbstractVector{<:AbstractBurstBackend}) -> Vector{BurstRaster}

The samples of several bursts of one subswath, opening each file they lie in once.

Bursts delivered as separate files have nothing to share and are asked one at a time; bursts stacked in a
subswath's raster share it, so a backend whose bursts do gives a method that opens it once.
"""
burst_rasters(sources::AbstractVector{<:AbstractBurstBackend}) = map(burst_raster, sources)

"""
    Amplitude{T,P} <: AbstractMatrix{T}

The magnitudes of a complex array, taken as they are read.

Holds no samples of its own, so wrapping a whole subswath costs nothing and a window taken from it
costs the window. `eltype` is `Float32`: the samples are 16-bit integers, whose magnitude is exact in
it, and it is what a correlator works in.
"""
struct Amplitude{T,P<:AbstractArray} <: AbstractMatrix{T}
    parent::P
end

"""
    amplitude(A) -> Amplitude
    amplitude(s::SLC) -> Amplitude

The magnitude of each sample, computed as it is read.

This is what the reference pipeline correlates: it takes `abs` of the samples and discards the phase,
so a feature tracker sees magnitudes. The phase is still in the product — [`pixels`](@ref) returns it
— for a consumer that correlates the complex samples instead.

```julia
ref = merge_bursts(bursts(safe; orbit = eof, swath = 2))
a = amplitude(ref)      # nothing read yet
a[1:512, 1:512]         # this window, and no more of the product
```
"""
amplitude(A::AbstractArray{<:Complex}) = Amplitude{Float32,typeof(A)}(A)
amplitude(A::AbstractArray{<:Real}) = Amplitude{Float32,typeof(A)}(A)
amplitude(s::AbstractSLC) = amplitude(pixels(s))

Base.parent(a::Amplitude) = a.parent
Base.size(a::Amplitude) = size(a.parent)
Base.IndexStyle(::Type{<:Amplitude{<:Any,P}}) where {P} = IndexStyle(P)

# `abs` of a `Complex{Int16}` promotes to `Float64` and goes through `hypot`, which guards against an
# overflow that cannot happen here: a sample's parts are 16-bit, so their squares and the sum of them
# are exact in `Float32`. Taking the root directly is twice as fast and gives the same answer.
_magnitude(::Type{T}, z::Complex) where {T} = sqrt(T(real(z))^2 + T(imag(z))^2)
_magnitude(::Type{T}, x::Real) where {T} = T(abs(x))

Base.@propagate_inbounds Base.getindex(a::Amplitude{T}, i::Int) where {T} =
    _magnitude(T, a.parent[i])
Base.@propagate_inbounds Base.getindex(a::Amplitude{T}, i::Int, j::Int) where {T} =
    _magnitude(T, a.parent[i, j])

# A window of magnitudes, taken from a window of the parent rather than sample by sample: the parent's
# range indexing is what makes reading a merged subswath one read per burst.
Base.@propagate_inbounds function Base.getindex(a::Amplitude{T},
                                                rows::AbstractUnitRange{<:Integer},
                                                cols::AbstractUnitRange{<:Integer}) where {T}
    src = a.parent[rows, cols]
    out = similar(src, T)
    for i in eachindex(src, out)
        out[i] = _magnitude(T, src[i])
    end
    return out
end

validmask(a::Amplitude) = validmask(a.parent)

# **The parent's chunks, because taking a magnitude does not move a boundary.** `Amplitude` is what a
# consumer holds — `amplitude(pixels(product))` is the whole of how a correlator is handed a band — and
# without these it answers `DiskArrays`' generic fallback for an `AbstractArray`: `Unchunked()`, and one chunk
# covering the array. A caller wanting to align its windows to storage would then see no boundaries to align
# to, and the alignment is worth 1.2x to 1.5x on a filtered NISAR band.
#
# Forwarded rather than inherited: `Amplitude` also wraps arrays that are already in memory, and routing
# those through `DiskArrays`' machinery would add indirection to a read that is a `getindex` on a `Matrix`.
DiskArrays.haschunks(a::Amplitude) = DiskArrays.haschunks(a.parent)
DiskArrays.eachchunk(a::Amplitude) = DiskArrays.eachchunk(a.parent)
