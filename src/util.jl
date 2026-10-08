# Every scalar string in a NISAR product is stored as fixed-length bytes, which Zarr.jl surfaces as a
# fixed-width string type, or as a character for a one-byte string. Trailing NULs are padding.
#
# The return type is annotated because the callers pass the result of a read whose type is only known
# once the product is open: without it every string field of a product would infer as `Any` and each
# use of one would dispatch at runtime.
_string(x::String)::String = String(rstrip(x, '\0'))
_string(x::AbstractString)::String = _string(String(x))
_string(x::AbstractVector{UInt8})::String = _string(String(x))
_string(x::AbstractChar)::String = _string(string(Char(x)))
_string(x::AbstractArray)::String = _string(only(x))
