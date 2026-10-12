# Strings read from a product — attributes, and the scalar strings of a NISAR product, stored as
# fixed-length bytes — may carry trailing NULs as padding.
#
# The return type is annotated because some callers pass a value whose type is only known once the
# product is read: without it a string field would infer as `Any` and each use of one would dispatch at
# runtime.
_string(x::String)::String = String(rstrip(x, '\0'))
_string(x::AbstractString)::String = _string(String(x))
_string(x::AbstractVector{UInt8})::String = _string(String(x))
_string(x::AbstractChar)::String = _string(string(Char(x)))
_string(x::AbstractArray)::String = _string(only(x))
