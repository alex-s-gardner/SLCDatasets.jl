# The one time representation both sensors report against, and the parsing of the timestamp forms the
# products use.
#
# A product records times to microseconds (Sentinel-1) or nanoseconds (NISAR), both finer than
# `DateTime`'s milliseconds. So an instant is carried as a `DateTime` truncated to the whole second plus
# a `Float64` residual, and differences are taken as the exact whole-second difference plus the
# residuals. Nothing is rounded into the sub-second digits on the way through.

"""
    UtcTime(datetime, seconds)

A UTC instant to better than millisecond resolution: `datetime` truncated to the whole second, plus
`seconds` of residual.

The residual is normally in `[0, 1)`, but is allowed past one second so a time can be advanced by an
arbitrary offset without renormalizing; [`seconds_between`](@ref) and `_utc_string` both handle that.
"""
struct UtcTime
    datetime::DateTime
    seconds::Float64
end

"""
    seconds_between(from::UtcTime, to::UtcTime) -> Float64

Seconds from one UTC instant to another.

The whole-second difference is exact in `Int64` milliseconds and the sub-second parts are added
afterwards, so the result does not inherit the rounding a `DateTime`-only subtraction would introduce.
"""
function seconds_between(from::UtcTime, to::UtcTime)
    whole = Dates.value(to.datetime - from.datetime) / 1000
    return whole + (to.seconds - from.seconds)
end

Base.isless(a::UtcTime, b::UtcTime) =
    a.datetime == b.datetime ? a.seconds < b.seconds : a.datetime < b.datetime

Base.:(==)(a::UtcTime, b::UtcTime) = a.datetime == b.datetime && a.seconds == b.seconds

_advance(t::UtcTime, seconds::Real) = UtcTime(t.datetime, t.seconds + seconds)

function _digit(b, i, s)
    c = b[i]
    (UInt8('0') <= c <= UInt8('9')) || _bad_utc(s)
    return Int(c - UInt8('0'))
end

_isdigit(c::UInt8) = UInt8('0') <= c <= UInt8('9')

@noinline _bad_utc(s) = throw(ArgumentError(
    "cannot read a UTC time from \"$s\"; expected `YYYY-MM-DDThh:mm:ss[.ffffff]`"))

"""
    parse_utc(s) -> UtcTime

The instant an ISO-8601 timestamp denotes, to the precision it is written with.

`YYYY-MM-DDThh:mm:ss` with an optional fractional part and an optional trailing zone designator, which
is ignored: these products are UTC throughout. The fractional digits are kept in the residual rather
than folded into the `DateTime`, so none are lost.
"""
function parse_utc(s::AbstractString)
    b = codeunits(s)
    n = length(b)
    # The separators are fixed by the format, so a wrong one means this is not a timestamp at all.
    (n >= 19 && b[5] == UInt8('-') && b[8] == UInt8('-') &&
     b[11] in (UInt8('T'), UInt8(' ')) && b[14] == UInt8(':') && b[17] == UInt8(':')) || _bad_utc(s)
    year = 1000 * _digit(b, 1, s) + 100 * _digit(b, 2, s) + 10 * _digit(b, 3, s) + _digit(b, 4, s)
    month = 10 * _digit(b, 6, s) + _digit(b, 7, s)
    day = 10 * _digit(b, 9, s) + _digit(b, 10, s)
    hour = 10 * _digit(b, 12, s) + _digit(b, 13, s)
    minute = 10 * _digit(b, 15, s) + _digit(b, 16, s)
    second = 10 * _digit(b, 18, s) + _digit(b, 19, s)
    frac = 0.0
    if n >= 21 && b[20] == UInt8('.')
        scale = 0.1
        i = 21
        while i <= n && _isdigit(b[i])
            frac += _digit(b, i, s) * scale
            scale /= 10
            i += 1
        end
    end
    return UtcTime(DateTime(year, month, day, hour, minute, second), frac)
end

"""
    _utc_string(t::UtcTime) -> String

`t` as `YYYY-MM-DDThh:mm:ss.ffffff`, the form the products write.

A residual carried past one second reads as the time it denotes rather than as an out-of-range field.
"""
function _utc_string(t::UtcTime)
    extra = floor(Int, t.seconds)
    frac = t.seconds - extra
    d = t.datetime + Dates.Second(extra)
    return _zpad(Dates.year(d), 4) * "-" * _zpad(Dates.month(d), 2) * "-" * _zpad(Dates.day(d), 2) * "T" *
           _zpad(Dates.hour(d), 2) * ":" * _zpad(Dates.minute(d), 2) * ":" * _zpad(Dates.second(d), 2) * "." *
           _zpad(round(Int, frac * 1_000_000), 6)
end

# `n` in decimal, zero-padded to `width` digits. Written out because `lpad` measures its padding with
# `textwidth`.
function _zpad(n::Int, width::Int)
    s = string(n)
    return ncodeunits(s) >= width ? s : String(fill(UInt8('0'), width - ncodeunits(s))) * s
end

# A trailing zone designator, if any, starts here. These products are UTC and `DateTime` carries no
# zone, so the designator is dropped rather than applied.
function _strip_zone(s::AbstractString)
    t = strip(s)
    isempty(t) && return SubString(t, 1, 0)
    (t[end] == 'Z' || t[end] == 'z') && return SubString(t, 1, prevind(t, lastindex(t)))
    # A `+hh:mm`/`-hhmm` offset, distinguished from the date's own hyphens by lying past the time.
    k = lastindex(t)
    while k >= firstindex(t) && !(t[k] == '+' || t[k] == '-')
        k = prevind(t, k)
    end
    k > 11 && return SubString(t, 1, prevind(t, k))
    return SubString(t, 1, lastindex(t))
end

"""
    parse_cf_epoch(units) -> DateTime

The reference instant of a CF-convention units string such as
`"seconds since 2025-10-28T00:00:00"`.

Only seconds are accepted: the products this reads record times in seconds, and silently rescaling a
different unit would corrupt every time by a constant factor.
"""
function parse_cf_epoch(units::AbstractString)
    s = _string(units)
    m = match(r"^\s*(\w+)\s+since\s+(.+?)\s*$", s)
    m === nothing && throw(ArgumentError(
        "cannot read a reference epoch from the units string \"$s\"; expected " *
        "\"<unit> since <timestamp>\""))
    unit = lowercase(m[1])
    unit in ("second", "seconds") || throw(ArgumentError(
        "time units are \"$unit\" in \"$s\", but only seconds are supported"))
    return DateTime(_strip_zone(replace(m[2], " " => "T")))
end

# Products record `zeroDopplerStartTime` to nanoseconds, which is finer than `DateTime`'s
# milliseconds. `Identification` keeps the string so nothing is lost; this converts only where
# millisecond resolution is enough, dropping the extra digits rather than rounding into them.
function _truncated_datetime(s::AbstractString)
    t = _strip_zone(_string(s))
    dot = findlast('.', t)
    dot === nothing && return DateTime(t)
    # Keep at most three fractional digits; `DateTime` rejects more.
    stop = min(lastindex(t), dot + 3)
    return DateTime(SubString(t, 1, stop))
end
