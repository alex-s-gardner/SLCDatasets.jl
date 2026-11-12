# The time representation both sensors report against, and the parsing of the timestamp forms the
# products use.
#
# A product records times to microseconds (Sentinel-1) or nanoseconds (NISAR), both finer than
# `DateTime`'s milliseconds. An instant is a `NanoDates.NanoDate`, which keeps whole nanoseconds, and a
# difference is taken in integer nanoseconds before it becomes seconds, so nothing is rounded into the
# sub-second digits on the way through.

"""
    seconds_between(from::NanoDate, to::NanoDate) -> Float64

Seconds from one UTC instant to another, from their exact difference in nanoseconds.
"""
seconds_between(from::NanoDate, to::NanoDate) = Dates.value(to - from) / 1.0e9

# Seconds from a whole-second epoch, as a `RadarGeometry` or `StateVectors` gives one, to an instant.
seconds_between(epoch::DateTime, to::NanoDate) = seconds_between(NanoDate(epoch), to)

# `t` advanced by `seconds`, to the nanosecond.
_advance(t::NanoDate, seconds::Real) = t + Dates.Nanosecond(round(Int64, seconds * 1.0e9))

function _digit(b, i, s)
    c = b[i]
    (UInt8('0') <= c <= UInt8('9')) || _bad_utc(s)
    return Int(c - UInt8('0'))
end

_isdigit(c::UInt8) = UInt8('0') <= c <= UInt8('9')

@noinline _bad_utc(s) = throw(ArgumentError(
    "cannot read a UTC time from \"$s\"; expected `YYYY-MM-DDThh:mm:ss[.fffffffff]`"))

"""
    parse_utc(s) -> NanoDate

The instant an ISO-8601 timestamp denotes, to the nanosecond.

`YYYY-MM-DDThh:mm:ss` with an optional fractional part and an optional trailing zone designator, which
is ignored: these products are UTC throughout. Fractional digits past the ninth are dropped.
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
    nanos = 0
    if n >= 21 && b[20] == UInt8('.')
        scale = 100_000_000
        i = 21
        while i <= n && _isdigit(b[i])
            nanos += _digit(b, i, s) * scale
            scale ÷= 10
            i += 1
        end
    end
    return NanoDate(DateTime(year, month, day, hour, minute, second)) + Dates.Nanosecond(nanos)
end

# `t` to the whole second, which a product's epoch and its file names are given at.
_whole_second(t::NanoDate) = floor(DateTime(t), Dates.Second)

"""
    _utc_string(t::NanoDate) -> String

`t` as `YYYY-MM-DDThh:mm:ss.ffffff`, the form the products write.
"""
function _utc_string(t::NanoDate)
    d = t.datetime
    us = 1000 * Dates.millisecond(t) + Dates.microsecond(t)
    return _zpad(Dates.year(d), 4) * "-" * _zpad(Dates.month(d), 2) * "-" * _zpad(Dates.day(d), 2) * "T" *
           _zpad(Dates.hour(d), 2) * ":" * _zpad(Dates.minute(d), 2) * ":" * _zpad(Dates.second(d), 2) * "." *
           _zpad(us, 6)
end

# `n` in decimal, zero-padded to `width` digits. Written out because `lpad` measures its padding with
# `textwidth`.
function _zpad(n::Int, width::Int)
    s = string(n)
    return ncodeunits(s) >= width ? s : String(fill(UInt8('0'), width - ncodeunits(s))) * s
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
    return DateTime(parse_utc(m[2]))
end
