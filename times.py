# Name:        times
# Purpose: several time conversions and stuff, formerly found in
#          time_array_conversion, timestamp or my_globals
#
# Author:      Thomas Pfaff, Magdalena Eder, Dirk Schlabing
#
# Created:     17.11.2011
# Copyright:   (c) guttenberg 2011
# Licence:     <who cares?>
#!/usr/bin/env python
"""
Time helpers and conversions (:mod:`times`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Helper functions to convert arrays containing time information.
All the functions with a ``2`` in the name provide the advertised conversions
for scalar as well as for *array* input.

.. currentmodule:: vg.times

.. autosummary::
   :nosignatures:
   :toctree: generated/

    cwr2datetime
    cwr2str
    cwr2unix
    datetime2doy
    datetime2ordinal
    datetime2str
    datetime2unix
    doy2datetime
    iso2datetime
    iso2unix
    ordinal2datetime
    str2datetime
    str2ordinal
    str2unix
    unix2cwr
    unix2datetime
    unix2ordinal
    unix2str

    time_part
    time_part_sort
    timestamp2index
    index2timestamp
"""

from datetime import datetime, timedelta
import itertools
import time
import calendar
import numpy as np
##from vg import helpers as my


def cwr2datetime(cwr_string):
    """Converts a CWR-timesting into a datetime object.
    >>> cwr2datetime("2040001.00")
    datetime.datetime(2040, 1, 1, 0, 0)
    """
    return unix2datetime(cwr2unix(cwr_string))


def cwr2str(cwr_string, d_format="%d.%m.%Y %H:%M:%S"):
    """Converts a CWR-timestring into a human-readable timestring.

    >>> cwr2str("2040001.00")
    '01.01.2040 00:00:00'
    """
    dt = cwr2datetime(cwr_string)
    return datetime2str(dt, d_format)


def cwr2unix(cwr_string):
    """Converts a CWR-timestring to a unix timestamp.

    >>> cwr2unix("2040001.00")
    2208988800.0
    """
    def cwr2unix_single(cwr_string_single):
##        cwr_string_single = str(cwr_string_single)
        cwr_string_single = '%013.5f' % float(cwr_string_single)
        yearday = int(float(cwr_string_single[4:]))  # round down to whole days
        # this saves us from calculating the month and day-of-month
        timestamp = str2unix(cwr_string_single[:7], "%Y%j")

        seconds = (float(cwr_string_single[4:]) - yearday) * 86400
        return timestamp + seconds

    try:
        return cwr2unix_single(cwr_string)
    except (ValueError, TypeError):
        return np.array([cwr2unix_single(cwr) for cwr in cwr_string])


def datetime2ordinal(dt):
    """Converts a datetime object into an ordinal.

    >>> import datetime
    >>> dt = datetime.datetime(2040, 1, 1, 0, 0, 0, 500000)
    >>> datetime2ordinal(dt)
    744730
    >>> datetime2ordinal(datetime.datetime(1, 1, 1, 0, 0))
    1
    """
    try:
        return np.array([sub_dt.toordinal() for sub_dt in dt])
    except TypeError:
        return dt.toordinal()


def datetime2cwr(dt):
    """Converts a datetime object into a CWR/Julian timestamp.

    >>> import datetime
    >>> datetime2cwr(datetime.datetime(2040, 1, 1, 0, 0, 0, 0))
    2040001.0
    """
    return unix2cwr(datetime2unix(dt))


def datetime2str(dt, d_format="%d.%m.%Y %H:%M:%S"):
    """Converts a datetime object into a human-readable string.

    >>> import datetime
    >>> datetime2str(datetime.datetime(2040, 1, 1, 0, 0, 0, 500000))
    '01.01.2040 00:00:00'
    """
    try:
        return np.array([sub_datetime.strftime(d_format)
                         for sub_datetime in dt])
    except TypeError:
        return dt.strftime(d_format)


def datetime2unix(dt):
    """Converts a dateime object into a unix-timestamp.

    >>> import datetime
    >>> dt = datetime.datetime(2040, 1, 1, 0, 0, 0, 500000)
    >>> datetime2unix(dt)
    2208988800.5
    >>> datetime2unix(datetime.datetime(1970, 1, 1, 0, 0))
    0.0
    """
    try:
        tzinfo = dt.tzinfo if hasattr(dt, "tzinfo") else dt[0].tzinfo
    except IndexError:
        tzinfo = None
    diff = dt - datetime(1970, 1, 1, 0, 0, tzinfo=tzinfo)
    try:
        return np.array([sub_diff.days * 86400 + sub_diff.seconds +
                         sub_diff.microseconds / 1e6
                         for sub_diff in diff])
    except TypeError:
        return diff.days * 86400 + diff.seconds + diff.microseconds / 1e6


def datetime2doy(dt):
    """ Extracts the day of year as a float from the given datetimes.

    >>> import datetime
    >>> dt = datetime.datetime(2040, 1, 1, 12, 30, 30, 500000)
    >>> datetime2doy(dt)
    1.5211863425925927
    """
    def datetime2doy_single(dt):
        return (time_part(dt, "%j") + dt.hour / 24. + dt.minute / (24. * 60) +
                dt.second / (24. * 60 ** 2) +
                dt.microsecond / (24. * 60 ** 2 * 1e6))
    try:
        return np.array([datetime2doy_single(sub_dt) for sub_dt in dt])
    except TypeError:
        return datetime2doy_single(dt)


def doy2datetime(doy, year=2000):
    """Constructs a datetime object from given day of year.

    Parameters
    ----------
    year : int, optional
        The year of the day of the year.

    Examples
    --------
    >>> import datetime
    >>> doy2datetime(1.5211863425925927, 2040)
    datetime.datetime(2040, 1, 1, 12, 30, 30, 500000)"""
    def doy2datetime_single(doy, year):
        return datetime(year, month=1, day=1) + timedelta(days=(doy - 1))
    # TODO: interpret 365, 0 as a new year
    try:
        if type(year) is int:
            years = itertools.cycle([year])
        else:
            years = year
        return np.array([doy2datetime_single(sub_doy, sub_year)
                         for sub_doy, sub_year in zip(doy, years)])
    except TypeError:
        return doy2datetime_single(doy, year)


# these two are virtually the same

def iso2datetime(iso):
    """Converts an ISO formated time string to a datetime object."""
    try:
        return str2datetime(iso, "%Y-%m-%dT%H:%M:%S.%f")
    except (ValueError, TypeError):
        return str2datetime(iso, "%Y-%m-%dT%H:%M:%S")


def datetimefromisoformat(ts, fmt='%Y-%m-%dT%H:%M:%S'):
    return datetime.strptime(ts, fmt)
# (same same but different)


def iso2unix(iso):
    """Converts an ISO formated time string to a unix time stamp."""
    try:
        return str2unix(iso, "%Y-%m-%dT%H:%M:%S.%f")
    except ValueError:
        # sometimes the microsecond is missing
        return str2unix(iso, "%Y-%m-%dT%H:%M:%S")


def ordinal2datetime(ord_):
    """Converts an ordinal to a datetime object.

    >>> ordinal2datetime(744730)
    datetime.datetime(2040, 1, 1, 0, 0)
    """
    try:
        return np.array([datetime.fromordinal(sub_ord) for sub_ord in ord_])
    except TypeError:
        return datetime.fromordinal(ord_)


def str2datetime(str_, d_format="%d.%m.%Y %H:%M:%S"):
    """Converts a human readable time string tinto a datetime object.

    >>> str2datetime("01.01.2040 00:00:00")
    datetime.datetime(2040, 1, 1, 0, 0)
    """
    try:
        return np.array([datetime.strptime(sub_str, d_format)
                         for sub_str in str_])
    except (TypeError, ValueError):
        return datetime.strptime(str_, d_format)


def str2ordinal(str_, d_format="%d.%m.%Y %H:%M:%S"):
    """Converts a human readable time string to an ordinal.

    >>> str2ordinal("01.01.2040 00:00:00")
    744730
    >>> str2ordinal("01.01.0001 00:00:00")
    1
    """
    return datetime2ordinal(str2datetime(str_, d_format))


def str2unix(str_, d_format="%d.%m.%Y %H:%M:%S"):
    """Converts a time-string of the given d_format (default:
    "dd.mm.yyyy HH:MM:SS") into a unix-timestamp.

    >>> str2unix("01.01.2040 00:00:00")
    2208988800.0
    >>> str2unix("01.01.1970 00:00:00")
    0.0
    """
    return datetime2unix(str2datetime(str_, d_format))


def unix2cwr(timestamp):
    """Converts a unix-timestamp to a CWR-timestring.

    >>> unix2cwr(2208988800.0)
    2040001.0
    """
    def unix2cwr_single(timestamp_single):
        src_tuple = time.gmtime(timestamp_single)
        # we are looking for "{year}{yearday}.{seconds expressed in days}
        year = str(src_tuple[0])
        yearday = time.strftime("%j", src_tuple)
        year_yearday_tuple = time.strptime(year + yearday, "%Y%j")
        diff_seconds = timestamp_single - calendar.timegm(year_yearday_tuple)
        seconds_as_days = float(diff_seconds) / 86400
        return "%s%03d%s" % (year, int(yearday),
                             str("%f" % seconds_as_days)[1:])

    try:
        return float(unix2cwr_single(timestamp))
    except TypeError:
        return np.array([float(unix2cwr_single(ts)) for ts in timestamp])


def unix2datetime(timestamp):
    """Convert a unix time stamp to a datetime object.

    >>> import datetime
    >>> unix2datetime(2208988800.5)
    datetime.datetime(2040, 1, 1, 0, 0, 0, 500000)
    >>> unix2datetime(0)
    datetime.datetime(1970, 1, 1, 0, 0)
    """
    try:
        return np.array([datetime(1970, 1, 1) + timedelta(seconds=float(stamp))
                         for stamp in timestamp])
    except TypeError:
        return datetime(1970, 1, 1) + timedelta(seconds=float(timestamp))


def unix2ordinal(timestamp):
    """Converts a unix time stamp to an ordinal.

    >>> unix2ordinal(2208988800.0)
    744730
    """
    try:
        np.array([datetime2ordinal(unix2datetime(stamp))
                  for stamp in timestamp])
    except TypeError:
        return datetime2ordinal(unix2datetime(timestamp))


def unix2str(timestamp, d_format="%d.%m.%Y %H:%M:%S"):
    """Converts a unix-timestamp into a time-string of the given d_format
    (default: "dd.mm.yyyy HH:MM:SS").

    >>> unix2str(2208988800.5)
    '01.01.2040 00:00:00'
    >>> unix2str(0)
    '01.01.1970 00:00:00'
    """
    return datetime2str(unix2datetime(timestamp), d_format)


def timedelta2seconds(dt):
    return dt.days * 86400 + dt.seconds + dt.microseconds * 1e-6


def timedelta2slice(delta, dt, offset=0, step=None):
    start = offset
    stop = None if delta is None else offset + int(timedelta2seconds(delta) /
                                                   timedelta2seconds(dt))
    return slice(start, stop, step)


def timestamps2slice(startts=None, endts=None, dt=None, refts=None, step=None):
    """Ask Thomas.

    >>> timestamps2slice('2010-03-01T23:54:00', '2010-03-07T23:52:00', \
                         timedelta(days=1))
    slice(0, 5, None)
    >>> timestamps2slice('2010-03-01T23:54:00', '2010-03-07T23:52:00', \
                         timedelta(days=1), '2010-03-04T00:00:00')
    slice(-2, 3, None)
    >>> timestamps2slice('2010-03-01T23:54:00', '2010-03-07T23:52:00', \
                         timedelta(days=1), '2010-01-04T00:00:00')
    slice(56, 61, None)
    """
    # Dirk: Ist das Kunst oder kann das weg?
#    if dt is None:
#        tdelta = timedelta(seconds=1)

    if startts is None:
        tstart = None
    if type(startts) == str:
        tstart = datetimefromisoformat(startts)
    else:
        tstart = startts

    if endts is None:
        tend = None
    if type(endts) == str:
        tend = datetimefromisoformat(endts)
    else:
        tend = endts

    if refts is None:
        tref = tstart
    elif type(refts) == str:
        tref = datetimefromisoformat(refts)
    else:
        tref = refts

    if (tstart is None) or (tref is None):
        slicestart = None
    else:
        slicestart = int(timedelta2seconds(tstart - tref) /
                         timedelta2seconds(dt))

    return timedelta2slice(None if ((tend is None)or(tstart is None))
                           else tend - tstart, dt, slicestart, step)


def timestamp2index(ts, dt, refts, **kwargs):
    """Calculates the array index for a certain time in an equidistant
    time-series given the reference time (where the index would be 0)
    and the time discretization.
    If any of the input parameters contains timezone information, all others
    also need to contain timezone information.

    Parameters
    ----------
    ts        : str or datetime-object
                The timestamp to determine the index for
                If it is a string, it will be converted to datetime using the
                function _datetimefromisoformat Formatting keywords may be
                passed to this function

    dt        : str or timedelta object
                The discretization of the time series (the amount of time that
                elapsed between indices)
                If used as a string, it needs to be given in the format
                "keyword1=value1,keyword2=value2". Keywords must be understood
                by the timedelta constructor (like days, hours,
                minutes, seconds) and the values may only be integers.

    refts     : str or datetime-object
                The timestamp to determine the index for
                If it is a string, it will be converted to datetime using the
                function _datetimefromisoformat Formatting keywords may be
                passed to this function

    Returns
    -------
    index    : integer
               The index of a discrete time series array of the given
               parameters

    Example
    -------
    >>> timestr1, timestr2 = '2008-06-01T00:00:00', '2007-01-01T00:00:00'
    >>> timestamp2index(timestr1, 'minutes=5', timestr2)
    148896
    >>> timestamp2index(timestr1, 'hours=1,minutes=5',timestr2)
    11453
    >>> timestamp2index(timestr1, timedelta(hours=1, minutes=5), timestr2)
    11453
    """
    if not isinstance(ts, datetime):
        _ts = datetimefromisoformat(ts, **kwargs)
    else:
        _ts = ts
    if not isinstance(refts, datetime):
        _refts = datetimefromisoformat(refts, **kwargs)
    else:
        _refts = refts
    if not isinstance(dt, timedelta):
        kwargs = dict([(sp[0], int(sp[1]))
                       for sp in [item.split('=') for item in dt.split(',')]])
        _dt = timedelta(**kwargs)
    else:
        _dt = dt
    return int(timedelta2seconds(_ts - _refts) / timedelta2seconds(_dt))


def isoformat2unix(ts):
    dt = timedelta(seconds=1)
    tstart = datetimefromisoformat('1970-01-01T00:00:00')
    return timestamp2index(datetimefromisoformat(ts), dt, tstart)


def index2timestamp(idx, dt, refts, **kwargs):
    """Calculates the ISOstring timestamp for a certain index in an equidistant
    time-series given the reference time (where the index would be 0)
    and the time discretization.
    If any of the input parameters contains timezone information, all others
    also need to contain timezone information.

    Parameters
    ----------
    idx        : str or datetime-object
                The timestamp to determine the index for
                If it is a string, it will be converted to datetime using the
                function _datetimefromisoformat Formatting keywords may be
                passed to this function

    dt        : str or timedelta object
                The discretization of the time series (the amount of time that
                elapsed between indices)
                If used as a string, it needs to be given in the format
                "keyword1=value1,keyword2=value2". Keywords must be understood
                by the timedelta constructor (like days, hours,
                minutes, seconds) and the values may only be integers.

    refts     : str or datetime-object
                The timestamp to determine the index for
                If it is a string, it will be converted to datetime using the
                function _datetimefromisoformat Formatting keywords may be
                passed to this function

    Returns
    -------
    ISO-timestamp    : string
               The ISO-timestamp of a discrete time series array of the given
               parameters in this format: '%Y-%m-%dT%H:%M:%S'

    Examples
    --------
    >>> index2timestamp(25637, 'seconds=10', '2010-09-25T00:00:10')
    '2010-09-27T23:13:00'
    >>> index2timestamp(365, 'days=1', '2010-09-25T00:00:10')
    '2011-09-25T00:00:10'
    """
    if not isinstance(refts, datetime):
        _refts = datetimefromisoformat(refts, **kwargs)
    else:
        _refts = refts
    if not isinstance(dt, timedelta):
        kwargs = dict([(sp[0], int(sp[1]))
                       for sp in [item.split('=') for item in dt.split(',')]])
        _dt = timedelta(**kwargs)
    else:
        _dt = dt

    return str(datetime.isoformat(_refts + idx * _dt))


# ### Time functions from Dirk's my_globals:

def build_diff_timestring(diff):
    """Split time difference given in seconds into a string representation of
    days, hours, minutes and seconds.
    """
    days = hours = minutes = 0
    seconds = abs(diff)
    while seconds >= 60:
        minutes += 1
        seconds -= 60
    while minutes >= 60:
        hours += 1
        minutes -= 60
    while hours >= 24:
        days += 1
        hours -= 24

    diff_str = ""
    if days > 0:
        diff_str += "%d d " % days
    if hours > 0:
        diff_str += "%d h " % hours
    if minutes > 0:
        diff_str += "%d m " % minutes
    diff_str += "%f s" % seconds

    return diff_str


def time_part(timestamps_or_datetimes, sub_format_str):
    """Returns the "time_part" of a timestamp or datetimes.
    The time_part has to be convertible into an integer.
    This is useful to group values according to months, weeks and so on.
    """
    def single_time_part_unix(stamp, d_format):
        return (int(unix2str(stamp, d_format))
                if np.isfinite(stamp) else None)

    def single_time_part_date(time_, d_format):
        return int(time_.strftime(d_format))

    times_ = np.asarray(timestamps_or_datetimes)

    if times_.dtype == object:
        single_time_part = single_time_part_date
    else:
        single_time_part = single_time_part_unix

    try:
        return np.array([single_time_part(timestamp, sub_format_str)
                         for timestamp in times_])
    except TypeError:
        return single_time_part(np.asscalar(times_), sub_format_str)


def time_part_(datetimes, date_part):
    try:
        return np.array([getattr(date, date_part) for date in datetimes])
    except TypeError:
        return getattr(datetimes, date_part)


def time_part_sort_(datetimes, values, date_part):
    assert len(datetimes) == len(values)
    time_of_values = np.array([getattr(date, date_part) for date in datetimes])
    grouped_values = []
    times = sorted(set(time_of_values))
    for time_key in times:
        grouped_values += \
            [values[np.where(time_of_values == time_key)[0]]]
    return times, grouped_values


def time_part_sort(timestamps, values, sub_format_str):
    """Groups "values" as a nested list according to the "sub_format_str" of
    the "timestamps"."""
    assert len(timestamps) == len(values)
    time_of_values = time_part(timestamps, sub_format_str)
    grouped_values = []
    times = sorted(set(time_of_values))
    for time_key in times:
        grouped_values += \
            [values[np.where(time_of_values == time_key)[0]]]
    return times, grouped_values


def expand_timeseries(timestamps, repeats=4, values=None):
    """Interpolates timestamps and repeats according values.  Assumes constant
    length of time-step.

    Examples
    --------
    >>> import numpy as np
    >>> timestamps = str2unix(["01.01.1980 00:00:00", "02.01.1980 00:00:00"])
    >>> values = np.array([1, 2])
    >>> timestamps, values = expand_timeseries(timestamps, 2, values)
    >>> print unix2str(timestamps), values
    ['01.01.1980 00:00:00' '01.01.1980 12:00:00' '02.01.1980 00:00:00'
     '02.01.1980 12:00:00'] [1 1 2 2]
    """
    dt = (timestamps[1] - timestamps[0]) / repeats
    dts = np.arange(0, repeats * dt, dt)
    old_len = len(timestamps)
    timestamps = timestamps.repeat(repeats).reshape((old_len, repeats))
    timestamps += dts
    timestamps = timestamps.ravel()
    if values is not None:
        return timestamps, values.repeat(repeats)
    return timestamps

def xls2datetime(xldate, datemode=0):
    """ Here's the bare-knuckle no-seat-belts use-at-own-risk version
    posted by John Machin in http://stackoverflow.com
    datemode: 0 for 1900-based, 1 for 1904-based
    """
    def xls2datetime_single(xldate, datemode):
        return (
            datetime(1899, 12, 30) #30.Dez, weil 1900 in Excel ein Schaltjahr ist --> ab dem 1. Maerz stimmts
            + timedelta(days=xldate + 1462 * datemode) )
    try:
        return xls2datetime_single(xldate, datemode)
    except (ValueError, TypeError):
        return np.array([xls2datetime_single(xl, datemode) for xl in xldate])


if __name__ == '__main__':
    import doctest
    doctest.testmod()
