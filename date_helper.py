from datetime import datetime, timedelta

from constants import T0

class DateHelper:
    def gregorian_to_julian_date(self, year, month, day, hour=0, minute=0, second=0):
        """Convert Gregorian date to Julian Date."""
        if month <= 2:
            year -= 1
            month += 12

        A = year // 100
        B = 2 - A + (A // 4)

        JD = (365.25 * (year + 4716)) // 1 + (30.6001 * (month + 1)) // 1 + day + B - 1524.5
        JD += (hour + minute / 60 + second / 3600) / 24

        return JD

    def julian_date_to_gregorian(self, jd):
        """Convert Julian Date to Gregorian Date (UTC)."""
        UNIX_EPOCH_JD = 2440587.5  # Julian Date for 1970-01-01 00:00 UTC

        # Convert JD to Unix timestamp
        unix_seconds = (jd - UNIX_EPOCH_JD) * 86400.0
        
        # Convert to UTC datetime
        dt = datetime(1970, 1, 1) + timedelta(seconds=unix_seconds)
        
        year = dt.year
        month = dt.month
        date = dt.day
        hour = dt.hour
        minute = dt.minute
        second = dt.second + dt.microsecond / 1_000_000
        
        return year, month, date, hour, minute, second

    def add_integration_time(self, t):
        jd = T0 + (t / 86400)
        
        return jd
    