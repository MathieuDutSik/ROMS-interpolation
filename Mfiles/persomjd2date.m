function [iYear, iMonth, iDay, iHour, iMin, iSec]=persomjd2date(TheVal)
[a, b, c, d, e, f]=mjd2date(TheVal+date2mjd(1968, 5, 23, 0, 0, 0));
[iYear, iMonth, iDay, iHour, iMin, iSec]=DateRoundingOperation(...
    a, b, c, d, e, f);
