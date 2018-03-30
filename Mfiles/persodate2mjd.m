function mjd = persodate2mjd(varargin)
%DATE2MJD Modified Julian day number.
%
%   MJD = DATE2MJD(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the
%   modified Julian day number of the given date plus a fractional part
%   depending on the day and time of day.
%
%   Any missing MONTH or DAY will be replaced by ones.  Any missing HOUR,
%   MINUTE or SECOND will be replaced by zeros.
%
%   Beginning of time is  1968-05-23
%
%   This is necessary for fitting in the ROMS model.

   mjd = date2mjd(varargin{:}) - date2mjd(1968, 5, 23, 0, 0, 0);
