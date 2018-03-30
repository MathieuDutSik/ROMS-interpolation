function str=DATE_ConvertMjd2mystringPres(TheVal)
[a, b, c, d, e, f]=persomjd2date(TheVal);
str=DATE_ConvertSix2mystringPres(a, b, c, d, e, f);
