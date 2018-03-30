function str=DATE_ConvertSix2mystringPres(...
    eYear, eMonth, eDay, eHour, eMin, eSec)

strYear=StringNumber(eYear, 4);
strMonth=StringNumber(eMonth, 2);
strDay=StringNumber(eDay, 2);
strHour=StringNumber(eHour, 2);
strMin=StringNumber(eMin, 2);
strSec=StringNumber(eSec, 2);
str=[strYear '-' strMonth '-' strDay ...
     ' ' strHour ':' strMin ':' strSec];
