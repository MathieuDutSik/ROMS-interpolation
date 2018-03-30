function TheRes=DATE_ConvertVect2mjd(TheVect)

eYear=TheVect(1,1);
eMonth=TheVect(1,2);
eDay=TheVect(1,3);
eHour=TheVect(1,4);
eMin=TheVect(1,5);
eSec=TheVect(1,6);

TheRes=persodate2mjd(eYear, eMonth, eDay, eHour, eMin, eSec);
