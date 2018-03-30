function monthlen=GetMonthLength(iYear, iMonth)
ListMonthLen=[31,0,31,  30,31,30, ...
	      31,31,30, 31,30,31];

if (abs(iMonth - 2)>0)
  monthlen=ListMonthLen(1,iMonth);
else
  if (isint(iYear/4) == 0)
    monthlen=28;
  else
    if (isint(iYear/100) == 0)
      monthlen=29;
    else
      if (isint(iYear/400) == 1)
	monthlen=29;
      else
	monthlen=28;
      end;
    end;
  end;
end;
