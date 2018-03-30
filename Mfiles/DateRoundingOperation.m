function [a, b, c, d, e, f]=DateRoundingOperation(...
    a1, b1, c1, d1, e1, f1)


H=round(f1);

if (H < 60)
  a=a1;
  b=b1;
  c=c1;
  d=d1;
  e=e1;
  f=H;
else
  H=e1+1;
  f=0;
  if (H < 60)
    a=a1;
    b=b1;
    c=c1;
    d=d1;
    e=H;
  else
    H=d1+1;
    e=0;
    if (H < 24)
      a=a1;
      b=b1;
      c=c1;
      d=H;
    else
      H=c1+1;
      d=0;
      monthlen=GetMonthLength(a1, b1);
      if (H < monthlen+1)
	a=a1;
	b=b1;
	c=H;
      else
	H=b1+1;
	c=1;
	if (H < 12+1)
	  a=a1;
	  b=H;
	else
	  a=a1+1;
	  b=1;
	end;
      end;
    end;
  end;
end;
