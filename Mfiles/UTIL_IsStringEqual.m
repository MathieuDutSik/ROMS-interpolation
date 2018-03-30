function test=UTIL_IsStringEqual(str1, str2)

nbChar1=size(str1, 2);
nbChar2=size(str2, 2);

if (nbChar1 == nbChar2)
  tes=str1==str2;
  Tsum=sum(tes(:));
  if (Tsum == nbChar1)
    test=1;
  else
    test=0;
  end;
else
  test=0;
end;
