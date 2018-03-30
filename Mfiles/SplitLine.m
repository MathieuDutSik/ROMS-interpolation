function result=SplitLine(A, pat)
%
% result=SplitLine(A, pat)
% This function splits the string A into blocks according
% to the pattern A much in the same way as the perl
% @A=split($str1, $str2);
% 
clear('result');
nbChar=size(A, 2);
nbCharPat=size(pat, 2);
if (nbChar < nbCharPat)
  result=zeros(0,0);
else
  nbSearch=1+nbChar-nbCharPat;
  LSM=zeros(0,1);
  ipos=0;
  for iChar=1:nbSearch
    TheSelect=A(1, iChar:iChar+nbCharPat-1);
    test1=TheSelect == pat;
    test2=sum(test1) == nbCharPat;
    if (test2 == 1)
      ipos=ipos+1;
      LSM(ipos, 1)=iChar;
    end;
  end;
  if (ipos == 0)
    result{1}=A;
  else
    nbREC=0;
    if (LSM(1,1) > 1)
      iBeg=1;
      iEnd=LSM(1,1)-1;
      TheSelect=A(1, iBeg:iEnd);
      nbREC=nbREC+1;
      result{nbREC}=TheSelect;
    end;
    for iter=2:ipos
      iBeg=LSM(iter-1,1)+nbCharPat;
      iEnd=LSM(iter,1)-1;
      TheSelect=A(1, iBeg:iEnd);
      nbREC=nbREC+1;
      result{nbREC}=TheSelect;
    end;
    if (LSM(ipos,1) < nbSearch)
      iBeg=LSM(ipos,1)+nbCharPat;
      iEnd=nbSearch;
      TheSelect=A(1, iBeg:iEnd);
      nbREC=nbREC+1;
      result{nbREC}=TheSelect;
    end;
  end;
end;
