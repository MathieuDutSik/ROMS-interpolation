function [test, reason]=IsMakeableFile(TheFile)

res=SplitLine(TheFile, '/');
nbComp=size(res,2);
if (nbComp == 1)
  test=1;
  reason='all ok';
else
  str='';
  len=nbComp-1;
  for i=1:len
    str=[str res{i} '/'];
  end;
  if (IsExistingFile(str) == 1)
    test=1;
    reason='all ok';
  else
    test=0;
    reason='directory non-existent';
  end;
end;
