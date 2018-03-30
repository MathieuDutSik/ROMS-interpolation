function [ListTimeHistory, ListIFile, ListIRecord, nbFile, ...
	  ListAbsoluteIndex]=...
    ROMShistoryGetInfo2(PrefixHis, BeginTimeSec, EndTimeSec)

% It reads from an existing list of record
%
iFileBegin=0;
while(1)
  iFileBegin=iFileBegin+1;
  TheHisFile=[PrefixHis StringNumber(iFileBegin, 4) '.nc'];
  if (IsExistingFile(TheHisFile) == 1)
    break;
  end;
  if (iFileBegin == 9999)
    disp(['maybe you specified wrong PrefixHis']);
    disp([' there is no files ' PrefixHis '????.nc']);
    keyboard;
  end;
end;
iFileEnd=iFileBegin;
while(1)
  TheHisFile=[PrefixHis StringNumber(iFileEnd+1, 4) '.nc'];
  if (IsExistingFile(TheHisFile) == 0)
    break;
  end;
  iFileEnd=iFileEnd+1;
end;
%disp(['Beginning index=' num2str(iFileBegin) ...
%      '  ending index=' num2str(iFileEnd)]);


TheHisFile=[PrefixHis StringNumber(iFileBegin, 4) '.nc'];
nc=netcdf(TheHisFile, 'nowrite');
LTime=nc{'ocean_time'}(:);
nbRec=size(LTime,1);
if (iFileBegin == 1)
  nbPerFile=nbRec-1;
else
  nbPerFile=nbRec;
end;
TheBeginTimeHis=LTime(1,1);
DeltaTime=LTime(2,1)-LTime(1,1);
close(nc);
num2str(['DeltaTime = ' num2str(DeltaTime)]);


ListTimeHistory=zeros(0,1);
ListIFile=zeros(0, 1);
ListIRecord=zeros(0, 1);
ListAbsoluteIndex=zeros(0, 1);
currentTime=LTime(1,1);
idx=0;
iAbsIndex=0;
for iFile=iFileBegin:iFileEnd
  if (iFile == 1)
    nbRecord=nbPerFile+1;
  else
    nbRecord=nbPerFile;
  end;
  for iRecord=1:nbRecord;
    iAbsIndex=iAbsIndex+1;
    if (currentTime >= BeginTimeSec && currentTime <= EndTimeSec)
      idx=idx+1;
      ListIFile(idx,1)=iFile;
      ListIRecord(idx,1)=iRecord;
      ListTimeHistory(idx,1)=currentTime;
      ListAbsoluteIndex(idx,1)=iAbsIndex;
    end;
    currentTime=currentTime+DeltaTime;
  end;
end;
nbFile=1+iFileEnd-iFileBegin;
