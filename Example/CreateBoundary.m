BigGridFile='ncom2_376x136_TOP7samp_simple_r0_20.nc';
SmallGridFile='ROVINJ_700m_grd_V4.nc';


AddiRecordInfo.BigThetaS=3;
AddiRecordInfo.BigThetaB=0.35;
AddiRecordInfo.Nbig=20;
AddiRecordInfo.Bighc=3;
%
AddiRecordInfo.SmaThetaS=3;
AddiRecordInfo.SmaThetaB=0.35;
AddiRecordInfo.Nsma=20;
AddiRecordInfo.Smahc=3;
%
AddiRecordInfo.DoEast=0;
AddiRecordInfo.DoWest=1;
AddiRecordInfo.DoNorth=1;
AddiRecordInfo.DoSouth=1;

Day2sec=24*3600;
BeginTime=[2008 3 1 0 0 0];
EndTime=[2008 3 11 0 0 0];
BeginTimeDay=DATE_ConvertVect2mjd(BeginTime);
EndTimeDay=DATE_ConvertVect2mjd(EndTime);
BeginTimeSec=BeginTimeDay*Day2sec;
EndTimeSec=EndTimeDay*Day2sec;

UseSta=0;
UseHis=0;

% Note that the netcdf file output of ROMS
% are missing here. 
% For the example to be complte, we would have needed
% PrefixHis????.nc (history file of the roms model)
% StationFileName (station file output of the roms model)






if (UseSta == 1)
  %
  % creation of the station.in file
  TheRecord=OW_RecordBoundaryPoints(SmaGridFile);
  OW_OutputStations(...
      StationFile, ...
      TheRecord.LonCoord, TheRecord.LatCoord, ...
      TheRecord.ListComments);
  %
  % creation of the boundaries of the domain.
  TheBoundaryFile='BoundSta.nc';
  TheOffset=0;
  NEST_CreateBoundaryFromStation(...
      StationFileName, TheOffset, ...
      SmallGridFile, AddiRecordInfo, ...
      BeginTimeSec, EndTimeSec, TheBoundaryFile);
end;
if (UseHis == 1)
  TheBoundaryFile='BoundHis.nc';
  UseSpMat=1;
  NestingTotalArray=NEST_CreateNestingTotalArray(...
      BigGridFile, SmallGridFile, AddiRecordInfo, UseSpMat);
  NEST_CreateBoundaryFromHistory(...
      PrefixHis, NestingTotalArray, ...
      BeginTimeSec, EndTimeSec, TheBoundaryFile);
end;
