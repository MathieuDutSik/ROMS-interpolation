function NEST_CreateBoundaryFromStation(...
    StationFileName, TheOffset, ...
    SmallGridFile, AddiRecordInfo, ...
    BeginTimeSec, EndTimeSec, BoundaryFileName)
%
% --StationFileName is the output by ROMS of the run
% --TheOffset is the move in index that is needed to be
%   made (station file are for multiple purposes)
%   TheOffSet=0 if only for this purpose or it is at the beginning.
% --SmallGridFile is the small grid nested in the big one.
% --AddiRecordInfo is the record of additional info:
%   AddiRecordInfo.BigThetaS (big grid thetaS)
%   AddiRecordInfo.BigThetaB (big grid thetaB)
%   AddiRecordInfo.Bighc     (big grid hc)
%   AddiRecordInfo.Nbig      (big grid N)
%   AddiRecordInfo.SmaThetaS (small grid thetaS)
%   AddiRecordInfo.SmaThetaB (small grid thetaB)
%   AddiRecordInfo.Smahc     (small grid hc)
%   AddiRecordInfo.Nsma      (small grid N)
%   AddiRecordInfo.DoEast    (1 for doing it for small grid)
%   AddiRecordInfo.DoWest    (1 for doing it for small grid)
%   AddiRecordInfo.DoNorth   (1 for doing it for small grid)
%   AddiRecordInfo.DoSouth   (1 for doing it for small grid)
% --BeginTimeSec, and EndTimeSec are the sought
%   time interval in the boundary file
% --BoundaryFileName is the name of the netcdf boundary file.



TheRecord=NEST_RecordBoundaryPoints(SmallGridFile);
[BigSc_w, BigCs_w, BigSc_r, BigCs_r]=GRID_GetSc_Cs(...
    AddiRecordInfo.Nbig, AddiRecordInfo.BigThetaS, AddiRecordInfo.BigThetaB);
Bigcff_r=AddiRecordInfo.Bighc*(BigSc_r-BigCs_r);
[SmaSc_w, SmaCs_w, SmaSc_r, SmaCs_r]=GRID_GetSc_Cs(...
    AddiRecordInfo.Nsma, AddiRecordInfo.SmaThetaS, AddiRecordInfo.SmaThetaB);
Smacff_r=AddiRecordInfo.Smahc*(SmaSc_r-SmaCs_r);

ncSma=netcdf(SmallGridFile, 'nowrite');
eta_rho=length(ncSma('eta_rho'));
xi_rho=length(ncSma('xi_rho'));
eta_u=length(ncSma('eta_u'));
xi_u=length(ncSma('xi_u'));
eta_v=length(ncSma('eta_v'));
xi_v=length(ncSma('xi_v'));
ANGLE_rho=ncSma{'angle'}(:);
close(ncSma);
ANGLE_u=(ANGLE_rho(:, 1:xi_u)+ANGLE_rho(:, 2:xi_rho))/2;
ANGLE_v=(ANGLE_rho(1:eta_v, :)+ANGLE_rho(2:eta_rho, :))/2;

[MSKsta_tot, MSKsta_rho, MSKsta_u, MSKsta_v]= ...
    NEST_DetermineStationValidity(StationFileName);

ncSta=netcdf(StationFileName, 'nowrite');
nbTime=length(ncSta('ocean_time'));
LTime=ncSta{'ocean_time'}(:);
%
K=find(LTime == BeginTimeSec);
nbK=size(K,1);
if (nbK == 1)
  iTimeBegin=K(1,1);
else
  iTimeBegin=-1;
  for iTime=1:nbTime
    if (LTime(iTime,1) < BeginTimeSec)
      iTimeBegin=iTime;
    end;
  end;
  if (iTimeBegin == -1)
    disp(['Sorry but the begin time of station file is further']);
    disp(['from the asked begin time']);
    strStaBegin=DATE_ConvertMjd2mystringPres(LTime(1,1)/(24*3600));
    strBegin=DATE_ConvertMjd2mystringPres(BeginTimeSec/(24*3600));
    disp(['Starting time station=' strStaBegin]);
    disp(['  Asked starting time=' strBegin]);
    keyboard;
  end;
end;
disp(['iTimeBegin=' num2str(iTimeBegin)]);
%
K=find(LTime == EndTimeSec);
nbK=size(K,1);
if (nbK == 1)
  iTimeEnd=K(1,1);
else
  iTimeEnd=-1;
  for iTime=1:nbTime
    jTime=nbTime+1-iTime;
    if (LTime(jTime,1) > EndTimeSec)
      iTimeEnd=jTime;
    end;
  end;
  if (iTimeEnd == -1)
    disp(['Sorry but the end time of station file is short of']);
    disp(['the asked end time']);
    strStaEnd=DATE_ConvertMjd2mystringPres(LTime(nbTime,1)/(24*3600));
    strEnd=DATE_ConvertMjd2mystringPres(EndTimeSec/(24*3600));
    disp(['Ending time station=' strStaEnd]);
    disp(['  Asked ending time=' strEnd]);
    keyboard;
  end;
end;
disp(['iTimeEnd=' num2str(iTimeEnd)]);
%
nbTimeTot=iTimeEnd-(iTimeBegin-1);
disp(['nbTimeTot=', num2str(nbTimeTot)]);


angleSta=ncSta{'angle'}(:);
depSta=ncSta{'h'}(:);
LON_sta=ncSta{'lon_rho'}(:);
LAT_sta=ncSta{'lat_rho'}(:);
nbStationRel=size(TheRecord.lon_nest,1);
depStaSel=depSta(TheOffset+1:TheOffset+nbStationRel);
LON_staSel=LON_sta(TheOffset+1:TheOffset+nbStationRel);
LAT_staSel=LAT_sta(TheOffset+1:TheOffset+nbStationRel);
%
MSKsta_tot_rel=MSKsta_tot(TheOffset+1:TheOffset+nbStationRel, 1);
Kcorr=find(MSKsta_tot_rel == 1);
LON_staRel=LON_staSel(Kcorr);
LAT_staRel=LAT_staSel(Kcorr);

lon_nest=TheRecord.lon_nest;
lat_nest=TheRecord.lat_nest;
RelevantStation=zeros(nbStationRel,1);
for iStation=1:nbStationRel
  if (MSKsta_tot_rel(iStation,1) == 0)
    [index,distance]=nearxy(...
	LON_staRel,LAT_staRel,...
	lon_nest(iStation,1), ...
	lat_nest(iStation,1));
    RelevantStation(iStation,1)=Kcorr(index);
    disp(['We find station ' num2str(iStation) ' to be wrong']);
  else
    RelevantStation(iStation,1)=iStation;
  end;
  disp(['iStation=' num2str(iStation) '  RelSta=' ...
	num2str(RelevantStation(iStation,1))]);
end;

s_rho=AddiRecordInfo.Nsma;
[test, reason]=IsMakeableFile(BoundaryFileName);
if (test == 0)
  disp(['We cannot create file ' BoundaryFileName]);
  disp(['For reason : ' reason]);
  keyboard;
end;
ncBound=netcdf(BoundaryFileName, 'clobber');
ncBound('eta_rho')=eta_rho;
ncBound('xi_rho')=xi_rho;
ncBound('eta_u')=eta_u;
ncBound('xi_u')=xi_u;
ncBound('eta_v')=eta_v;
ncBound('xi_v')=xi_v;
ncBound('zeta_time')=nbTimeTot;
ncBound('v2d_time')=nbTimeTot;
ncBound('v3d_time')=nbTimeTot;
ncBound('temp_time')=nbTimeTot;
ncBound('salt_time')=nbTimeTot;
ncBound('s_rho')=s_rho;
ncBound('one')=1;
%
ncBound{'zeta_time'}=ncdouble('zeta_time');
ncBound{'v2d_time'}=ncdouble('v2d_time');
ncBound{'v3d_time'}=ncdouble('v3d_time');
ncBound{'temp_time'}=ncdouble('temp_time');
ncBound{'salt_time'}=ncdouble('salt_time');
ncBound{'has_debug'}=ncdouble('one');
ncBound{'has_debug'}(:)=0;
%
if (AddiRecordInfo.DoEast == 1)
  ncBound{'zeta_east'}=ncfloat('zeta_time', 'eta_rho');
  ncBound{'temp_east'}=ncfloat('temp_time', 's_rho', 'eta_rho');
  ncBound{'salt_east'}=ncfloat('salt_time', 's_rho', 'eta_rho');
  ncBound{'ubar_east'}=ncfloat('v2d_time', 'eta_u');
  ncBound{'vbar_east'}=ncfloat('v2d_time', 'eta_v');
  ncBound{'u_east'}=ncfloat('v3d_time', 's_rho', 'eta_u');
  ncBound{'v_east'}=ncfloat('v3d_time', 's_rho', 'eta_v');
end;
if (AddiRecordInfo.DoWest == 1)
  ncBound{'zeta_west'}=ncfloat('zeta_time', 'eta_rho');
  ncBound{'temp_west'}=ncfloat('temp_time', 's_rho', 'eta_rho');
  ncBound{'salt_west'}=ncfloat('salt_time', 's_rho', 'eta_rho');
  ncBound{'ubar_west'}=ncfloat('v2d_time', 'eta_u');
  ncBound{'vbar_west'}=ncfloat('v2d_time', 'eta_v');
  ncBound{'u_west'}=ncfloat('v3d_time', 's_rho', 'eta_u');
  ncBound{'v_west'}=ncfloat('v3d_time', 's_rho', 'eta_v');
end;
if (AddiRecordInfo.DoNorth == 1)
  ncBound{'zeta_north'}=ncfloat('zeta_time', 'xi_rho');
  ncBound{'temp_north'}=ncfloat('temp_time', 's_rho', 'xi_rho');
  ncBound{'salt_north'}=ncfloat('salt_time', 's_rho', 'xi_rho');
  ncBound{'ubar_north'}=ncfloat('v2d_time', 'xi_u');
  ncBound{'vbar_north'}=ncfloat('v2d_time', 'xi_v');
  ncBound{'u_north'}=ncfloat('v3d_time', 's_rho', 'xi_u');
  ncBound{'v_north'}=ncfloat('v3d_time', 's_rho', 'xi_v');
end;
if (AddiRecordInfo.DoSouth == 1)
  ncBound{'zeta_south'}=ncfloat('zeta_time', 'xi_rho');
  ncBound{'temp_south'}=ncfloat('temp_time', 's_rho', 'xi_rho');
  ncBound{'salt_south'}=ncfloat('salt_time', 's_rho', 'xi_rho');
  ncBound{'ubar_south'}=ncfloat('v2d_time', 'xi_u');
  ncBound{'vbar_south'}=ncfloat('v2d_time', 'xi_v');
  ncBound{'u_south'}=ncfloat('v3d_time', 's_rho', 'xi_u');
  ncBound{'v_south'}=ncfloat('v3d_time', 's_rho', 'xi_v');
end;

Day2sec=24*3600;
for iRecordTime=iTimeBegin:iTimeEnd
  iTime=iRecordTime-(iTimeBegin-1);
  disp(['iTime=' num2str(iTime) '  iRecordTime=' num2str(iRecordTime)]);
  eTime=ncSta{'ocean_time'}(iRecordTime);
  ncBound{'zeta_time'}(iTime)=eTime/Day2sec;
  ncBound{'v2d_time'}(iTime)=eTime/Day2sec;
  ncBound{'v3d_time'}(iTime)=eTime/Day2sec;
  ncBound{'temp_time'}(iTime)=eTime/Day2sec;
  ncBound{'salt_time'}(iTime)=eTime/Day2sec;
  TempField=ncSta{'temp'}(...
      iRecordTime,TheOffset+1:TheOffset+nbStationRel,:);
  SaltField=ncSta{'salt'}(...
      iRecordTime,TheOffset+1:TheOffset+nbStationRel,:);
  ZetaField=ncSta{'zeta'}(iRecordTime, TheOffset+1:TheOffset+nbStationRel);
  UbarField=ncSta{'ubar'}(iRecordTime, TheOffset+1:TheOffset+nbStationRel);
  VbarField=ncSta{'vbar'}(iRecordTime, TheOffset+1:TheOffset+nbStationRel);
  UField=ncSta{'u'}(iRecordTime,TheOffset+1:TheOffset+nbStationRel,:);
  VField=ncSta{'v'}(iRecordTime,TheOffset+1:TheOffset+nbStationRel,:);
  if (AddiRecordInfo.DoEast == 1)
    iBeginRho=TheRecord.rho_east_begin;
    iBeginU=TheRecord.u_east_begin;
    iBeginV=TheRecord.v_east_begin;
    ZETAeast=zeros(eta_rho, 1);
    TEMPeast=zeros(s_rho, eta_rho);
    SALTeast=zeros(s_rho, eta_rho);
    UBAReast=zeros(eta_u, 1);
    VBAReast=zeros(eta_v, 1);
    Ueast=zeros(s_rho, eta_u);
    Veast=zeros(s_rho, eta_v);
    for iEta=1:eta_rho
      if (TheRecord.rho_east_info(iEta, 1) == 1)
	iSta=RelevantStation(iBeginRho,1);
	depStaFind=depSta(iSta);
	ZETAeast(iEta, 1)=ZetaField(1, iSta);
	BigDepth=Bigcff_r+BigCs_r*depSta(iSta);
	SmaDepth=Smacff_r+SmaCs_r*TheRecord.dep_nest(iSta);
	BigTempField=TempField(iSta, :);
	BigSaltField=SaltField(iSta, :);
	SmaTempField=interp_vertical(BigDepth, BigTempField, SmaDepth);
	SmaSaltField=interp_vertical(BigDepth, BigSaltField, SmaDepth);
	TEMPeast(:, iEta)=SmaTempField(:);
	SALTeast(:, iEta)=SmaSaltField(:);
	iBeginRho=iBeginRho+1;
      end;
    end;
    for iEta=1:eta_u
      if (TheRecord.u_east_info(iEta, 1) == 1)
	iSta=RelevantStation(iBeginU,1);
	depStaFind=depSta(iSta);
	BigAngle=angleSta(iSta,1);
	SmaAngle=TheRecord.ang_nest(iSta);
	deltaAng=BigAngle-SmaAngle;
	%
	AngBigUBAReast=UbarField(1, iSta);
	AngBigVBAReast=VbarField(1, iSta);
	AngSmaUBAReast=cos(deltaAng)*AngBigUBAReast-...
	    sin(deltaAng)*AngBigVBAReast;
	AngSmaVBAReast=sin(deltaAng)*AngBigUBAReast+...
	    cos(deltaAng)*AngBigVBAReast;
	UBAReast(iEta, 1)=AngSmaUBAReast;
	%
	AngBigUeast=UField(iSta, :);
	AngBigVeast=VField(iSta, :);
	AngSmaUeast=cos(deltaAng)*AngBigUeast-...
	    sin(deltaAng)*AngBigVeast;
	AngSmaVeast=sin(deltaAng)*AngBigUeast+...
	    cos(deltaAng)*AngBigVeast;
	Ueast(:, iEta)=AngSmaUeast;
	iBeginU=iBeginU+1;
      end;
    end;
    for iEta=1:eta_v
      if (TheRecord.v_east_info(iEta, 1) == 1)
	iSta=RelevantStation(iBeginV,1);
	depStaFind=depSta(iSta);
	BigAngle=angleSta(iSta,1);
	SmaAngle=TheRecord.ang_nest(iSta);
	deltaAng=BigAngle-SmaAngle;
	%
	AngBigUBAReast=UbarField(1, iSta);
	AngBigVBAReast=VbarField(1, iSta);
	AngSmaUBAReast=cos(deltaAng)*AngBigUBAReast-...
	    sin(deltaAng)*AngBigVBAReast;
	AngSmaVBAReast=sin(deltaAng)*AngBigUBAReast+...
	    cos(deltaAng)*AngBigVBAReast;
	VBAReast(iEta, 1)=AngSmaVBAReast;
	%
	AngBigUeast=UField(iSta, :);
	AngBigVeast=VField(iSta, :);
	AngSmaUeast=cos(deltaAng)*AngBigUeast-...
	    sin(deltaAng)*AngBigVeast;
	AngSmaVeast=sin(deltaAng)*AngBigUeast+...
	    cos(deltaAng)*AngBigVeast;
	Veast(:, iEta)=AngSmaVeast;
	iBeginV=iBeginV+1;
      end;
    end;
    ncBound{'zeta_east'}(iTime, :)=ZETAeast;
    ncBound{'temp_east'}(iTime, :, :)=TEMPeast;
    ncBound{'salt_east'}(iTime, :, :)=SALTeast;
    ncBound{'ubar_east'}(iTime, :)=UBAReast;
    ncBound{'vbar_east'}(iTime, :)=VBAReast;
    ncBound{'u_east'}(iTime, :, :)=Ueast;
    ncBound{'v_east'}(iTime, :, :)=Veast;
  end;
  %
  if (AddiRecordInfo.DoWest == 1)
    iBeginRho=TheRecord.rho_west_begin;
    iBeginU=TheRecord.u_west_begin;
    iBeginV=TheRecord.v_west_begin;
    ZETAwest=zeros(eta_rho, 1);
    TEMPwest=zeros(s_rho, eta_rho);
    SALTwest=zeros(s_rho, eta_rho);
    UBARwest=zeros(eta_u, 1);
    VBARwest=zeros(eta_v, 1);
    Uwest=zeros(s_rho, eta_u);
    Vwest=zeros(s_rho, eta_v);
    for iEta=1:eta_rho
      if (TheRecord.rho_west_info(iEta, 1) == 1)
	iSta=RelevantStation(iBeginRho,1);
	ZETAwest(iEta, 1)=ZetaField(1, iSta);
	BigDepth=Bigcff_r+BigCs_r*depSta(iSta);
	SmaDepth=Smacff_r+SmaCs_r*TheRecord.dep_nest(iSta);
	BigTempField=TempField(iSta, :);
	BigSaltField=SaltField(iSta, :);
	SmaTempField=interp_vertical(BigDepth, BigTempField, SmaDepth);
	SmaSaltField=interp_vertical(BigDepth, BigSaltField, SmaDepth);
	TEMPwest(:, iEta)=SmaTempField(:);
	SALTwest(:, iEta)=SmaSaltField(:);
	iBeginRho=iBeginRho+1;
      end;
    end;
    for iEta=1:eta_u
      if (TheRecord.u_west_info(iEta, 1) == 1)
	iSta=RelevantStation(iBeginU,1);
	BigAngle=angleSta(iSta,1);
	SmaAngle=TheRecord.ang_nest(iSta);
	deltaAng=BigAngle-SmaAngle;
	%
	AngBigUBARwest=UbarField(1, iSta);
	AngBigVBARwest=VbarField(1, iSta);
	AngSmaUBARwest=cos(deltaAng)*AngBigUBARwest-...
	    sin(deltaAng)*AngBigVBARwest;
	AngSmaVBARwest=sin(deltaAng)*AngBigUBARwest+...
	    cos(deltaAng)*AngBigVBARwest;
	UBARwest(iEta, 1)=AngSmaUBARwest;
	%
	AngBigUwest=UField(iSta, :);
	AngBigVwest=VField(iSta, :);
	AngSmaUwest=cos(deltaAng)*AngBigUwest-...
	    sin(deltaAng)*AngBigVwest;
	AngSmaVwest=sin(deltaAng)*AngBigUwest+...
	    cos(deltaAng)*AngBigVwest;
	Uwest(:, iEta)=AngSmaUwest;
	iBeginU=iBeginU+1;
      end;
    end;
    for iEta=1:eta_v
      if (TheRecord.v_west_info(iEta, 1) == 1)
	iSta=RelevantStation(iBeginV,1);
	BigAngle=angleSta(iSta,1);
	SmaAngle=TheRecord.ang_nest(iSta);
	deltaAng=BigAngle-SmaAngle;
	%
	AngBigUBARwest=UbarField(1, iSta);
	AngBigVBARwest=VbarField(1, iSta);
	AngSmaUBARwest=cos(deltaAng)*AngBigUBARwest-...
	    sin(deltaAng)*AngBigVBARwest;
	AngSmaVBARwest=sin(deltaAng)*AngBigUBARwest+...
	    cos(deltaAng)*AngBigVBARwest;
	VBARwest(iEta, 1)=AngSmaVBARwest;
	%
	AngBigUwest=UField(iSta, :);
	AngBigVwest=VField(iSta, :);
	AngSmaUwest=cos(deltaAng)*AngBigUwest-...
	    sin(deltaAng)*AngBigVwest;
	AngSmaVwest=sin(deltaAng)*AngBigUwest+...
	    cos(deltaAng)*AngBigVwest;
	Vwest(:, iEta)=AngSmaVwest;
	iBeginV=iBeginV+1;
      end;
    end;
    ncBound{'zeta_west'}(iTime, :)=ZETAwest;
    ncBound{'temp_west'}(iTime, :, :)=TEMPwest;
    ncBound{'salt_west'}(iTime, :, :)=SALTwest;
    ncBound{'ubar_west'}(iTime, :)=UBARwest;
    ncBound{'vbar_west'}(iTime, :)=VBARwest;
    ncBound{'u_west'}(iTime, :, :)=Uwest;
    ncBound{'v_west'}(iTime, :, :)=Vwest;
  end;
  %
  if (AddiRecordInfo.DoSouth == 1)
    iBeginRho=TheRecord.rho_south_begin;
    iBeginU=TheRecord.u_south_begin;
    iBeginV=TheRecord.v_south_begin;
    ZETAsouth=zeros(xi_rho, 1);
    TEMPsouth=zeros(s_rho, xi_rho);
    SALTsouth=zeros(s_rho, xi_rho);
    UBARsouth=zeros(xi_u, 1);
    VBARsouth=zeros(xi_v, 1);
    Usouth=zeros(s_rho, xi_u);
    Vsouth=zeros(s_rho, xi_v);
    for iXi=1:xi_rho
      if (TheRecord.rho_south_info(iXi, 1) == 1)
	iSta=RelevantStation(iBeginRho,1);
	ZETAsouth(iXi, 1)=ZetaField(1, iSta);
	BigDepth=Bigcff_r+BigCs_r*depSta(iSta);
	SmaDepth=Smacff_r+SmaCs_r*TheRecord.dep_nest(iSta);
	BigTempField=TempField(iSta, :);
	BigSaltField=SaltField(iSta, :);
	SmaTempField=interp_vertical(BigDepth, BigTempField, SmaDepth);
	SmaSaltField=interp_vertical(BigDepth, BigSaltField, SmaDepth);
	TEMPsouth(:, iXi)=SmaTempField(:);
	SALTsouth(:, iXi)=SmaSaltField(:);
	iBeginRho=iBeginRho+1;
      end;
    end;
    for iXi=1:xi_u
      if (TheRecord.u_south_info(iXi, 1) == 1)
	iSta=RelevantStation(iBeginU,1);
	BigAngle=angleSta(iSta,1);
	SmaAngle=TheRecord.ang_nest(iSta);
	deltaAng=BigAngle-SmaAngle;
	%
	AngBigUBARsouth=UbarField(1, iSta);
	AngBigVBARsouth=VbarField(1, iSta);
	AngSmaUBARsouth=cos(deltaAng)*AngBigUBARsouth-...
	    sin(deltaAng)*AngBigVBARsouth;
	AngSmaVBARsouth=sin(deltaAng)*AngBigUBARsouth+...
	    cos(deltaAng)*AngBigVBARsouth;
	UBARsouth(iXi, 1)=AngSmaUBARsouth;
	%
	AngBigUsouth=UField(iSta, :);
	AngBigVsouth=VField(iSta, :);
	AngSmaUsouth=cos(deltaAng)*AngBigUsouth-...
	    sin(deltaAng)*AngBigVsouth;
	AngSmaVsouth=sin(deltaAng)*AngBigUsouth+...
	    cos(deltaAng)*AngBigVsouth;
	Usouth(:, iXi)=AngSmaUsouth;
	iBeginU=iBeginU+1;
      end;
    end;
    for iXi=1:xi_v
      if (TheRecord.v_south_info(iXi, 1) == 1)
	iSta=RelevantStation(iBeginV,1);
	BigAngle=angleSta(iSta,1);
	SmaAngle=TheRecord.ang_nest(iSta);
	deltaAng=BigAngle-SmaAngle;
	%
	AngBigUBARsouth=UbarField(1, iSta);
	AngBigVBARsouth=VbarField(1, iSta);
	AngSmaUBARsouth=cos(deltaAng)*AngBigUBARsouth-...
	    sin(deltaAng)*AngBigVBARsouth;
	AngSmaVBARsouth=sin(deltaAng)*AngBigUBARsouth+...
	    cos(deltaAng)*AngBigVBARsouth;
	VBARsouth(iXi, 1)=AngSmaVBARsouth;
	%
	AngBigUsouth=UField(iSta, :);
	AngBigVsouth=VField(iSta, :);
	AngSmaUsouth=cos(deltaAng)*AngBigUsouth-...
	    sin(deltaAng)*AngBigVsouth;
	AngSmaVsouth=sin(deltaAng)*AngBigUsouth+...
	    cos(deltaAng)*AngBigVsouth;
	Vsouth(:, iXi)=AngSmaVsouth;
	iBeginV=iBeginV+1;
      end;
    end;
    ncBound{'zeta_south'}(iTime, :)=ZETAsouth;
    ncBound{'temp_south'}(iTime, :, :)=TEMPsouth;
    ncBound{'salt_south'}(iTime, :, :)=SALTsouth;
    ncBound{'ubar_south'}(iTime, :)=UBARsouth;
    ncBound{'vbar_south'}(iTime, :)=VBARsouth;
    ncBound{'u_south'}(iTime, :, :)=Usouth;
    ncBound{'v_south'}(iTime, :, :)=Vsouth;
  end;
  %
  if (AddiRecordInfo.DoNorth == 1)
    iBeginRho=TheRecord.rho_north_begin;
    iBeginU=TheRecord.u_north_begin;
    iBeginV=TheRecord.v_north_begin;
    ZETAnorth=zeros(xi_rho, 1);
    TEMPnorth=zeros(s_rho, xi_rho);
    SALTnorth=zeros(s_rho, xi_rho);
    UBARnorth=zeros(xi_u, 1);
    VBARnorth=zeros(xi_v, 1);
    Unorth=zeros(s_rho, xi_u);
    Vnorth=zeros(s_rho, xi_v);
    for iXi=1:xi_rho
      if (TheRecord.rho_north_info(iXi, 1) == 1)
	iSta=RelevantStation(iBeginRho,1);
	ZETAnorth(iXi, 1)=ZetaField(1, iSta);
	BigDepth=Bigcff_r+BigCs_r*depSta(iSta);
	SmaDepth=Smacff_r+SmaCs_r*TheRecord.dep_nest(iSta);
	BigTempField=TempField(iSta, :);
	BigSaltField=SaltField(iSta, :);
	SmaTempField=interp_vertical(BigDepth, BigTempField, SmaDepth);
	SmaSaltField=interp_vertical(BigDepth, BigSaltField, SmaDepth);
	TEMPnorth(:, iXi)=SmaTempField(:);
	SALTnorth(:, iXi)=SmaSaltField(:);
	iBeginRho=iBeginRho+1;
      end;
    end;
    for iXi=1:xi_u
      if (TheRecord.u_north_info(iXi, 1) == 1)
	iSta=RelevantStation(iBeginU,1);
	BigAngle=angleSta(iSta,1);
	SmaAngle=TheRecord.ang_nest(iSta);
	deltaAng=BigAngle-SmaAngle;
	%
	AngBigUBARnorth=UbarField(1, iSta);
	AngBigVBARnorth=VbarField(1, iSta);
	AngSmaUBARnorth=cos(deltaAng)*AngBigUBARnorth-...
	    sin(deltaAng)*AngBigVBARnorth;
	AngSmaVBARnorth=sin(deltaAng)*AngBigUBARnorth+...
	    cos(deltaAng)*AngBigVBARnorth;
	UBARnorth(iXi, 1)=AngSmaUBARnorth;
	%
	AngBigUnorth=UField(iSta, :);
	AngBigVnorth=VField(iSta, :);
	AngSmaUnorth=cos(deltaAng)*AngBigUnorth-...
	    sin(deltaAng)*AngBigVnorth;
	AngSmaVnorth=sin(deltaAng)*AngBigUnorth+...
	    cos(deltaAng)*AngBigVnorth;
	Unorth(:, iXi)=AngSmaUnorth;
	iBeginU=iBeginU+1;
      end;
    end;
    for iXi=1:xi_v
      if (TheRecord.v_north_info(iXi, 1) == 1)
	iSta=RelevantStation(iBeginV,1);
	BigAngle=angleSta(iSta,1);
	SmaAngle=TheRecord.ang_nest(iSta);
	deltaAng=BigAngle-SmaAngle;
	%
	AngBigUBARnorth=UbarField(1, iSta);
	AngBigVBARnorth=VbarField(1, iSta);
	AngSmaUBARnorth=cos(deltaAng)*AngBigUBARnorth-...
	    sin(deltaAng)*AngBigVBARnorth;
	AngSmaVBARnorth=sin(deltaAng)*AngBigUBARnorth+...
	    cos(deltaAng)*AngBigVBARnorth;
	VBARnorth(iXi, 1)=AngSmaVBARnorth;
	%
	AngBigUnorth=UField(iSta, :);
	AngBigVnorth=VField(iSta, :);
	AngSmaUnorth=cos(deltaAng)*AngBigUnorth-...
	    sin(deltaAng)*AngBigVnorth;
	AngSmaVnorth=sin(deltaAng)*AngBigUnorth+...
	    cos(deltaAng)*AngBigVnorth;
	Vnorth(:, iXi)=AngSmaVnorth;
	iBeginV=iBeginV+1;
      end;
    end;
    ncBound{'zeta_north'}(iTime, :)=ZETAnorth;
    ncBound{'temp_north'}(iTime, :, :)=TEMPnorth;
    ncBound{'salt_north'}(iTime, :, :)=SALTnorth;
    ncBound{'ubar_north'}(iTime, :)=UBARnorth;
    ncBound{'vbar_north'}(iTime, :)=VBARnorth;
    ncBound{'u_north'}(iTime, :, :)=Unorth;
    ncBound{'v_north'}(iTime, :, :)=Vnorth;
  end;
end;
close(ncBound);
disp(['NEST_CreateBoundaryFromStation finished']);