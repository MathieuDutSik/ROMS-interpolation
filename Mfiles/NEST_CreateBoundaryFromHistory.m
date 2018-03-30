function NEST_CreateBoundaryFromHistory(...
    PrefixHis, NestingTotalArray, ...
    BeginTimeSec, EndTimeSec, BoundaryFileName)
%
% --StationFileName is the output by ROMS of the run
% --TheOffset is the move in index that is needed to be
%   made (station file are for multiple purposes)
%   TheOffSet=0 if only for this purpose or it is at the beginning.

eta_rho_sma=NestingTotalArray.eta_rho_sma;
xi_rho_sma=NestingTotalArray.xi_rho_sma;
eta_u_sma=NestingTotalArray.eta_u_sma;
xi_u_sma=NestingTotalArray.xi_u_sma;
eta_v_sma=NestingTotalArray.eta_v_sma;
xi_v_sma=NestingTotalArray.xi_v_sma;
%
[ListTimeHistory, ListIFile, ListIRecord, nbFile, ListAbsoluteIndex]=...
    ROMShistoryGetInfo2(PrefixHis, BeginTimeSec, EndTimeSec);
nbTime=size(ListTimeHistory,1);
%
s_rho=NestingTotalArray.AddiRecordInfo.Nsma;
[test, reason]=IsMakeableFile(BoundaryFileName);
if (test == 0)
  disp(['We cannot create file ' BoundaryFileName]);
  disp(['For reason : ' reason]);
  keyboard;
end;
ncBound=netcdf(BoundaryFileName, 'clobber');
ncBound('eta_rho')=eta_rho_sma;
ncBound('xi_rho')=xi_rho_sma;
ncBound('eta_u')=eta_u_sma;
ncBound('xi_u')=xi_u_sma;
ncBound('eta_v')=eta_v_sma;
ncBound('xi_v')=xi_v_sma;
ncBound('zeta_time')=nbTime;
ncBound('v2d_time')=nbTime;
ncBound('v3d_time')=nbTime;
ncBound('temp_time')=nbTime;
ncBound('salt_time')=nbTime;
ncBound('s_rho')=s_rho;
%
ncBound{'zeta_time'}=ncdouble('zeta_time');
ncBound{'v2d_time'}=ncdouble('v2d_time');
ncBound{'v3d_time'}=ncdouble('v3d_time');
ncBound{'temp_time'}=ncdouble('temp_time');
ncBound{'salt_time'}=ncdouble('salt_time');
%
if (NestingTotalArray.AddiRecordInfo.DoEast == 1)
  ncBound{'zeta_east'}=ncfloat('zeta_time', 'eta_rho');
  ncBound{'temp_east'}=ncfloat('temp_time', 's_rho', 'eta_rho');
  ncBound{'salt_east'}=ncfloat('salt_time', 's_rho', 'eta_rho');
  ncBound{'ubar_east'}=ncfloat('v2d_time', 'eta_u');
  ncBound{'vbar_east'}=ncfloat('v2d_time', 'eta_v');
  ncBound{'u_east'}=ncfloat('v3d_time', 's_rho', 'eta_u');
  ncBound{'v_east'}=ncfloat('v3d_time', 's_rho', 'eta_v');
end;
if (NestingTotalArray.AddiRecordInfo.DoWest == 1)
  ncBound{'zeta_west'}=ncfloat('zeta_time', 'eta_rho');
  ncBound{'temp_west'}=ncfloat('temp_time', 's_rho', 'eta_rho');
  ncBound{'salt_west'}=ncfloat('salt_time', 's_rho', 'eta_rho');
  ncBound{'ubar_west'}=ncfloat('v2d_time', 'eta_u');
  ncBound{'vbar_west'}=ncfloat('v2d_time', 'eta_v');
  ncBound{'u_west'}=ncfloat('v3d_time', 's_rho', 'eta_u');
  ncBound{'v_west'}=ncfloat('v3d_time', 's_rho', 'eta_v');
end;
if (NestingTotalArray.AddiRecordInfo.DoNorth == 1)
  ncBound{'zeta_north'}=ncfloat('zeta_time', 'xi_rho');
  ncBound{'temp_north'}=ncfloat('temp_time', 's_rho', 'xi_rho');
  ncBound{'salt_north'}=ncfloat('salt_time', 's_rho', 'xi_rho');
  ncBound{'ubar_north'}=ncfloat('v2d_time', 'xi_u');
  ncBound{'vbar_north'}=ncfloat('v2d_time', 'xi_v');
  ncBound{'u_north'}=ncfloat('v3d_time', 's_rho', 'xi_u');
  ncBound{'v_north'}=ncfloat('v3d_time', 's_rho', 'xi_v');
end;
if (NestingTotalArray.AddiRecordInfo.DoSouth == 1)
  ncBound{'zeta_south'}=ncfloat('zeta_time', 'xi_rho');
  ncBound{'temp_south'}=ncfloat('temp_time', 's_rho', 'xi_rho');
  ncBound{'salt_south'}=ncfloat('salt_time', 's_rho', 'xi_rho');
  ncBound{'ubar_south'}=ncfloat('v2d_time', 'xi_u');
  ncBound{'vbar_south'}=ncfloat('v2d_time', 'xi_v');
  ncBound{'u_south'}=ncfloat('v3d_time', 's_rho', 'xi_u');
  ncBound{'v_south'}=ncfloat('v3d_time', 's_rho', 'xi_v');
end;

Day2sec=24*3600;
for iTime=1:nbTime
  iFile=ListIFile(iTime,1);
  iRecord=ListIRecord(iTime,1);
  HistoryFile=[PrefixHis StringNumber(iFile, 4) '.nc'];
  TheState=ReadSingleHistoryRecord(...
      HistoryFile, NestingTotalArray.GrdArrBig, iRecord);
  disp(['iTime=' num2str(iTime)]);
  eTime=TheState.eTime;
  ncBound{'zeta_time'}(iTime)=eTime/Day2sec;
  ncBound{'v2d_time'}(iTime)=eTime/Day2sec;
  ncBound{'v3d_time'}(iTime)=eTime/Day2sec;
  ncBound{'temp_time'}(iTime)=eTime/Day2sec;
  ncBound{'salt_time'}(iTime)=eTime/Day2sec;
  if (NestingTotalArray.UseSpMat == 1)
    ZETAsma=InterpolSpMat_R2R_2Dfield(...
	NestingTotalArray.ArrayBigSma2D, ...
	TheState.ZETA);
    TEMPsma=InterpolSpMat_R2R_3Dfield(...
	NestingTotalArray.ArrayBigSma3D, ...
	TheState.TEMP);
    SALTsma=InterpolSpMat_R2R_3Dfield(...
	NestingTotalArray.ArrayBigSma3D, ...
	TheState.SALT);
    [UBARsma, VBARsma]=InterpolSpMat_R2R_2Duvfield(...
	NestingTotalArray.ArrayBigSma2Duv, ...
	TheState.UBAR, TheState.VBAR);
    [Usma, Vsma]=InterpolSpMat_R2R_3Duvfield(...
	NestingTotalArray.ArrayBigSma3Duv, ...
	TheState.U, TheState.V);
  else
    ZETAsma=InterpolMemEff_R2R_2Dfield(...
	NestingTotalArray, TheState.ZETA);
    TEMPsma=InterpolMemEff_R2R_3Dfield(...
	NestingTotalArray, NestingTotalArray.AddiRecordInfo, ...
	TheState.TEMP);
    SALTsma=InterpolMemEff_R2R_3Dfield(...
	NestingTotalArray, NestingTotalArray.AddiRecordInfo, ...
	TheState.SALT);
    [UBARsma, VBARsma]=InterpolMemEff_R2R_2Duvfield(...
	NestingTotalArray, ...
	TheState.UBAR, TheState.VBAR);
    [Usma, Vsma]=InterpolMemEff_R2R_3Duvfield(...
	NestingTotalArray, NestingTotalArray.AddiRecordInfo, ...
	TheState.U, TheState.V);
  end;
  if (NestingTotalArray.AddiRecordInfo.DoEast == 1)
    ZETAeast=squeeze(ZETAsma(:, xi_rho_sma));
    TEMPeast=squeeze(TEMPsma(:, :, xi_rho_sma));
    SALTeast=squeeze(SALTsma(:, :, xi_rho_sma));
    UBAReast=squeeze(UBARsma(:, xi_u_sma));
    VBAReast=squeeze(VBARsma(:, xi_v_sma));
    Ueast=squeeze(Usma(:, :, xi_u_sma));
    Veast=squeeze(Vsma(:, :, xi_v_sma));
    ncBound{'zeta_east'}(iTime, :)=ZETAeast;
    ncBound{'temp_east'}(iTime, :, :)=TEMPeast;
    ncBound{'salt_east'}(iTime, :, :)=SALTeast;
    ncBound{'ubar_east'}(iTime, :)=UBAReast;
    ncBound{'vbar_east'}(iTime, :)=VBAReast;
    ncBound{'u_east'}(iTime, :, :)=Ueast;
    ncBound{'v_east'}(iTime, :, :)=Veast;
  end;
  %
  if (NestingTotalArray.AddiRecordInfo.DoWest == 1)
    ZETAwest=squeeze(ZETAsma(:, 1));
    TEMPwest=squeeze(TEMPsma(:, :, 1));
    SALTwest=squeeze(SALTsma(:, :, 1));
    UBARwest=squeeze(UBARsma(:, 1));
    VBARwest=squeeze(VBARsma(:, 1));
    Uwest=squeeze(Usma(:, :, 1));
    Vwest=squeeze(Vsma(:, :, 1));
    ncBound{'zeta_west'}(iTime, :)=ZETAwest;
    ncBound{'temp_west'}(iTime, :, :)=TEMPwest;
    ncBound{'salt_west'}(iTime, :, :)=SALTwest;
    ncBound{'ubar_west'}(iTime, :)=UBARwest;
    ncBound{'vbar_west'}(iTime, :)=VBARwest;
    ncBound{'u_west'}(iTime, :, :)=Uwest;
    ncBound{'v_west'}(iTime, :, :)=Vwest;
  end;
  %
  if (NestingTotalArray.AddiRecordInfo.DoSouth == 1)
    ZETAsouth=squeeze(ZETAsma(1, :));
    TEMPsouth=squeeze(TEMPsma(:, 1, :));
    SALTsouth=squeeze(SALTsma(:, 1, :));
    UBARsouth=squeeze(UBARsma(1, :));
    VBARsouth=squeeze(VBARsma(1, :));
    Usouth=squeeze(Usma(:, 1, :));
    Vsouth=squeeze(Vsma(:, 1, :));
    ncBound{'zeta_south'}(iTime, :)=ZETAsouth;
    ncBound{'temp_south'}(iTime, :, :)=TEMPsouth;
    ncBound{'salt_south'}(iTime, :, :)=SALTsouth;
    ncBound{'ubar_south'}(iTime, :)=UBARsouth;
    ncBound{'vbar_south'}(iTime, :)=VBARsouth;
    ncBound{'u_south'}(iTime, :, :)=Usouth;
    ncBound{'v_south'}(iTime, :, :)=Vsouth;
  end;
  %
  if (NestingTotalArray.AddiRecordInfo.DoNorth == 1)
    ZETAnorth=squeeze(ZETAsma(eta_rho_sma,:));
    TEMPnorth=squeeze(TEMPsma(:, eta_rho_sma,:));
    SALTnorth=squeeze(SALTsma(:, eta_rho_sma,:));
    UBARnorth=squeeze(UBARsma(eta_u_sma,:));
    VBARnorth=squeeze(VBARsma(eta_v_sma,:));
    Unorth=squeeze(Usma(:, eta_u_sma,:));
    Vnorth=squeeze(Vsma(:, eta_v_sma,:));
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
