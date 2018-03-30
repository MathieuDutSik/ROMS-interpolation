function TheState=ReadSingleHistoryRecord(HistoryFile, GrdArr, iRecord)
%
% TheState=ReadSingleHistoryRecord(HistoryFile, GrdArr, iRecord)
% return the state of the ROMS history file HistoryFile
% computed with the ROMS grid array GrdArr
% at iRecord
%
nc=netcdf(HistoryFile, 'nowrite');
eTime=nc{'ocean_time'}(iRecord);
xy_rho=length(nc('xy_rho'));
if (xy_rho > 0)
  % we have a write water file
  xy_rho=length(nc('xy_rho'));
  xyz_rho=length(nc('xyz_rho'));
  s_rho=xyz_rho/xy_rho;
  ZETAwp=nc{'zeta'}(iRecord, :);
  TEMPwp=nc{'temp'}(iRecord, :);
  SALTwp=nc{'salt'}(iRecord, :);
  UBARwp=nc{'ubar'}(iRecord, :);
  VBARwp=nc{'vbar'}(iRecord, :);
  Uwp=nc{'u'}(iRecord, :);
  Vwp=nc{'v'}(iRecord, :);
  ZETA=Field2D_wp2usual(GrdArr.MSK_rho, ZETAwp);
  UBAR=Field2D_wp2usual(GrdArr.MSK_u, UBARwp);
  VBAR=Field2D_wp2usual(GrdArr.MSK_v, VBARwp);
  TEMP=Field3D_wp2usual(GrdArr.MSK_rho, s_rho, TEMPwp);
  SALT=Field3D_wp2usual(GrdArr.MSK_rho, s_rho, SALTwp);
  U=Field3D_wp2usual(GrdArr.MSK_u, s_rho, Uwp);
  V=Field3D_wp2usual(GrdArr.MSK_v, s_rho, Vwp);
else
  ZETA=nc{'zeta'}(iRecord, :, :);
  TEMP=nc{'temp'}(iRecord, :, :, :);
  SALT=nc{'salt'}(iRecord, :, :, :);
  UBAR=nc{'ubar'}(iRecord, :, :);
  VBAR=nc{'vbar'}(iRecord, :, :);
  U=nc{'u'}(iRecord, :, :, :);
  V=nc{'v'}(iRecord, :, :, :);
  s_rho=length(nc('s_rho'));
  ZETA(GrdArr.KlandRho)=0;
  UBAR(GrdArr.KlandU)=0;
  VBAR(GrdArr.KlandV)=0;
  TEMP=Field3D_PutToFillVal(GrdArr.MSK_rho, TEMP, 0);
  SALT=Field3D_PutToFillVal(GrdArr.MSK_rho, SALT, 0);
  U=Field3D_PutToFillVal(GrdArr.MSK_u, U, 0);
  V=Field3D_PutToFillVal(GrdArr.MSK_v, V, 0);
end;
close(nc);
TheState.ZETA=ZETA;
TheState.TEMP=TEMP;
TheState.SALT=SALT;
TheState.UBAR=UBAR;
TheState.VBAR=VBAR;
TheState.U=U;
TheState.V=V;
TheState.eTime=eTime;
