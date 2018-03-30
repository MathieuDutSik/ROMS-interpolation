function [MSKsta_tot, MSKsta_rho, MSKsta_u, MSKsta_v]=...
    NEST_DetermineStationValidity(StationOutFile)
nc=netcdf(StationOutFile, 'nowrite');
DEP=nc{'h'}(:);
UBAR=nc{'ubar'}(1, :);
VBAR=nc{'vbar'}(1, :);
nbStation=length(nc('station'));
close(nc);
MSKsta_rho=zeros(nbStation, 1);
MSKsta_u=zeros(nbStation, 1);
MSKsta_v=zeros(nbStation, 1);
for iStation=1:nbStation
  if (DEP(iStation, 1) > 10000)
    MSKsta_rho(iStation,1)=0;
  else
    MSKsta_rho(iStation,1)=1;
  end;
  if (UBAR(1, iStation) > 10000)
    MSKsta_u(iStation,1)=0;
  else
    MSKsta_u(iStation,1)=1;
  end;
  if (VBAR(1, iStation) > 10000)
    MSKsta_v(iStation,1)=0;
  else
    MSKsta_v(iStation,1)=1;
  end;
end;
MSKsta_tot=MSKsta_rho.*MSKsta_u.*MSKsta_v;
