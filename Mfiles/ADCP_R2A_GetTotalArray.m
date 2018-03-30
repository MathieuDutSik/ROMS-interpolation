function TotalArray=ADCP_R2A_GetTotalArray(...
    GridFile, AddiRecordVertDisc, ...
    ListLonADCP, ListLatADCP, ADCPverticalLevels)

%
% ---GridFile is the ROMS grid file
% ---AddiRecordVertDisc is the description of the vertical
%    discretization
%    AddiRecordVertDisc.ThetaS
%    AddiRecordVertDisc.ThetaB
%    AddiRecordVertDisc.N
%    AddiRecordVertDisc.hc
% ---ListLonADCP, ListLatADCP are arrays of ADCP in the format
%    real(nbADCP, 1)
% ---ADCPverticalLevels is a matrix of depth of the form
%    real(nbADCP, nbVert)
%    if no depth is available (numbr of levels can vary from
%    one ADCP to another, then please put NaN)


nc=netcdf(GridFile, 'nowrite');
LON_rho=nc{'lon_rho'}(:);
LAT_rho=nc{'lat_rho'}(:);
LON_u=nc{'lon_u'}(:);
LAT_u=nc{'lat_u'}(:);
LON_v=nc{'lon_v'}(:);
LAT_v=nc{'lat_v'}(:);
MSK_rho=nc{'mask_rho'}(:);
MSK_u=nc{'mask_u'}(:);
MSK_v=nc{'mask_v'}(:);
DEP_rho=nc{'h'}(:);
ANG_rho=nc{'angle'}(:);
close(nc);

nbADCP=size(ListLonADCP, 1);
ListMsk=ones(nbADCP, 1);

[eta_rho_big, xi_rho_big]=size(DEP_rho);
DEPbig_u=(DEP_rho(:, 1:xi_rho_big-1)+DEP_rho(:, 2:xi_rho_big))/2;
DEPbig_v=(DEP_rho(1:eta_rho_big-1, :)+DEP_rho(2:eta_rho_big, :))/2;

[ListRelETA_u, ListRelXI_u, ListRelCoeff_u]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
	LON_u, LAT_u, MSK_u, ...
	ListLonADCP, ListLatADCP, ListMsk);
[ListRelETA_v, ListRelXI_v, ListRelCoeff_v]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
	LON_v, LAT_v, MSK_v, ...
	ListLonADCP, ListLatADCP, ListMsk);
[ListRelETA_rho, ListRelXI_rho, ListRelCoeff_rho]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
	LON_rho, LAT_rho, MSK_rho, ...
	ListLonADCP, ListLatADCP, ListMsk);
ListAng=zeros(nbADCP, 1);
for iADCP=1:nbADCP
  eSum=0;
  eSumWeight=0;
  for idx=1:4
    eWeight=ListRelCoeff_rho(iADCP,1,idx);
    if (eWeight > 0)
      iEtaBig=ListRelETA_rho(iADCP, 1, idx);
      iXiBig=ListRelXI_rho(iADCP, 1, idx);
      eSum=eSum+ANG_rho(iEtaBig, iXiBig)*eWeight;
      eSumWeight=eSumWeight+eWeight;
    end;
  end;
  ListAng(iADCP,1)=eSum/eSumWeight;
end;
clear('TotalArray');
TotalArray.nbADCP=nbADCP;
TotalArray.ListRelETA_u=ListRelETA_u;
TotalArray.ListRelXI_u=ListRelXI_u;
TotalArray.ListRelCoeff_u=ListRelCoeff_u;

TotalArray.ListRelETA_v=ListRelETA_v;
TotalArray.ListRelXI_v=ListRelXI_v;
TotalArray.ListRelCoeff_v=ListRelCoeff_v;

TotalArray.ListRelETA_rho=ListRelETA_rho;
TotalArray.ListRelXI_rho=ListRelXI_rho;
TotalArray.ListRelCoeff_rho=ListRelCoeff_rho;

TotalArray.ListAng=ListAng;
TotalArray.DEPbig_rho=DEP_rho;
TotalArray.DEPbig_u=DEPbig_u;
TotalArray.DEPbig_v=DEPbig_v;

TotalArray.AddiRecordVertDisc=AddiRecordVertDisc;
TotalArray.ADCPverticalLevels=ADCPverticalLevels;



