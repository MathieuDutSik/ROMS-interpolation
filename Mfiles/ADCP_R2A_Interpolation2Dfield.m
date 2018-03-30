function ZETAadcp=ADCP_R2A_Interpolation2Dfield(TotalArray, ZETAroms)
%
% ZETAadcp=ADCP_R2A_Interpolation2Dfield(TotalArray, ZETAroms)
%
% ZETAadcp is the height at the specified ADCP in output
% TotalArray is obtained by the function ADCP_R2A_GetTotalArray
% ZETAroms is the sea surface height of the ROMS model.
%

nbADCP=TotalArray.nbADCP;
ZETAadcp=zeros(nbADCP, 1);
for iADCP=1:nbADCP
  eSum=0;
  eSumWeight=0;
  for idx=1:4
    eWeight=TotalArray.ListRelCoeff_rho(iADCP, 1, idx);
    if (eWeight > 0)
      iEtaBig=TotalArray.ListRelETA_rho(iADCP, 1, idx);
      iXiBig=TotalArray.ListRelXI_rho(iADCP, 1, idx);
      eSum=eSum+ZETAroms(iEtaBig, iXiBig)*eWeight;
      eSumWeight=eSumWeight+eWeight;
    end;
  end;
  ZETAadcp(iADCP,1)=eSum/eSumWeight;
end;
