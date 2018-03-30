function [UVbar_east, UVbar_north]=ADCP_R2A_Interpolation2Duvfield(...
    TotalArray, UBARroms, VBARroms)
%
% [UVbar_east, UVbar_north]=ADCP_R2A_Interpolation2Duvfield(...
%    TotalArray, UBARroms, VBARroms)
%
% UVbar_east is the eastern component of barotropic velocity
% UVbar_north is the northern component of barotropic velocity
% 
% TotalArray is obtained by the function ADCP_R2A_GetTotalArray
% UBARroms is the u-component of barotropic velocity
% VBARroms is the v-component of barotropic velocity
%

nbADCP=TotalArray.nbADCP;

ListRelETA_u=TotalArray.ListRelETA_u;
ListRelXI_u=TotalArray.ListRelXI_u;
ListRelCoeff_u=TotalArray.ListRelCoeff_u;

ListRelETA_v=TotalArray.ListRelETA_v;
ListRelXI_v=TotalArray.ListRelXI_v;
ListRelCoeff_v=TotalArray.ListRelCoeff_v;


UVbar_east=zeros(nbADCP,1);
UVbar_north=zeros(nbADCP,1);
for iADCP=1:nbADCP
  eSum=0;
  eSumWeight=0;
  for idx=1:4
    eWeight=TotalArray.ListRelCoeff_u(iADCP, 1, idx);
    if (eWeight > 0)
      iEtaBig=TotalArray.ListRelETA_u(iADCP, 1, idx);
      iXiBig=TotalArray.ListRelXI_u(iADCP, 1, idx);
      eSum=eSum+UBARroms(iEtaBig, iXiBig)*eWeight;
      eSumWeight=eSumWeight+eWeight;
    end;
  end;
  UVBAR_u=eSum/eSumWeight;
  eSum=0;
  eSumWeight=0;
  for idx=1:4
    eWeight=TotalArray.ListRelCoeff_v(iADCP, 1, idx);
    if (eWeight > 0)
      iEtaBig=TotalArray.ListRelETA_v(iADCP, 1, idx);
      iXiBig=TotalArray.ListRelXI_v(iADCP, 1, idx);
      eSum=eSum+VBARroms(iEtaBig, iXiBig)*eWeight;
      eSumWeight=eSumWeight+eWeight;
    end;
  end;
  UVBAR_v=eSum/eSumWeight;
  TheAng=TotalArray.ListAng(iADCP,1);
  UVbarEast=cos(TheAng)*UVBAR_u-sin(TheAng)*UVBAR_v;
  UVbarNorth=sin(TheAng)*UVBAR_u+cos(TheAng)*UVBAR_v;
  UVbar_east(iADCP,1)=UVbarEast;
  UVbar_north(iADCP,1)=UVbarNorth;
end;
