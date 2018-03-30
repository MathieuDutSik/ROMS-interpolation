function [nbPositive, nbNegative]=GRID_GridOrientationInfo(LON_rho, LAT_rho)

% GRID_GridOrientationInfo(LON_rho, LAT_rho)
%
% this is a diagnostic function that prints result on the
% orientation of the grid (presumably all positive or all negative)
[eta_rho, xi_rho]=size(LON_rho);
nbPositive=0;
nbNegative=0;
eta_psi=eta_rho-1;
xi_psi=xi_rho-1;
TheMat=zeros(2,2);
for iEta=1:eta_psi
  for iXi=1:xi_psi
    TheMat(1,1)=LON_rho(iEta+1, iXi)-LON_rho(iEta, iXi);
    TheMat(2,1)=LAT_rho(iEta+1, iXi)-LAT_rho(iEta, iXi);
    TheMat(1,2)=LON_rho(iEta, iXi+1)-LON_rho(iEta, iXi);
    TheMat(2,2)=LAT_rho(iEta, iXi+1)-LAT_rho(iEta, iXi);
    TheDet=det(TheMat);
    if (TheDet>0)
      nbPositive=nbPositive+1;
    else
      nbNegative=nbNegative+1;
    end;
  end;
end;
