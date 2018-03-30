function RadiusInflMatrix=...
    GRID_GetMaxRadiusInfluence_kernel(...
	LON_rho, LAT_rho, MSK_rho, strDistType)
[eta_rho, xi_rho]=size(MSK_rho);
RadiusInflMatrix=zeros(eta_rho, xi_rho);
if (UTIL_IsStringEqual(strDistType, 'spherical') == 1)
  TheChoice=1;
elseif (UTIL_IsStringEqual(strDistType, 'euclidean') == 1)
  TheChoice=2;
else
  disp('You did not choose rightly! death on you');
  keyboard;
end;
%
ListDir=[1 1;
	 1 -1;
	 -1 1;
	 -1 -1];
for iEta=1:eta_rho
  for iXi=1:xi_rho
    ePair1(1)=LON_rho(iEta, iXi);
    ePair1(2)=LAT_rho(iEta, iXi);
    MaxDist=0;
    for ineigh=1:4
      iEtaN=iEta+ListDir(ineigh,1);
      iXiN=iXi+ListDir(ineigh,2);
      if (iEtaN <= eta_rho && iEtaN >= 1 && ...
	  iXiN <= xi_rho && iXiN >= 1)
	ePair2(1)=LON_rho(iEtaN, iXiN);
	ePair2(2)=LAT_rho(iEtaN, iXiN);
	if (TheChoice == 1)
	  dist=GeodesicDistance(ePair1, ePair2);
	else
	  dist=sqrt((ePair1(1)-ePair2(1))^2+(ePair1(2)-ePair2(2))^2);
	end;
	if (dist > MaxDist)
	  MaxDist=dist;
	end;
      end;
    end;
    RadiusInflMatrix(iEta, iXi)=MaxDist;
  end;
end;
