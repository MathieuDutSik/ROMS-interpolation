function [str1, str2, str3, str4, ...
	  str5, str6, str7, str8]=GRID_GetEWNSinfo(MSK_rho)
	 

[eta_rho, xi_rho]=size(MSK_rho);

nbWetSouth=0;
nbWetNorth=0;
for iXi=1:xi_rho
  if (MSK_rho(1, iXi) == 1)
    nbWetSouth=nbWetSouth+1;
  end;
  if (MSK_rho(eta_rho, iXi) == 1)
    nbWetNorth=nbWetNorth+1;
  end;
end;
nbWetEast=0;
nbWetWest=0;
for iEta=1:eta_rho
  if (MSK_rho(iEta, 1) == 1)
    nbWetWest=nbWetWest+1;
  end;
  if (MSK_rho(iEta, xi_rho) == 1)
    nbWetEast=nbWetEast+1;
  end;
end;
str1=['(:,xi_rho)      East : ' num2str(nbWetEast) ' wet points'];
str2=['(:,1)           West : ' num2str(nbWetWest) ' wet points'];
str3=['(eta_rho,:)    North : ' num2str(nbWetNorth) ' wet points'];
str4=['(1,:)          South : ' num2str(nbWetSouth) ' wet points'];
%
pos=0;
while(1)
  nbSea=sum(MSK_rho(1+pos,:));
  if (nbSea > 0)
    break;
  end;
  if (pos == eta_rho-1)
    break;
  end;
  pos=pos+1;
end;
str5=['boundary of size ' num2str(pos) ' in the south'];
%
pos=0;
while(1)
  nbSea=sum(MSK_rho(:,1+pos));
  if (nbSea > 0)
    break;
  end;
  if (pos == xi_rho-1)
    break;
  end;
  pos=pos+1;
end;
str6=['boundary of size ' num2str(pos) ' in the east'];
%
pos=0;
while(1)
  nbSea=sum(MSK_rho(eta_rho-pos,:));
  if (nbSea > 0)
    break;
  end;
  if (pos == eta_rho-1)
    break;
  end;
  pos=pos+1;
end;
str7=['boundary of size ' num2str(pos) ' in the north'];
%
pos=0;
while(1)
  nbSea=sum(MSK_rho(:,xi_rho-pos));
  if (nbSea > 0)
    break;
  end;
  if (pos == xi_rho-1)
    break;
  end;
  pos=pos+1;
end;
str8=['boundary of size ' num2str(pos) ' in the west'];
