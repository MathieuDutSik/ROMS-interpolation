function [LON_psi2, LAT_psi2, MSK_psi2]=GRID_ExtendedPsiSec(...
    MSK_rho, LON_rho, LAT_rho, LON_psi, LAT_psi, ...
    LON_u, LAT_u, LON_v, LAT_v)

[eta_rho, xi_rho]=size(MSK_rho);
eta_psi=eta_rho-1;
xi_psi=xi_rho-1;
eta_psi2=eta_psi+2;
xi_psi2=xi_psi+2;
MSK_psi=MSK_rho(1:eta_rho-1,1:xi_rho-1).*...
        MSK_rho(2:eta_rho,1:xi_rho-1).*...
        MSK_rho(1:eta_rho-1,2:xi_rho).*...
        MSK_rho(2:eta_rho,2:xi_rho);



LAT_psi2=zeros(eta_psi+2, xi_psi+2);
LON_psi2=zeros(eta_psi+2, xi_psi+2);
for iEta=2:eta_psi2-1
  for iXi=2:xi_psi2-1
    LON_psi2(iEta, iXi)=LON_psi(iEta-1, iXi-1);
    LAT_psi2(iEta, iXi)=LAT_psi(iEta-1, iXi-1);
  end;
end;

for iEta=2:eta_psi+1
  LON_psi2(iEta, 1)=2*LON_v(iEta-1, 1)-LON_psi(iEta-1,1);
  LAT_psi2(iEta, 1)=2*LAT_v(iEta-1, 1)-LAT_psi(iEta-1,1);
end;
for iEta=2:eta_psi+1
  LON_psi2(iEta, xi_psi+2)=2*LON_v(iEta-1, xi_rho)-LON_psi(iEta-1,xi_psi);
  LAT_psi2(iEta, xi_psi+2)=2*LAT_v(iEta-1, xi_rho)-LAT_psi(iEta-1,xi_psi);
end;

for iXi=2:xi_psi+1
  LON_psi2(1, iXi)=2*LON_u(1,iXi-1)-LON_psi(1,iXi-1);
  LAT_psi2(1, iXi)=2*LAT_u(1,iXi-1)-LAT_psi(1,iXi-1);
end;
for iXi=2:xi_psi+1
  LON_psi2(eta_psi+2, iXi)=2*LON_u(eta_rho,iXi-1)-LON_psi(eta_psi,iXi-1);
  LAT_psi2(eta_psi+2, iXi)=2*LAT_u(eta_rho,iXi-1)-LAT_psi(eta_psi,iXi-1);
end;

LON_psi2(1,1)=2*LON_psi2(2,1)-LON_psi2(3,1);
LAT_psi2(1,1)=2*LAT_psi2(2,1)-LAT_psi2(3,1);

LON_psi2(eta_psi+2,1)=2*LON_psi2(eta_psi+1,1)-LON_psi2(eta_psi,1);
LAT_psi2(eta_psi+2,1)=2*LAT_psi2(eta_psi+1,1)-LAT_psi2(eta_psi,1);

LON_psi2(1,xi_psi+2)=2*LON_psi2(2,xi_psi+2)-LON_psi2(3,xi_psi+2);
LAT_psi2(1,xi_psi+2)=2*LAT_psi2(2,xi_psi+2)-LAT_psi2(3,xi_psi+2);

LON_psi2(eta_psi+2,xi_psi+2)=2*LON_psi2(eta_psi+1,xi_psi+2)-...
    LON_psi2(eta_psi,xi_psi+2);
LAT_psi2(eta_psi+2,xi_psi+2)=2*LAT_psi2(eta_psi+1,xi_psi+2)-...
    LAT_psi2(eta_psi,xi_psi+2);

MSK_psi2=zeros(eta_psi2, xi_psi2);
for iEta=2:eta_psi+1
  for iXi=2:xi_psi+1
    MSK_psi2(iEta, iXi)=MSK_psi(iEta-1, iXi-1);
  end;
end;
%
for iEta=2:eta_psi+1
  if (MSK_psi(iEta-1,1) == 1)
    MSK_psi2(iEta,1)=1;
  end;
  if (MSK_psi(iEta-1, xi_psi) == 1)
    MSK_psi2(iEta, xi_psi2)=1;
  end;
end;
for iXi=2:xi_psi+1
  if (MSK_psi(1, iXi-1) == 1)
    MSK_psi2(1,iXi)=1;
  end;
  if (MSK_psi(eta_psi, iXi-1) == 1)
    MSK_psi2(eta_psi2, iXi)=1;
  end;
end;
%
if (MSK_psi2(2,1) == 1 && MSK_psi2(1,2) == 1)
  MSK_psi2(1,1)=1;
end;
if (MSK_psi2(eta_psi2,2) == 1 && MSK_psi2(eta_psi2-1,1) == 1)
  MSK_psi2(eta_psi2,1)=1;
end;
if (MSK_psi2(eta_psi2, xi_psi2-1) == 1 && MSK_psi2(eta_psi2-1, xi_psi2) == 1)
  MSK_psi2(eta_psi2, xi_psi2)=1;
end;
if (MSK_psi2(1, xi_psi2-1) == 1 && MSK_psi2(2, xi_psi2) == 1)
  MSK_psi2(1, xi_psi2)=1;
end;
