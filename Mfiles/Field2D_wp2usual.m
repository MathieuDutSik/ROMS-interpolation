function ZETAusual=Field2D_wp2usual(MSK_rho, ZETAwp)

[eta_rho, xi_rho]=size(MSK_rho);
ZETAusual=zeros(eta_rho, xi_rho);
idx=0;
for iEta=1:eta_rho
  for iXi=1:xi_rho
    if (MSK_rho(iEta, iXi) == 1)
      idx=idx+1;
      ZETAusual(iEta,iXi)=ZETAwp(1,idx);
    end;
  end;
end;
