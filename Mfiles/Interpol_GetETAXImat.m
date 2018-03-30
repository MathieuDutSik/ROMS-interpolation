function [ETAmat, XImat]=...
    Interpol_GetETAXImat(eta_rho, xi_rho)

ETAmat=zeros(eta_rho, xi_rho);
XImat=zeros(eta_rho, xi_rho);
for iEta=1:eta_rho
  for iXi=1:xi_rho
    ETAmat(iEta, iXi)=iEta;
    XImat(iEta, iXi)=iXi;
  end;
end;
