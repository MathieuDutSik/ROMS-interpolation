function [ListRelETA, ListRelXI, ListRelCoeff]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
	LONbig_rho, LATbig_rho, MSKbig_rho, ...
	LONsma_rho, LATsma_rho, MSKsma_rho);

KseaBig=find(MSKbig_rho == 1);
[eta_rho_big, xi_rho_big]=size(MSKbig_rho);
[ETAmat_big, XImat_big]=Interpol_GetETAXImat(...
    eta_rho_big, xi_rho_big);
LONbig_rho_Sel=LONbig_rho(KseaBig);
LATbig_rho_Sel=LATbig_rho(KseaBig);
ETAmat_big_Sel=ETAmat_big(KseaBig);
XImat_big_Sel=XImat_big(KseaBig);

[eta_rho_sma, xi_rho_sma]=size(MSKsma_rho);
[ListETA, ListXI, ListCoeff]=Interpol_R2R_GetCoefficients(...
    LONbig_rho, LATbig_rho, ...
    LONsma_rho, LATsma_rho);
ListBlock=[0 0;
	   1 0;
	   0 1;
	   1 1];
ListRelETA=zeros(eta_rho_sma, xi_rho_sma, 4);
ListRelXI=zeros(eta_rho_sma, xi_rho_sma, 4);
ListRelCoeff=zeros(eta_rho_sma, xi_rho_sma, 4);
for iEtaSma=1:eta_rho_sma
  for iXiSma=1:xi_rho_sma
    if (MSKsma_rho(iEtaSma, iXiSma) == 1)
      eLonSma=LONsma_rho(iEtaSma, iXiSma);
      eLatSma=LATsma_rho(iEtaSma, iXiSma);
      iEtaBig=ListETA(iEtaSma, iXiSma);
      iXiBig=ListXI(iEtaSma, iXiSma);
      if (iEtaBig == -1 && iXiBig == -1)
	disp(['Sorry, but you are trying to interpolate on wrong' ...
	      ' place']);
	keyboard;
      end;
      idx=0;
      eSumWeight=0;
      for iblock=1:4
	iEtaN=iEtaBig+ListBlock(iblock,1);
	iXiN=iXiBig+ListBlock(iblock,2);
	if (MSKbig_rho(iEtaN, iXiN) == 1)
	  idx=idx+1;
	  eWeight=ListCoeff(iEtaSma, iXiSma, iblock);
	  ListRelETA(iEtaSma, iXiSma, idx)=iEtaN;
	  ListRelXI(iEtaSma, iXiSma, idx)=iXiN;
	  eSumWeight=eSumWeight+eWeight;
	  ListRelCoeff(iEtaSma, iXiSma, idx)=eWeight;
	end;
      end;
      for i=1:idx
	ListRelCoeff(iEtaSma, iXiSma, i)=...
	    ListRelCoeff(iEtaSma, iXiSma, i)/eSumWeight;
      end;
      if (idx == 0)
	idx=nearxy(LONbig_rho_Sel, LATbig_rho_Sel, ...
		   eLonSma, eLatSma);
	iEta=ETAmat_big_Sel(idx,1);
	iXi=XImat_big_Sel(idx,1);
	ListRelETA(iEtaSma, iXiSma, 1)=iEta;
	ListRelXI(iEtaSma, iXiSma, 1)=iXi;
	ListRelCoeff(iEtaSma, iXiSma, 1)=1;
      end;
    end;
  end;
end;
