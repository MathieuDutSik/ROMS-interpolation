function ZETAsma=InterpolMemEff_R2R_2Dfield(TotalArray, ZETAbig)
MSKbig_rho=TotalArray.MSKbig_rho;
MSKsma_rho=TotalArray.MSKsma_rho;
ListRelETA_rr=TotalArray.ListRelETA_rr;
ListRelXI_rr=TotalArray.ListRelXI_rr;
ListRelCoeff_rr=TotalArray.ListRelCoeff_rr;
[eta_rho_big, xi_rho_big]=size(MSKbig_rho);
[eta_rho_big1, xi_rho_big1]=size(ZETAbig);
if (eta_rho_big1 ~= eta_rho_big || xi_rho_big1 ~= xi_rho_big)
  disp('We have a size error between TotalArray and ZETAbig');
  keyboard;
end;
[eta_rho_sma, xi_rho_sma]=size(MSKsma_rho);
ZETAsma=zeros(eta_rho_sma, xi_rho_sma);
for iEtaSma=1:eta_rho_sma
  for iXiSma=1:xi_rho_sma
    if (MSKsma_rho(iEtaSma, iXiSma) == 1)
      eSum=0;
      for idx=1:4
	eWeight=ListRelCoeff_rr(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_rr(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_rr(iEtaSma, iXiSma, idx);
	  eSum=eSum+ZETAbig(iEtaBig, iXiBig)*eWeight;
	end;
      end;
      ZETAsma(iEtaSma, iXiSma)=eSum;
    end;
  end;
end;
