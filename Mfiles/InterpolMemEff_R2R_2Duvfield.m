function [UBARsma, VBARsma]=InterpolMemEff_R2R_2Duvfield(TotalArray, ...
    UBARbig, VBARbig)
MSKbig_rho=TotalArray.MSKbig_rho;
MSKbig_u=TotalArray.MSKbig_u;
MSKbig_v=TotalArray.MSKbig_v;
ANGbig_rho=TotalArray.ANGbig_rho;
MSKsma_rho=TotalArray.MSKsma_rho;
MSKsma_u=TotalArray.MSKsma_u;
MSKsma_v=TotalArray.MSKsma_v;
ANGsma_rho=TotalArray.ANGsma_rho;
ANGsma_u_big=TotalArray.ANGsma_u_big;
ANGsma_v_big=TotalArray.ANGsma_v_big;
ListRelETA_uu=TotalArray.ListRelETA_uu;
ListRelXI_uu=TotalArray.ListRelXI_uu;
ListRelCoeff_uu=TotalArray.ListRelCoeff_uu;
ListRelETA_uv=TotalArray.ListRelETA_uv;
ListRelXI_uv=TotalArray.ListRelXI_uv;
ListRelCoeff_uv=TotalArray.ListRelCoeff_uv;
ListRelETA_vu=TotalArray.ListRelETA_vu;
ListRelXI_vu=TotalArray.ListRelXI_vu;
ListRelCoeff_vu=TotalArray.ListRelCoeff_vu;
ListRelETA_vv=TotalArray.ListRelETA_vv;
ListRelXI_vv=TotalArray.ListRelXI_vv;
ListRelCoeff_vv=TotalArray.ListRelCoeff_vv;

[eta_rho_big,xi_rho_big]=size(MSKbig_rho);
[eta_u_big,xi_u_big]=size(MSKbig_u);
[eta_v_big,xi_v_big]=size(MSKbig_v);
ANGbig_u=(ANGbig_rho(:, 1:xi_rho_big-1)+ANGbig_rho(:, 2:xi_rho_big))/2;
ANGbig_v=(ANGbig_rho(1:eta_rho_big-1, :)+ANGbig_rho(2:eta_rho_big, :))/2;
[eta_u_big1, xi_u_big1]=size(UBARbig);
[eta_v_big1, xi_v_big1]=size(VBARbig);
if (eta_u_big ~= eta_u_big1 || xi_u_big ~= xi_u_big1)
  disp('size error between TotalArray.MSKbig_u and UBARbig');
  keyboard;
end;
if (eta_v_big ~= eta_v_big1 || xi_v_big ~= xi_v_big1)
  disp('size error between TotalArray.MSKbig_v and VBARbig');
  keyboard;
end;
[eta_rho_sma,xi_rho_sma]=size(MSKsma_rho);
[eta_u_sma,xi_u_sma]=size(MSKsma_u);
[eta_v_sma,xi_v_sma]=size(MSKsma_v);
ANGsma_u=(ANGsma_rho(:, 1:xi_rho_sma-1)+ANGsma_rho(:, 2:xi_rho_sma))/2;
ANGsma_v=(ANGsma_rho(1:eta_rho_sma-1, :)+ANGsma_rho(2:eta_rho_sma, :))/2;
%
UBARsma=zeros(eta_u_sma, xi_u_sma);
for iEtaSma=1:eta_u_sma
  for iXiSma=1:xi_u_sma
    if (MSKsma_u(iEtaSma, iXiSma) == 1)
      UBAR_u=0;
      for idx=1:4
	eWeight=ListRelCoeff_uu(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_uu(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_uu(iEtaSma, iXiSma, idx);
	  UBAR_u=UBAR_u+UBARbig(iEtaBig, iXiBig)*eWeight;
	end;
      end;
      VBAR_u=0;
      for idx=1:4
	eWeight=ListRelCoeff_vu(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_vu(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_vu(iEtaSma, iXiSma, idx);
	  VBAR_u=VBAR_u+VBARbig(iEtaBig, iXiBig)*eWeight;
	end;
      end;
      deltaAng=ANGsma_u(iEtaSma, iXiSma)-...
	  ANGsma_u_big(iEtaSma, iXiSma);
      UBARsma(iEtaSma, iXiSma)=cos(deltaAng)*UBAR_u+...
	  sin(deltaAng)*VBAR_u;
    end;
  end;
end;
VBARsma=zeros(eta_v_sma, xi_v_sma);
for iEtaSma=1:eta_v_sma
  for iXiSma=1:xi_v_sma
    if (MSKsma_v(iEtaSma, iXiSma) == 1)
      UBAR_v=0;
      for idx=1:4
	eWeight=ListRelCoeff_uv(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_uv(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_uv(iEtaSma, iXiSma, idx);
	  UBAR_v=UBAR_v+UBARbig(iEtaBig, iXiBig)*eWeight;
	end;
      end;
      VBAR_v=0;
      for idx=1:4
	eWeight=ListRelCoeff_vv(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_vv(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_vv(iEtaSma, iXiSma, idx);
	  VBAR_v=VBAR_v+VBARbig(iEtaBig, iXiBig)*eWeight;
	end;
      end;
      deltaAng=ANGsma_v(iEtaSma, iXiSma)-...
	  ANGsma_v_big(iEtaSma, iXiSma);
      VBARsma(iEtaSma, iXiSma)=cos(deltaAng)*VBAR_v-...
	  sin(deltaAng)*UBAR_v;
    end;
  end;
end;
