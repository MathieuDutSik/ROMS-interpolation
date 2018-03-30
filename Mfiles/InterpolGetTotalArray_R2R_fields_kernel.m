function TotalArray=InterpolGetTotalArray_R2R_fields_kernel(...
    GrdArrBig, GrdArrSma)

LONbig_rho=GrdArrBig.LON_rho;
LATbig_rho=GrdArrBig.LAT_rho;
MSKbig_rho=GrdArrBig.MSK_rho;
LONbig_u=GrdArrBig.LON_u;
LATbig_u=GrdArrBig.LAT_u;
MSKbig_u=GrdArrBig.MSK_u;
LONbig_v=GrdArrBig.LON_v;
LATbig_v=GrdArrBig.LAT_v;
MSKbig_v=GrdArrBig.MSK_v;
ANGbig_rho=GrdArrBig.ANG_rho;
DEPbig_rho=GrdArrBig.DEP_rho;
%
LONsma_rho=GrdArrSma.LON_rho;
LATsma_rho=GrdArrSma.LAT_rho;
MSKsma_rho=GrdArrSma.MSK_rho;
LONsma_u=GrdArrSma.LON_u;
LATsma_u=GrdArrSma.LAT_u;
MSKsma_u=GrdArrSma.MSK_u;
LONsma_v=GrdArrSma.LON_v;
LATsma_v=GrdArrSma.LAT_v;
MSKsma_v=GrdArrSma.MSK_v;
ANGsma_rho=GrdArrSma.ANG_rho;
DEPsma_rho=GrdArrSma.DEP_rho;
%
[ListRelETA_rr, ListRelXI_rr, ListRelCoeff_rr]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
        LONbig_rho, LATbig_rho, MSKbig_rho, ...
        LONsma_rho, LATsma_rho, MSKsma_rho);
%
[ListRelETA_uu, ListRelXI_uu, ListRelCoeff_uu]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
        LONbig_u, LATbig_u, MSKbig_u, ...
        LONsma_u, LATsma_u, MSKsma_u);
[ListRelETA_uv, ListRelXI_uv, ListRelCoeff_uv]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
        LONbig_u, LATbig_u, MSKbig_u, ...
        LONsma_v, LATsma_v, MSKsma_v);
[ListRelETA_vu, ListRelXI_vu, ListRelCoeff_vu]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
        LONbig_v, LATbig_v, MSKbig_v, ...
        LONsma_u, LATsma_u, MSKsma_u);
[ListRelETA_vv, ListRelXI_vv, ListRelCoeff_vv]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
        LONbig_v, LATbig_v, MSKbig_v, ...
        LONsma_v, LATsma_v, MSKsma_v);
%
[ListRelETA_ru, ListRelXI_ru, ListRelCoeff_ru]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
        LONbig_rho, LATbig_rho, MSKbig_rho, ...
        LONsma_u, LATsma_u, MSKsma_u);
[ListRelETA_rv, ListRelXI_rv, ListRelCoeff_rv]=...
    Interpol_R2R_FindCoefficient2Dgrid(...
        LONbig_rho, LATbig_rho, MSKbig_rho, ...
        LONsma_v, LATsma_v, MSKsma_v);
%
[eta_rho_sma, xi_rho_sma]=size(LONsma_rho);
[eta_u_sma,xi_u_sma]=size(LONsma_u);
[eta_v_sma,xi_v_sma]=size(LONsma_v);
ANGsma_u_big=zeros(eta_u_sma, xi_u_sma);
for iEtaSma=1:eta_u_sma
  for iXiSma=1:xi_u_sma
    eSum=0;
    for idx=1:4
      eWeight=ListRelCoeff_ru(iEtaSma, iXiSma, idx);
      if (eWeight > 0)
	iEtaBig=ListRelETA_ru(iEtaSma, iXiSma, idx);
	iXiBig=ListRelXI_ru(iEtaSma, iXiSma, idx);
	eSum=eSum+ANGbig_rho(iEtaBig, iXiBig)*eWeight;
      end;
    end;
    ANGsma_u_big(iEtaSma, iXiSma)=eSum;
  end;
end;
ANGsma_v_big=zeros(eta_v_sma, xi_v_sma);
for iEtaSma=1:eta_v_sma
  for iXiSma=1:xi_v_sma
    eSum=0;
    for idx=1:4
      eWeight=ListRelCoeff_rv(iEtaSma, iXiSma, idx);
      if (eWeight > 0)
	iEtaBig=ListRelETA_rv(iEtaSma, iXiSma, idx);
	iXiBig=ListRelXI_rv(iEtaSma, iXiSma, idx);
	eSum=eSum+ANGbig_rho(iEtaBig, iXiBig)*eWeight;
      end;
    end;
    ANGsma_v_big(iEtaSma, iXiSma)=eSum;
  end;
end;
%
TotalArray.ListRelETA_rr=ListRelETA_rr;
TotalArray.ListRelXI_rr=ListRelXI_rr;
TotalArray.ListRelCoeff_rr=ListRelCoeff_rr;
TotalArray.ListRelETA_uu=ListRelETA_uu;
TotalArray.ListRelXI_uu=ListRelXI_uu;
TotalArray.ListRelCoeff_uu=ListRelCoeff_uu;
TotalArray.ListRelETA_uv=ListRelETA_uv;
TotalArray.ListRelXI_uv=ListRelXI_uv;
TotalArray.ListRelCoeff_uv=ListRelCoeff_uv;
TotalArray.ListRelETA_vu=ListRelETA_vu;
TotalArray.ListRelXI_vu=ListRelXI_vu;
TotalArray.ListRelCoeff_vu=ListRelCoeff_vu;
TotalArray.ListRelETA_vv=ListRelETA_vv;
TotalArray.ListRelXI_vv=ListRelXI_vv;
TotalArray.ListRelCoeff_vv=ListRelCoeff_vv;
TotalArray.ListRelETA_ru=ListRelETA_ru;
TotalArray.ListRelXI_ru=ListRelXI_ru;
TotalArray.ListRelCoeff_ru=ListRelCoeff_ru;
TotalArray.ListRelETA_rv=ListRelETA_rv;
TotalArray.ListRelXI_rv=ListRelXI_rv;
TotalArray.ListRelCoeff_rv=ListRelCoeff_rv;
TotalArray.MSKbig_rho=MSKbig_rho;
TotalArray.MSKbig_u=MSKbig_u;
TotalArray.MSKbig_v=MSKbig_v;
TotalArray.MSKsma_rho=MSKsma_rho;
TotalArray.MSKsma_u=MSKsma_u;
TotalArray.MSKsma_v=MSKsma_v;
TotalArray.ANGbig_rho=ANGbig_rho;
TotalArray.ANGsma_rho=ANGsma_rho;
TotalArray.ANGsma_u_big=ANGsma_u_big;
TotalArray.ANGsma_v_big=ANGsma_v_big;
TotalArray.DEPbig_rho=DEPbig_rho;
TotalArray.DEPsma_rho=DEPsma_rho;
TotalArray.eta_rho_sma=eta_rho_sma;
TotalArray.xi_rho_sma=xi_rho_sma;
TotalArray.eta_u_sma=eta_u_sma;
TotalArray.xi_u_sma=xi_u_sma;
TotalArray.eta_v_sma=eta_v_sma;
TotalArray.xi_v_sma=xi_v_sma;
