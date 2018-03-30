function [Usma, Vsma]=InterpolMemEff_R2R_3Duvfield(...
    TotalArray, AddiRecordInfo, Ubig, Vbig)
MSKbig_rho=TotalArray.MSKbig_rho;
DEPbig_rho=TotalArray.DEPbig_rho;
MSKbig_u=TotalArray.MSKbig_u;
MSKbig_v=TotalArray.MSKbig_v;
ANGbig_rho=TotalArray.ANGbig_rho;
MSKsma_rho=TotalArray.MSKsma_rho;
DEPsma_rho=TotalArray.DEPsma_rho;
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
%
Nbig=AddiRecordInfo.Nbig;
BigThetaS=AddiRecordInfo.BigThetaS;
BigThetaB=AddiRecordInfo.BigThetaB;
Nsma=AddiRecordInfo.Nsma;
SmaThetaS=AddiRecordInfo.SmaThetaS;
SmaThetaB=AddiRecordInfo.SmaThetaB;
%
[eta_rho_big,xi_rho_big]=size(MSKbig_rho);
[eta_u_big,xi_u_big]=size(MSKbig_u);
[eta_v_big,xi_v_big]=size(MSKbig_v);
%
[NbigU1, eta_u_big1, xi_u_big1]=size(Ubig);
[NbigV1, eta_v_big1, xi_v_big1]=size(Vbig);
if (eta_u_big ~= eta_u_big1 || xi_u_big ~= xi_u_big1 || ...
    NbigU1 ~= Nbig)
  disp('size error between TotalArray.MSKbig_u and Ubig');
  keyboard;
end;
if (eta_v_big ~= eta_v_big1 || xi_v_big ~= xi_v_big1 || ...
    NbigV1 ~= Nbig)
  disp('size error between TotalArray.MSKbig_v and Vbig');
  keyboard;
end;
ANGbig_u=(ANGbig_rho(:, 1:xi_rho_big-1)+ANGbig_rho(:, 2:xi_rho_big))/2;
ANGbig_v=(ANGbig_rho(1:eta_rho_big-1, :)+ANGbig_rho(2:eta_rho_big, :))/2;
DEPbig_u=(DEPbig_rho(:, 1:xi_rho_big-1)+DEPbig_rho(:, 2:xi_rho_big))/2;
DEPbig_v=(DEPbig_rho(1:eta_rho_big-1, :)+DEPbig_rho(2:eta_rho_big, :))/2;
[BigSc_w, BigCs_w, BigSc_r, BigCs_r]=GRID_GetSc_Cs(...
    Nbig, BigThetaS, BigThetaB);
Bigcff_r=AddiRecordInfo.Bighc*(BigSc_r-BigCs_r);
SIZ_u_big=eta_u_big*xi_u_big*Nbig;
SIZ_v_big=eta_v_big*xi_v_big*Nbig;
SIZ_big=SIZ_u_big+SIZ_v_big;
%
%
[eta_rho_sma,xi_rho_sma]=size(MSKsma_rho);
[eta_u_sma,xi_u_sma]=size(MSKsma_u);
[eta_v_sma,xi_v_sma]=size(MSKsma_v);
ANGsma_u=(ANGsma_rho(:, 1:xi_rho_sma-1)+ANGsma_rho(:, 2:xi_rho_sma))/2;
ANGsma_v=(ANGsma_rho(1:eta_rho_sma-1, :)+ANGsma_rho(2:eta_rho_sma, :))/2;
DEPsma_u=(DEPsma_rho(:, 1:xi_rho_sma-1)+DEPsma_rho(:, 2:xi_rho_sma))/2;
DEPsma_v=(DEPsma_rho(1:eta_rho_sma-1, :)+DEPsma_rho(2:eta_rho_sma, :))/2;
[SmaSc_w, SmaCs_w, SmaSc_r, SmaCs_r]=GRID_GetSc_Cs(...
    Nsma, SmaThetaS, SmaThetaB);
Smacff_r=AddiRecordInfo.Smahc*(SmaSc_r-SmaCs_r);
SIZ_u_sma=eta_u_sma*xi_u_sma*Nsma;
SIZ_v_sma=eta_v_sma*xi_v_sma*Nsma;
SIZ_sma=SIZ_u_sma+SIZ_v_sma;
%
%
ListDepth=zeros(4,1);
Usma=zeros(Nsma, eta_u_sma,xi_u_sma);
for iEtaSma=1:eta_u_sma
  for iXiSma=1:xi_u_sma
    if (MSKsma_u(iEtaSma, iXiSma) == 1)
      deltaAng=ANGsma_u(iEtaSma, iXiSma)-...
	       ANGsma_u_big(iEtaSma, iXiSma);
      ROMSdepthSma=Smacff_r+SmaCs_r*DEPsma_u(iEtaSma, iXiSma);
      for iNsma=1:Nsma
	eSumU_u=0;
	eSumWeight=0;
	MyDep=ROMSdepthSma(1,iNsma);
	for idx=1:4
	  ListDepth(idx,1)=0;
	  eWeight=ListRelCoeff_uu(iEtaSma, iXiSma, idx);
	  if (eWeight > 0)
	    iEtaBig=ListRelETA_uu(iEtaSma, iXiSma, idx);
	    iXiBig=ListRelXI_uu(iEtaSma, iXiSma, idx);
	    depbig=DEPbig_u(iEtaBig, iXiBig);
	    ListDepth(idx,1)=depbig;
	    ROMSdepthBig=Bigcff_r+BigCs_r*depbig;
	    test=0;
	    for iNbig=1:Nbig-1
	      dep1=ROMSdepthBig(1,iNbig);
	      dep2=ROMSdepthBig(1,iNbig+1);
	      if (dep1 <= MyDep && MyDep <= dep2)
		test=1;
		alpha1=(dep2-MyDep)/(dep2-dep1);
		alpha2=(MyDep-dep1)/(dep2-dep1);
		eSumU_u=eSumU_u+eWeight*alpha1*...
			Ubig(iNbig, iEtaBig, iXiBig);
		eSumU_u=eSumU_u+eWeight*alpha2*...
			Ubig(iNbig+1, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		eSumU_u=eSumU_u+eWeight*Ubig(Nbig, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  U_u=eSumU_u/eSumWeight;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_uu(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_uu(iEtaSma, iXiSma, idxdep);
	  U_u=Ubig(1, iEtaBig, iXiBig);
	end;
	%
	eSumWeight=0;
	eSumV_u=0;
	for idx=1:4
	  ListDepth(idx,1)=0;
	  eWeight=ListRelCoeff_vu(iEtaSma, iXiSma, idx);
	  if (eWeight > 0)
	    iEtaBig=ListRelETA_vu(iEtaSma, iXiSma, idx);
	    iXiBig=ListRelXI_vu(iEtaSma, iXiSma, idx);
	    depbig=DEPbig_v(iEtaBig, iXiBig);
	    ListDepth(idx,1)=depbig;
	    ROMSdepthBig=Bigcff_r+BigCs_r*depbig;
	    test=0;
	    for iNbig=1:Nbig-1
	      dep1=ROMSdepthBig(1,iNbig);
	      dep2=ROMSdepthBig(1,iNbig+1);
	      if (dep1 <= MyDep && MyDep <= dep2)
		test=1;
		alpha1=(dep2-MyDep)/(dep2-dep1);
		alpha2=(MyDep-dep1)/(dep2-dep1);
		eSumV_u=eSumV_u+eWeight*alpha1*...
			Vbig(iNbig, iEtaBig, iXiBig);
		eSumV_u=eSumV_u+eWeight*alpha2*...
			Vbig(iNbig+1, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		eSumV_u=eSumV_u+eWeight*Vbig(Nbig, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  V_u=eSumV_u/eSumWeight;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_vu(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_vu(iEtaSma, iXiSma, idxdep);
	  V_u=Vbig(1, iEtaBig, iXiBig);
	end;
	TheVal=U_u*cos(deltaAng)+V_u*sin(deltaAng);
	Usma(iNsma, iEtaSma, iXiSma)=TheVal;
      end;
    end;
  end;
end;
Vsma=zeros(Nsma, eta_v_sma,xi_v_sma);
for iEtaSma=1:eta_v_sma
  for iXiSma=1:xi_v_sma
    if (MSKsma_v(iEtaSma, iXiSma) == 1)
      deltaAng=ANGsma_v(iEtaSma, iXiSma)-...
	       ANGsma_v_big(iEtaSma, iXiSma);
      ROMSdepthSma=Smacff_r+SmaCs_r*DEPsma_v(iEtaSma, iXiSma);
      for iNsma=1:Nsma
	MyDep=ROMSdepthSma(1,iNsma);
	eSumU_v=0;
	eSumWeight=0;
	nbU_v=0;
	for idx=1:4
	  ListDepth(idx,1)=0;
	  eWeight=ListRelCoeff_uv(iEtaSma, iXiSma, idx);
	  if (eWeight > 0)
	    iEtaBig=ListRelETA_uv(iEtaSma, iXiSma, idx);
	    iXiBig=ListRelXI_uv(iEtaSma, iXiSma, idx);
	    depbig=DEPbig_u(iEtaBig, iXiBig);
	    ListDepth(idx,1)=depbig;
	    ROMSdepthBig=Bigcff_r+BigCs_r*depbig;
	    test=0;
	    for iNbig=1:Nbig-1
	      dep1=ROMSdepthBig(1,iNbig);
	      dep2=ROMSdepthBig(1,iNbig+1);
	      if (dep1 <= MyDep && MyDep <= dep2)
		test=1;
		alpha1=(dep2-MyDep)/(dep2-dep1);
		alpha2=(MyDep-dep1)/(dep2-dep1);
		eSumU_v=eSumU_v+eWeight*alpha1*...
			Ubig(iNbig, iEtaBig, iXiBig);
		eSumU_v=eSumU_v+eWeight*alpha2*...
			Ubig(iNbig+1, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		eSumU_v=eSumU_v+eWeight*Ubig(Nbig, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  U_v=eSumU_v/eSumWeight;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_uv(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_uv(iEtaSma, iXiSma, idxdep);
	  U_v=Ubig(1, iEtaBig, iXiBig);
	end;
	%
	eSumWeight=0;
	eSumV_v=0;
	for idx=1:4
	  ListDepth(idx,1)=0;
	  eWeight=ListRelCoeff_vv(iEtaSma, iXiSma, idx);
	  if (eWeight > 0)
	    iEtaBig=ListRelETA_vv(iEtaSma, iXiSma, idx);
	    iXiBig=ListRelXI_vv(iEtaSma, iXiSma, idx);
	    depbig=DEPbig_v(iEtaBig, iXiBig);
	    ListDepth(idx,1)=depbig;
	    ROMSdepthBig=Bigcff_r+BigCs_r*depbig;
	    test=0;
	    for iNbig=1:Nbig-1
	      dep1=ROMSdepthBig(1,iNbig);
	      dep2=ROMSdepthBig(1,iNbig+1);
	      if (dep1 <= MyDep && MyDep <= dep2)
		test=1;
		alpha1=(dep2-MyDep)/(dep2-dep1);
		alpha2=(MyDep-dep1)/(dep2-dep1);
		eSumV_v=eSumV_v+eWeight*alpha1*...
			Vbig(iNbig, iEtaBig, iXiBig);
		eSumV_v=eSumV_v+eWeight*alpha2*...
			Vbig(iNbig+1, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		eSumV_v=eSumV_v+eWeight*Vbig(Nbig, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  V_v=eSumV_v/eSumWeight;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_vv(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_vv(iEtaSma, iXiSma, idxdep);
	  V_v=Vbig(1, iEtaBig, iXiBig);
	end;
	Vsma(iNsma, iEtaSma, iXiSma)=-U_v*sin(deltaAng)+V_v*cos(deltaAng);
      end;
    end;
  end;
end;
