function ArrayBigSma3Duv=InterpolGetSpMat_R2R_3Duvfield(...
    MSKbig_rho, DEPbig_rho, MSKbig_u, ...
    MSKbig_v, ANGbig_rho, ...
    MSKsma_rho, DEPsma_rho, MSKsma_u, ...
    MSKsma_v, ANGsma_rho, ...
    ANGsma_u_big, ANGsma_v_big, ...
    ListRelETA_uu, ListRelXI_uu, ListRelCoeff_uu, ...
    ListRelETA_uv, ListRelXI_uv, ListRelCoeff_uv, ...
    ListRelETA_vu, ListRelXI_vu, ListRelCoeff_vu, ...
    ListRelETA_vv, ListRelXI_vv, ListRelCoeff_vv, ...
    AddiRecordInfo)
disp('Computing sparse matrix for 2D + 1D UV-field interpolation');
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
ANGbig_u=(ANGbig_rho(:, 1:xi_rho_big-1)+ANGbig_rho(:, 2:xi_rho_big))/2;
ANGbig_v=(ANGbig_rho(1:eta_rho_big-1, :)+ANGbig_rho(2:eta_rho_big, :))/2;
DEPbig_u=(DEPbig_rho(:, 1:xi_rho_big-1)+DEPbig_rho(:, 2:xi_rho_big))/2;
DEPbig_v=(DEPbig_rho(1:eta_rho_big-1, :)+DEPbig_rho(2:eta_rho_big, :))/2;
[BigSc_w, BigCs_w, BigSc_r, BigCs_r]=GRID_GetSc_Cs(...
    Nbig, BigThetaS, BigThetaB);
Bigcff_r=AddiRecordInfo.Bighc*(BigSc_r-BigCs_r);
CorrespMatBig_u=GetCorrespMatrix(3, [Nbig; eta_u_big; xi_u_big]);
CorrespMatBig_v=GetCorrespMatrix(3, [Nbig; eta_v_big; xi_v_big]);
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
CorrespMatSma_u=GetCorrespMatrix(3, [Nsma; eta_u_sma; xi_u_sma]);
CorrespMatSma_v=GetCorrespMatrix(3, [Nsma; eta_v_sma; xi_v_sma]);
SIZ_u_sma=eta_u_sma*xi_u_sma*Nsma;
SIZ_v_sma=eta_v_sma*xi_v_sma*Nsma;
SIZ_sma=SIZ_u_sma+SIZ_v_sma;
%
UpperEstPoint=eta_u_sma*xi_u_sma*16*Nsma+eta_v_sma*xi_v_sma*16*Nsma;
iList=zeros(UpperEstPoint, 1);
jList=zeros(UpperEstPoint, 1);
sList=zeros(UpperEstPoint, 1);
idxspm=0;
%
ListDepth=zeros(4,1);
iListU_u=zeros(8,1);
jListU_u=zeros(8,1);
sListU_u=zeros(8,1);
iListV_u=zeros(8,1);
jListV_u=zeros(8,1);
sListV_u=zeros(8,1);
iListU_v=zeros(8,1);
jListU_v=zeros(8,1);
sListU_v=zeros(8,1);
iListV_v=zeros(8,1);
jListV_v=zeros(8,1);
sListV_v=zeros(8,1);
for iEtaSma=1:eta_u_sma
  for iXiSma=1:xi_u_sma
    if (MSKsma_u(iEtaSma, iXiSma) == 1)
      deltaAng=ANGsma_u(iEtaSma, iXiSma)-...
	       ANGsma_u_big(iEtaSma, iXiSma);
      ROMSdepthSma=Smacff_r+SmaCs_r*DEPsma_u(iEtaSma, iXiSma);
      for iNsma=1:Nsma
	idxSma_u=CorrespMatSma_u(iNsma, iEtaSma, iXiSma);
	eSumWeight=0;
	MyDep=ROMSdepthSma(1,iNsma);
	nbU_u=0;
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
		idxbig1=CorrespMatBig_u(iNbig, iEtaBig, iXiBig);
		idxbig2=CorrespMatBig_u(iNbig+1, iEtaBig, iXiBig);
		nbU_u=nbU_u+1;
		iListU_u(nbU_u,1)=idxSma_u;
		jListU_u(nbU_u,1)=idxbig1;
		sListU_u(nbU_u,1)=eWeight*alpha1;
		nbU_u=nbU_u+1;
		iListU_u(nbU_u,1)=idxSma_u;
		jListU_u(nbU_u,1)=idxbig2;
		sListU_u(nbU_u,1)=eWeight*alpha2;
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		idxbig=CorrespMatBig_u(Nbig, iEtaBig, iXiBig);
		nbU_u=nbU_u+1;
		iListU_u(nbU_u,1)=idxSma_u;
		jListU_u(nbU_u,1)=idxbig;
		sListU_u(nbU_u,1)=eWeight;
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  for i=1:nbU_u
	    sListU_u(i,1)=sListU_u(i,1)/eSumWeight;
	  end;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_uu(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_uu(iEtaSma, iXiSma, idxdep);
	  idxbig=CorrespMatBig_u(1, iEtaBig, iXiBig);
	  nbU_u=nbU_u+1;
	  iListU_u(nbU_u,1)=idxSma_u;
	  jListU_u(nbU_u,1)=idxbig;
	  sListU_u(nbU_u,1)=1;
	end;
	%
	eSumWeight=0;
	nbV_u=0;
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
		idxbig1=CorrespMatBig_v(iNbig, iEtaBig, iXiBig);
		idxbig2=CorrespMatBig_v(iNbig+1, iEtaBig, iXiBig);
		nbV_u=nbV_u+1;
		iListV_u(nbV_u,1)=idxSma_u;
		jListV_u(nbV_u,1)=SIZ_u_big+idxbig1;
		sListV_u(nbV_u,1)=eWeight*alpha1;
		nbV_u=nbV_u+1;
		iListV_u(nbV_u,1)=idxSma_u;
		jListV_u(nbV_u,1)=SIZ_u_big+idxbig2;
		sListV_u(nbV_u,1)=eWeight*alpha2;
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		idxbig=CorrespMatBig_v(Nbig, iEtaBig, iXiBig);
		nbV_u=nbV_u+1;
		iListV_u(nbV_u,1)=idxSma_u;
		jListV_u(nbV_u,1)=SIZ_u_big+idxbig;
		sListV_u(nbV_u,1)=eWeight;
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  for i=1:nbV_u
	    sListV_u(i,1)=sListV_u(i,1)/eSumWeight;
	  end;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_vu(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_vu(iEtaSma, iXiSma, idxdep);
	  idxbig=CorrespMatBig_v(1, iEtaBig, iXiBig);
	  nbV_u=nbV_u+1;
	  iListV_u(nbV_u,1)=idxSma_u;
	  jListV_u(nbV_u,1)=SIZ_u_big+idxbig;
	  sListV_u(nbV_u,1)=1;
	end;
	for i=1:nbU_u
	  idxspm=idxspm+1;
	  iList(idxspm,1)=iListU_u(i,1);
	  jList(idxspm,1)=jListU_u(i,1);
	  sList(idxspm,1)=sListU_u(i,1)*cos(deltaAng);
	end;
	for i=1:nbV_u
	  idxspm=idxspm+1;
	  iList(idxspm,1)=iListV_u(i,1);
	  jList(idxspm,1)=jListV_u(i,1);
	  sList(idxspm,1)=sListV_u(i,1)*sin(deltaAng);
	end;
      end;
    end;
  end;
end;
for iEtaSma=1:eta_v_sma
  for iXiSma=1:xi_v_sma
    if (MSKsma_v(iEtaSma, iXiSma) == 1)
      deltaAng=ANGsma_v(iEtaSma, iXiSma)-...
	       ANGsma_v_big(iEtaSma, iXiSma);
      ROMSdepthSma=Smacff_r+SmaCs_r*DEPsma_v(iEtaSma, iXiSma);
      for iNsma=1:Nsma
	idxSma_v=CorrespMatSma_v(iNsma, iEtaSma, iXiSma);
	MyDep=ROMSdepthSma(1,iNsma);
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
		idxbig1=CorrespMatBig_u(iNbig, iEtaBig, iXiBig);
		idxbig2=CorrespMatBig_u(iNbig+1, iEtaBig, iXiBig);
		nbU_v=nbU_v+1;
		iListU_v(nbU_v,1)=SIZ_u_sma+idxSma_v;
		jListU_v(nbU_v,1)=idxbig1;
		sListU_v(nbU_v,1)=eWeight*alpha1;
		nbU_v=nbU_v+1;
		iListU_v(nbU_v,1)=SIZ_u_sma+idxSma_v;
		jListU_v(nbU_v,1)=idxbig2;
		sListU_v(nbU_v,1)=eWeight*alpha2;
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		idxbig=CorrespMatBig_u(Nbig, iEtaBig, iXiBig);
		nbU_v=nbU_v+1;
		iListU_v(nbU_v,1)=SIZ_u_sma+idxSma_v;
		jListU_v(nbU_v,1)=idxbig;
		sListU_v(nbU_v,1)=eWeight;
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  for i=1:nbU_v
	    sListU_v(nbU_v,1)=sListU_v(nbU_v,1)/eSumWeight;
	  end;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_uv(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_uv(iEtaSma, iXiSma, idxdep);
	  idxbig=CorrespMatBig_u(1, iEtaBig, iXiBig);
	  nbU_v=nbU_v+1;
	  iListU_v(nbU_v,1)=SIZ_u_sma+idxSma_v;
	  jListU_v(nbU_v,1)=idxbig;
	  sListU_v(nbU_v,1)=1;
	end;
	%
	eSumWeight=0;
	nbV_v=0;
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
		idxbig1=CorrespMatBig_v(iNbig, iEtaBig, iXiBig);
		idxbig2=CorrespMatBig_v(iNbig+1, iEtaBig, iXiBig);
		nbV_v=nbV_v+1;
		iListV_v(nbV_v,1)=SIZ_u_sma+idxSma_v;
		jListV_v(nbV_v,1)=SIZ_u_big+idxbig1;
		sListV_v(nbV_v,1)=eWeight*alpha1;
		nbV_v=nbV_v+1;
		iListV_v(nbV_v,1)=SIZ_u_sma+idxSma_v;
		jListV_v(nbV_v,1)=SIZ_u_big+idxbig2;
		sListV_v(nbV_v,1)=eWeight*alpha2;
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		idxbig=CorrespMatBig_v(Nbig, iEtaBig, iXiBig);
		nbV_v=nbV_v+1;
		iListV_v(nbV_v,1)=SIZ_u_sma+idxSma_v;
		jListV_v(nbV_v,1)=SIZ_u_big+idxbig;
		sListV_v(nbV_v,1)=eWeight;
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  for i=1:nbV_v
	    sListV_v(i,1)=sListV_v(i,1)/eSumWeight;
	  end;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_vv(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_vv(iEtaSma, iXiSma, idxdep);
	  idxbig=CorrespMatBig_v(1, iEtaBig, iXiBig);
	  nbV_v=nbV_v+1;
	  iListV_v(nbV_v,1)=SIZ_u_sma+idxSma_v;
	  jListV_v(nbV_v,1)=SIZ_u_big+idxbig;
	  sListV_v(nbV_v,1)=1;
	end;
	for i=1:nbU_v
	  idxspm=idxspm+1;
	  iList(idxspm,1)=iListU_v(i,1);
	  jList(idxspm,1)=jListU_v(i,1);
	  sList(idxspm,1)=-sListU_v(i,1)*sin(deltaAng);
	end;
	for i=1:nbV_v
	  idxspm=idxspm+1;
	  iList(idxspm,1)=iListV_v(i,1);
	  jList(idxspm,1)=jListV_v(i,1);
	  sList(idxspm,1)=sListV_v(i,1)*cos(deltaAng);
	end;
      end;
    end;
  end;
end;
iListRed=iList(1:idxspm,1);
jListRed=jList(1:idxspm,1);
sListRed=sList(1:idxspm,1);
SPmat=sparse(iListRed, jListRed, sListRed, SIZ_sma, SIZ_big);
ArrayBigSma3Duv.SPmat=SPmat;
ArrayBigSma3Duv.eta_u_sma=eta_u_sma;
ArrayBigSma3Duv.xi_u_sma=xi_u_sma;
ArrayBigSma3Duv.eta_v_sma=eta_v_sma;
ArrayBigSma3Duv.xi_v_sma=xi_v_sma;
ArrayBigSma3Duv.SIZ_u_sma=SIZ_u_sma;
ArrayBigSma3Duv.SIZ_v_sma=SIZ_v_sma;
ArrayBigSma3Duv.SIZ_sma=SIZ_sma;
ArrayBigSma3Duv.Nsma=Nsma;
