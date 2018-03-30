function ArrayBigSma3D=InterpolGetSpMat_R2R_3Dfield(...
    MSKbig_rho, MSKsma_rho, ...
    DEPbig_rho, DEPsma_rho, ...
    ListRelETA_rr, ListRelXI_rr, ListRelCoeff_rr, ...
    AddiRecordInfo)
disp('Computing sparse matrix for 2D + 1D interpolation');
%
[eta_rho_big, xi_rho_big]=size(MSKbig_rho);
[eta_rho_sma, xi_rho_sma]=size(MSKsma_rho);
%
Nbig=AddiRecordInfo.Nbig;
Nsma=AddiRecordInfo.Nsma;
CorrespMatBig=GetCorrespMatrix(3, [Nbig; eta_rho_big; xi_rho_big]);
CorrespMatSma=GetCorrespMatrix(3, [Nsma; eta_rho_sma; xi_rho_sma]);
SIZ_sma=eta_rho_sma*xi_rho_sma*Nsma;
SIZ_big=eta_rho_big*xi_rho_big*Nbig;
%
[BigSc_w, BigCs_w, BigSc_r, BigCs_r]=GRID_GetSc_Cs(...
    AddiRecordInfo.Nbig, AddiRecordInfo.BigThetaS, AddiRecordInfo.BigThetaB);
Bigcff_r=AddiRecordInfo.Bighc*(BigSc_r-BigCs_r);
[SmaSc_w, SmaCs_w, SmaSc_r, SmaCs_r]=GRID_GetSc_Cs(...
    AddiRecordInfo.Nsma, AddiRecordInfo.SmaThetaS, AddiRecordInfo.SmaThetaB);
Smacff_r=AddiRecordInfo.Smahc*(SmaSc_r-SmaCs_r);
%
ListDepth=zeros(4,1);
%
UpperEstPoint=eta_rho_sma*xi_rho_sma*8*Nsma;
iList=zeros(UpperEstPoint, 1);
jList=zeros(UpperEstPoint, 1);
sList=zeros(UpperEstPoint, 1);
iListPar=zeros(8,1);
jListPar=zeros(8,1);
sListPar=zeros(8,1);
idxspm=0;
for iEtaSma=1:eta_rho_sma
  for iXiSma=1:xi_rho_sma
    if (MSKsma_rho(iEtaSma, iXiSma) == 1)
      ROMSdepthSma=Smacff_r+SmaCs_r*DEPsma_rho(iEtaSma, iXiSma);
      for iNsma=1:Nsma
	idxSma=CorrespMatSma(iNsma, iEtaSma,iXiSma);
	eSumWeight=0;
	MyDep=ROMSdepthSma(1,iNsma);
	nbPar=0;
	for idx=1:4
	  ListDepth(idx,1)=0;
	  eWeight=ListRelCoeff_rr(iEtaSma, iXiSma, idx);
	  if (eWeight > 0)
	    iEtaBig=ListRelETA_rr(iEtaSma, iXiSma, idx);
	    iXiBig=ListRelXI_rr(iEtaSma, iXiSma, idx);
	    depbig=DEPbig_rho(iEtaBig, iXiBig);
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
		idxbig1=CorrespMatBig(iNbig, iEtaBig, iXiBig);
		idxbig2=CorrespMatBig(iNbig+1, iEtaBig, iXiBig);
		nbPar=nbPar+1;
		iListPar(nbPar,1)=idxSma;
		jListPar(nbPar,1)=idxbig1;
		sListPar(nbPar,1)=eWeight*alpha1;
		nbPar=nbPar+1;
		iListPar(nbPar,1)=idxSma;
		jListPar(nbPar,1)=idxbig2;
		sListPar(nbPar,1)=eWeight*alpha2;
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		idxbig=CorrespMatBig(Nbig, iEtaBig, iXiBig);
		nbPar=nbPar+1;
		iListPar(nbPar,1)=idxSma;
		jListPar(nbPar,1)=idxbig;
		sListPar(nbPar,1)=eWeight;
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  for i=1:nbPar
	    idxspm=idxspm+1;
	    iList(idxspm,1)=iListPar(i,1);
	    jList(idxspm,1)=jListPar(i,1);
	    sList(idxspm,1)=sListPar(i,1)/eSumWeight;
	  end;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_rr(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_rr(iEtaSma, iXiSma, idxdep);
	  idxbig=CorrespMatBig(1, iEtaBig, iXiBig);
	  idxspm=idxspm+1;
	  iList(idxspm,1)=idxSma;
	  jList(idxspm,1)=idxbig;
	  sList(idxspm,1)=1;
	end;
      end;
    end;
  end;
end;
iListRed=iList(1:idxspm,1);
jListRed=jList(1:idxspm,1);
sListRed=sList(1:idxspm,1);
SPmat=sparse(iListRed, jListRed, sListRed, SIZ_sma, SIZ_big);
ArrayBigSma3D.SPmat=SPmat;
ArrayBigSma3D.eta_rho_sma=eta_rho_sma;
ArrayBigSma3D.xi_rho_sma=xi_rho_sma;
ArrayBigSma3D.Nsma=Nsma;
ArrayBigSma3D.Nbig=Nbig;
