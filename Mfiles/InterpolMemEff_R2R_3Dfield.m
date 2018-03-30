function TEMPsma=InterpolMemEff_R2R_3Dfield(...
    TotalArray, AddiRecordInfo, TEMPbig)
MSKbig_rho=TotalArray.MSKbig_rho;
MSKsma_rho=TotalArray.MSKsma_rho;
DEPbig_rho=TotalArray.DEPbig_rho;
DEPsma_rho=TotalArray.DEPsma_rho;
ListRelETA_rr=TotalArray.ListRelETA_rr;
ListRelXI_rr=TotalArray.ListRelXI_rr;
ListRelCoeff_rr=TotalArray.ListRelCoeff_rr;
%
[eta_rho_big, xi_rho_big]=size(MSKbig_rho);
[eta_rho_sma, xi_rho_sma]=size(MSKsma_rho);
%
Nbig=AddiRecordInfo.Nbig;
Nsma=AddiRecordInfo.Nsma;
%
[Nbig1, eta_rho_big1, xi_rho_big1]=size(TEMPbig);
if (eta_rho_big1 ~= eta_rho_big || xi_rho_big1 ~= xi_rho_big || ...
    Nbig1 ~= Nbig)
  disp('We have a size error between TotalArray and TEMPbig');
  keyboard;
end;

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
TEMPsma=zeros(Nsma, eta_rho_sma, xi_rho_sma);
for iEtaSma=1:eta_rho_sma
  for iXiSma=1:xi_rho_sma
    if (MSKsma_rho(iEtaSma, iXiSma) == 1)
      ROMSdepthSma=Smacff_r+SmaCs_r*DEPsma_rho(iEtaSma, iXiSma);
      for iNsma=1:Nsma
	eSum=0;
	eSumWeight=0;
	MyDep=ROMSdepthSma(1,iNsma);
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
		eSum=eSum+eWeight*alpha1*...
		     TEMPbig(iNbig, iEtaBig, iXiBig);
		eSum=eSum+eWeight*alpha2*...
		     TEMPbig(iNbig+1, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 0)
	      if (ROMSdepthBig(1,Nbig) <= MyDep)
		test=1;
		eSum=eSum+eWeight*TEMPbig(Nbig, iEtaBig, iXiBig);
	      end;
	    end;
	    if (test == 1)
	      eSumWeight=eSumWeight+eWeight;
	    end;
	  end;
	end;
	if (eSumWeight > 0)
	  TEMPsma(iNsma, iEtaSma,iXiSma)=eSum/eSumWeight;
	else
	  Kdep=find(ListDepth == max(ListDepth));
	  idxdep=Kdep(1,1);
	  iEtaBig=ListRelETA_rr(iEtaSma, iXiSma, idxdep);
	  iXiBig=ListRelXI_rr(iEtaSma, iXiSma, idxdep);
	  TEMPsma(iNsma, iEtaSma,iXiSma)=TEMPbig(1, iEtaBig, iXiBig);
	end;
      end;
    end;
  end;
end;
