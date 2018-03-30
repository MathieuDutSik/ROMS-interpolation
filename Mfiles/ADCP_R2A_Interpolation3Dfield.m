function TEMPadcp=ADCP_R2A_Interpolation3Dfield(TotalArray, TEMProms)
%
% TEMPadcp=ADCP_R2A_Interpolation3Dfield(TotalArray, TEMProms)
%
% TEMPadcp is the temperature as interpolated from the temperature
% of TEMProms at the specified depth of ADCPverticalLEvels.
%
% TotalArray is obtained by the function ADCP_R2A_GetTotalArray
% TEMProms is the sea surface height of the ROMS model.


AddiRecordVertDisc=TotalArray.AddiRecordVertDisc;
ADCPverticalLevels=TotalArray.ADCPverticalLevels;

Nbig=AddiRecordVertDisc.N;

nbADCP=TotalArray.nbADCP;
ListNvert=size(nbADCP, 1);
Nvert=size(ADCPverticalLevels, 2);

[Sc_w, Cs_w, Sc_r, Cs_r]=GRID_GetSc_Cs(...
    AddiRecordVertDisc.N, AddiRecordVertDisc.ThetaS, AddiRecordVertDisc.ThetaB);
cff_r=AddiRecordVertDisc.hc*(Sc_r-Cs_r);

ListRelETA_rho=TotalArray.ListRelETA_rho;
ListRelXI_rho=TotalArray.ListRelXI_rho;
ListRelCoeff_rho=TotalArray.ListRelCoeff_rho;

TEMPadcp=zeros(nbADCP, Nvert);
ListDepth=zeros(4, 1);
for iADCP=1:nbADCP
  eSum=0;
  eSumWeight=0;
  for iVert=1:Nvert
    MyDep=ADCPverticalLevels(iADCP, iVert);
    if (isnan(MyDep) == 0)
      eSum=0;
      eSumWeight=0;
      for idx=1:4
	ListDepth(idx,1)=0;
	eWeight=ListRelCoeff_rho(iADCP, 1, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_rho(iADCP, 1, idx);
	  iXiBig=ListRelXI_rho(iADCP, 1, idx);
	  TheDep=TotalArray.DEPbig_rho(iEtaBig, iXiBig);
	  ListDepth(idx,1)=TheDep;
	  ROMSdepthBig=cff_r+Cs_r*TheDep;
	  test=0;
	  for iNbig=1:Nbig-1
	    dep1=ROMSdepthBig(1,iNbig);
	    dep2=ROMSdepthBig(1,iNbig+1);
	    if (dep1 <= MyDep && MyDep <= dep2)
	      test=1;
	      alpha1=(dep2-MyDep)/(dep2-dep1);
	      alpha2=(MyDep-dep1)/(dep2-dep1);
	      eSum=eSum+eWeight*alpha1*...
		   TEMProms(iNbig, iEtaBig, iXiBig);
	      eSum=eSum+eWeight*alpha2*...
		   TEMProms(iNbig+1, iEtaBig, iXiBig);
	    end;
	  end;
	  if (test == 0)
	    if (ROMSdepthBig(1,Nbig) <= MyDep)
	      test=1;
	      eSum=eSum+eWeight*TEMProms(Nbig, iEtaBig, iXiBig);
	    end;
	  end;
	  if (test == 1)
	    eSumWeight=eSumWeight+eWeight;
	  end;
	end;
      end;
      if (eSumWeight > 0)
	TEMPadcp(iADCP, iVert)=eSum/eSumWeight;
      else
	Kdep=find(ListDepth == max(ListDepth));
	idxdep=Kdep(1,1);
	iEtaBig=ListRelETA_rho(iADCP, 1, idxdep);
	iXiBig=ListRelXI_rho(iADCP, 1, idxdep);
	TEMPadcp(iADCP, iVert)=TEMProms(1, iEtaBig, iXiBig);
      end;
    else
      TEMPadcp(iADCP, iVert)=NaN;
    end;
  end;
end;
