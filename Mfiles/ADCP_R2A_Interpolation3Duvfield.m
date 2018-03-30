function [UV_east, UV_north]=...
    ADCP_R2A_Interpolation3Duvfield(TotalArray, Uroms, Vroms)
%
% [UV_east, UV_north]=...
%    ADCP_R2A_Interpolation3Duvfield(TotalArray, Uroms, Vroms)
%
% UV_east is the eastern component of velocity
% UV_north is the northern component of velocity
% 
% TotalArray is obtained by the function ADCP_R2A_GetTotalArray
% Uroms is the u-component of velocity
% Vroms is the v-component of velocity

AddiRecordVertDisc=TotalArray.AddiRecordVertDisc;
ADCPverticalLevels=TotalArray.ADCPverticalLevels;

Nbig=AddiRecordVertDisc.N;

nbADCP=TotalArray.nbADCP;
ListNvert=size(nbADCP, 1);
Nvert=size(ADCPverticalLevels, 2);

[Sc_w, Cs_w, Sc_r, Cs_r]=GRID_GetSc_Cs(...
    AddiRecordVertDisc.N, AddiRecordVertDisc.ThetaS, AddiRecordVertDisc.ThetaB);
cff_r=AddiRecordVertDisc.hc*(Sc_r-Cs_r);

ListRelETA_u=TotalArray.ListRelETA_u;
ListRelXI_u=TotalArray.ListRelXI_u;
ListRelCoeff_u=TotalArray.ListRelCoeff_u;

ListRelETA_v=TotalArray.ListRelETA_v;
ListRelXI_v=TotalArray.ListRelXI_v;
ListRelCoeff_v=TotalArray.ListRelCoeff_v;

DEPbig_u=TotalArray.DEPbig_u;
DEPbig_v=TotalArray.DEPbig_v;

UV_east=zeros(nbADCP, Nvert);
UV_north=zeros(nbADCP, Nvert);
for iADCP=1:nbADCP
  TheAng=TotalArray.ListAng(iADCP, 1);
  for iVert=1:Nvert
    MyDep=ADCPverticalLevels(iADCP, iVert);
    if (isnan(MyDep) == 0)
      eSumU_rho=0;
      eSumWeight=0;
      for idx=1:4
	ListDepth(idx,1)=0;
	eWeight=ListRelCoeff_u(iADCP, 1, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_u(iADCP, 1, idx);
	  iXiBig=ListRelXI_u(iADCP, 1, idx);
	  TheDep=DEPbig_u(iEtaBig, iXiBig);
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
	      eSumU_rho=eSumU_rho+eWeight*alpha1*...
			Uroms(iNbig, iEtaBig, iXiBig);
	      eSumU_rho=eSumU_rho+eWeight*alpha2*...
			Uroms(iNbig+1, iEtaBig, iXiBig);
	    end;
	  end;
	  if (test == 0)
	    if (ROMSdepthBig(1,Nbig) <= MyDep)
	      test=1;
	      eSumU_rho=eSumU_rho+eWeight*Uroms(Nbig, iEtaBig, iXiBig);
	    end;
	  end;
	  if (test == 1)
	    eSumWeight=eSumWeight+eWeight;
	  end;
	end;
      end;
      if (eSumWeight > 0)
	U_rho=eSumU_rho/eSumWeight;
      else
	Kdep=find(ListDepth == max(ListDepth));
	idxdep=Kdep(1,1);
	iEtaBig=ListRelETA_u(iADCP, 1, idxdep);
	iXiBig=ListRelXI_u(iADCP, 1, idxdep);
	U_rho=Uroms(1, iEtaBig, iXiBig);
      end;
      %
      eSumV_rho=0;
      eSumWeight=0;
      for idx=1:4
	ListDepth(idx,1)=0;
	eWeight=ListRelCoeff_v(iADCP, 1, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_v(iADCP, 1, idx);
	  iXiBig=ListRelXI_v(iADCP, 1, idx);
	  TheDep=DEPbig_v(iEtaBig, iXiBig);
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
	      eSumV_rho=eSumV_rho+eWeight*alpha1*...
			Vroms(iNbig, iEtaBig, iXiBig);
	      eSumV_rho=eSumV_rho+eWeight*alpha2*...
			Vroms(iNbig+1, iEtaBig, iXiBig);
	    end;
	  end;
	  if (test == 0)
	    if (ROMSdepthBig(1,Nbig) <= MyDep)
	      test=1;
	      eSumV_rho=eSumV_rho+eWeight*Vroms(Nbig, iEtaBig, iXiBig);
	    end;
	  end;
	  if (test == 1)
	    eSumWeight=eSumWeight+eWeight;
	  end;
	end;
      end;
      if (eSumWeight > 0)
	V_rho=eSumV_rho/eSumWeight;
      else
	Kdep=find(ListDepth == max(ListDepth));
	idxdep=Kdep(1,1);
	iEtaBig=ListRelETA_v(iADCP, 1, idxdep);
	iXiBig=ListRelXI_v(iADCP, 1, idxdep);
	V_rho=Vroms(1, iEtaBig, iXiBig);
      end;
      UVsingle_east=cos(TheAng)*U_rho-sin(TheAng)*V_rho;
      UVsingle_north=sin(TheAng)*U_rho+cos(TheAng)*V_rho;
      UV_east(iADCP, iVert)=UVsingle_east;
      UV_north(iADCP, iVert)=UVsingle_north;
    else
      UV_east(iADCP, iVert)=NaN;
      UV_north(iADCP, iVert)=NaN;
    end;
  end;
end;
