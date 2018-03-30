function ArrayBigSma2D=InterpolGetSpMat_R2R_2Dfield(TotalArray)
MSKbig_rho=TotalArray.MSKbig_rho;
MSKsma_rho=TotalArray.MSKsma_rho;
ListRelETA_rr=TotalArray.ListRelETA_rr;
ListRelXI_rr=TotalArray.ListRelXI_rr;
ListRelCoeff_rr=TotalArray.ListRelCoeff_rr;
%
disp('Computing sparse matrix for 2D interpolation');
%
[eta_rho_big, xi_rho_big]=size(MSKbig_rho);
SIZ_big=eta_rho_big*xi_rho_big;
%
[eta_rho_sma, xi_rho_sma]=size(MSKsma_rho);
SIZ_sma=eta_rho_sma*xi_rho_sma;
%
CorrespMatBig=GetCorrespMatrix(2, [eta_rho_big; xi_rho_big]);
CorrespMatSma=GetCorrespMatrix(2, [eta_rho_sma; xi_rho_sma]);
%
ListOccurBig_rho=zeros(eta_rho_big, xi_rho_big);
%
KseaRhoSma=find(MSKsma_rho == 1);
nbWetRhoSma=size(KseaRhoSma, 1);
UpperEstPoint=nbWetRhoSma*4;
iList=zeros(UpperEstPoint, 1);
jList=zeros(UpperEstPoint, 1);
sList=zeros(UpperEstPoint, 1);
idxspm=0;
iListPar=zeros(4,1);
jListPar=zeros(4,1);
sListPar=zeros(4,1);
for iEtaSma=1:eta_rho_sma
  for iXiSma=1:xi_rho_sma
    if (MSKsma_rho(iEtaSma, iXiSma) == 1)
      idxSma=CorrespMatSma(iEtaSma,iXiSma);
      nbPar=0;
      for idx=1:4
	eWeight=ListRelCoeff_rr(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_rr(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_rr(iEtaSma, iXiSma, idx);
	  idxBig=CorrespMatBig(iEtaBig,iXiBig);
	  ListOccurBig_rho(iEtaBig, iXiBig)=1;
	  nbPar=nbPar+1;
	  iListPar(nbPar,1)=idxSma;
	  jListPar(nbPar,1)=idxBig;
	  sListPar(nbPar,1)=eWeight;
	end;
      end;
      for i=1:nbPar
	idxspm=idxspm+1;
	iList(idxspm,1)=iListPar(i,1);
	jList(idxspm,1)=jListPar(i,1);
	sList(idxspm,1)=sListPar(i,1);
      end;
    end;
  end;
end;
iListRed=iList(1:idxspm,1);
jListRed=jList(1:idxspm,1);
sListRed=sList(1:idxspm,1);
SPmat=sparse(iListRed, jListRed, sListRed, SIZ_sma, SIZ_big);
%
K=find(ListOccurBig_rho == 1);
[ETAmat, XImat]=Interpol_GetETAXImat(eta_rho_big, xi_rho_big);
ListRhoPoints_eta=ETAmat(K);
ListRhoPoints_xi=XImat(K);
%
ArrayBigSma2D.SPmat=SPmat;
ArrayBigSma2D.eta_rho_sma=eta_rho_sma;
ArrayBigSma2D.xi_rho_sma=xi_rho_sma;
ArrayBigSma2D.SIZ_sma=SIZ_sma;
ArrayBigSma2D.ListOccurBig_rho=ListOccurBig_rho;
ArrayBigSma2D.ListRhoPoints_eta=ListRhoPoints_eta;
ArrayBigSma2D.ListRhoPoints_xi=ListRhoPoints_xi;
ArrayBigSma2D.eta_rho_big=eta_rho_big;
ArrayBigSma2D.xi_rho_big=xi_rho_big;
