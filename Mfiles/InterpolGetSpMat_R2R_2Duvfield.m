function ArrayBigSma2Duv=InterpolGetSpMat_R2R_2Duvfield(TotalArray)
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
disp('Computing sparse matrix for 2D UV-field interpolation');
%
[eta_rho_big,xi_rho_big]=size(MSKbig_rho);
[eta_u_big,xi_u_big]=size(MSKbig_u);
[eta_v_big,xi_v_big]=size(MSKbig_v);
ANGbig_u=(ANGbig_rho(:, 1:xi_rho_big-1)+ANGbig_rho(:, 2:xi_rho_big))/2;
ANGbig_v=(ANGbig_rho(1:eta_rho_big-1, :)+ANGbig_rho(2:eta_rho_big, :))/2;
%
[eta_rho_sma,xi_rho_sma]=size(MSKsma_rho);
[eta_u_sma,xi_u_sma]=size(MSKsma_u);
[eta_v_sma,xi_v_sma]=size(MSKsma_v);
ANGsma_u=(ANGsma_rho(:, 1:xi_rho_sma-1)+ANGsma_rho(:, 2:xi_rho_sma))/2;
ANGsma_v=(ANGsma_rho(1:eta_rho_sma-1, :)+ANGsma_rho(2:eta_rho_sma, :))/2;
%
CorrespMatBig_u=GetCorrespMatrix(2, [eta_u_big; xi_u_big]);
CorrespMatBig_v=GetCorrespMatrix(2, [eta_v_big; xi_v_big]);
CorrespMatSma_u=GetCorrespMatrix(2, [eta_u_sma; xi_u_sma]);
CorrespMatSma_v=GetCorrespMatrix(2, [eta_v_sma; xi_v_sma]);
ListOccurBig_u=zeros(eta_u_big, xi_u_big);
ListOccurBig_v=zeros(eta_v_big, xi_v_big);
%
SIZ_u_big=eta_u_big*xi_u_big;
SIZ_v_big=eta_v_big*xi_v_big;
SIZ_u_sma=eta_u_sma*xi_u_sma;
SIZ_v_sma=eta_v_sma*xi_v_sma;
%
SIZ_sma=SIZ_u_sma+SIZ_v_sma;
SIZ_big=SIZ_u_big+SIZ_v_big;
%
KseaUsma=find(MSKsma_u == 1);
nbWetUsma=size(KseaUsma,1);
KseaVsma=find(MSKsma_v == 1);
nbWetVsma=size(KseaVsma,1);
UpperEstPoint=nbWetUsma*8+nbWetVsma*8;
iList=zeros(UpperEstPoint, 1);
jList=zeros(UpperEstPoint, 1);
sList=zeros(UpperEstPoint, 1);
iListUBAR_u=zeros(4,1);
jListUBAR_u=zeros(4,1);
sListUBAR_u=zeros(4,1);
iListVBAR_u=zeros(4,1);
jListVBAR_u=zeros(4,1);
sListVBAR_u=zeros(4,1);
idxspm=0;
for iEtaSma=1:eta_u_sma
  for iXiSma=1:xi_u_sma
    if (MSKsma_u(iEtaSma, iXiSma) == 1)
      deltaAng=ANGsma_u(iEtaSma, iXiSma)-...
	  ANGsma_u_big(iEtaSma, iXiSma);
      idxSma_u=CorrespMatSma_u(iEtaSma,iXiSma);
      nbUBAR_u=0;
      for idx=1:4
	eWeight=ListRelCoeff_uu(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_uu(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_uu(iEtaSma, iXiSma, idx);
	  ListOccurBig_u(iEtaBig, iXiBig)=1;
	  idxBig_u=CorrespMatBig_u(iEtaBig, iXiBig);
	  nbUBAR_u=nbUBAR_u+1;
	  iListUBAR_u(nbUBAR_u,1)=idxSma_u;
	  jListUBAR_u(nbUBAR_u,1)=idxBig_u;
	  sListUBAR_u(nbUBAR_u,1)=eWeight;
	end;
      end;
      nbVBAR_u=0;
      for idx=1:4
	eWeight=ListRelCoeff_vu(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_vu(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_vu(iEtaSma, iXiSma, idx);
	  ListOccurBig_v(iEtaBig, iXiBig)=1;
	  idxBig_v=CorrespMatBig_v(iEtaBig,iXiBig);
	  nbVBAR_u=nbVBAR_u+1;
	  iListVBAR_u(nbVBAR_u,1)=idxSma_u;
	  jListVBAR_u(nbVBAR_u,1)=SIZ_u_big+idxBig_v;
	  sListVBAR_u(nbVBAR_u,1)=eWeight;
	end;
      end;
      for i=1:nbUBAR_u;
	idxspm=idxspm+1;
	iList(idxspm,1)=iListUBAR_u(i,1);
	jList(idxspm,1)=jListUBAR_u(i,1);
	sList(idxspm,1)=sListUBAR_u(i,1)*cos(deltaAng);
      end;
      for i=1:nbVBAR_u
	idxspm=idxspm+1;
	iList(idxspm,1)=iListVBAR_u(i,1);
	jList(idxspm,1)=jListVBAR_u(i,1);
	sList(idxspm,1)=sListVBAR_u(i,1)*sin(deltaAng);
      end;
    end;
  end;
end;
iListUBAR_v=zeros(4,1);
jListUBAR_v=zeros(4,1);
sListUBAR_v=zeros(4,1);
iListVBAR_v=zeros(4,1);
jListVBAR_v=zeros(4,1);
sListVBAR_v=zeros(4,1);
for iEtaSma=1:eta_v_sma
  for iXiSma=1:xi_v_sma
    if (MSKsma_v(iEtaSma, iXiSma) == 1)
      idxSma_v=CorrespMatSma_v(iEtaSma,iXiSma);
      deltaAng=ANGsma_v(iEtaSma, iXiSma)-...
	  ANGsma_v_big(iEtaSma, iXiSma);
      nbUBAR_v=0;
      for idx=1:4
	eWeight=ListRelCoeff_uv(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_uv(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_uv(iEtaSma, iXiSma, idx);
	  ListOccurBig_u(iEtaBig, iXiBig)=1;
	  idxBig_u=CorrespMatBig_u(iEtaBig,iXiBig);
	  nbUBAR_v=nbUBAR_v+1;
	  iListUBAR_v(nbUBAR_v,1)=SIZ_u_sma+idxSma_v;
	  jListUBAR_v(nbUBAR_v,1)=idxBig_u;
	  sListUBAR_v(nbUBAR_v,1)=eWeight;
	end;
      end;
      nbVBAR_v=0;
      for idx=1:4
	eWeight=ListRelCoeff_vv(iEtaSma, iXiSma, idx);
	if (eWeight > 0)
	  iEtaBig=ListRelETA_vv(iEtaSma, iXiSma, idx);
	  iXiBig=ListRelXI_vv(iEtaSma, iXiSma, idx);
	  ListOccurBig_v(iEtaBig, iXiBig)=1;
	  idxBig_v=CorrespMatBig_v(iEtaBig,iXiBig);
	  nbVBAR_v=nbVBAR_v+1;
	  iListVBAR_v(nbVBAR_v,1)=SIZ_u_sma+idxSma_v;
	  jListVBAR_v(nbVBAR_v,1)=SIZ_u_big+idxBig_v;
	  sListVBAR_v(nbVBAR_v,1)=eWeight;
	end;
      end;
      for i=1:nbUBAR_v
	idxspm=idxspm+1;
	iList(idxspm,1)=iListUBAR_v(i,1);
	jList(idxspm,1)=jListUBAR_v(i,1);
	sList(idxspm,1)=-sListUBAR_v(i,1)*sin(deltaAng);
      end;
      for i=1:nbVBAR_v
	idxspm=idxspm+1;
	iList(idxspm,1)=iListVBAR_v(i,1);
	jList(idxspm,1)=jListVBAR_v(i,1);
	sList(idxspm,1)=sListVBAR_v(i,1)*cos(deltaAng);
      end;
    end;
  end;
end;
iListRed=iList(1:idxspm,1);
jListRed=jList(1:idxspm,1);
sListRed=sList(1:idxspm,1);
SPmat=sparse(iListRed, jListRed, sListRed, SIZ_sma, SIZ_big);
%
K=find(ListOccurBig_u == 1);
[ETAmat, XImat]=Interpol_GetETAXImat(eta_u_big, xi_u_big);
ListUPoints_eta=ETAmat(K);
ListUPoints_xi=XImat(K);
K=find(ListOccurBig_v == 1);
[ETAmat, XImat]=Interpol_GetETAXImat(eta_v_big, xi_v_big);
ListVPoints_eta=ETAmat(K);
ListVPoints_xi=XImat(K);
%
ArrayBigSma2Duv.SPmat=SPmat;
ArrayBigSma2Duv.SIZ_sma=SIZ_sma;
ArrayBigSma2Duv.SIZ_u_sma=SIZ_u_sma;
ArrayBigSma2Duv.SIZ_v_sma=SIZ_v_sma;
ArrayBigSma2Duv.eta_u_sma=eta_u_sma;
ArrayBigSma2Duv.xi_u_sma=xi_u_sma;
ArrayBigSma2Duv.eta_v_sma=eta_v_sma;
ArrayBigSma2Duv.xi_v_sma=xi_v_sma;
ArrayBigSma2Duv.ListOccurBig_u=ListOccurBig_u;
ArrayBigSma2Duv.ListOccurBig_v=ListOccurBig_v;
ArrayBigSma2Duv.ListUPoints_eta=ListUPoints_eta;
ArrayBigSma2Duv.ListUPoints_xi=ListUPoints_xi;
ArrayBigSma2Duv.ListVPoints_eta=ListVPoints_eta;
ArrayBigSma2Duv.ListVPoints_xi=ListVPoints_xi;
ArrayBigSma2Duv.SIZ_big=SIZ_big;
ArrayBigSma2Duv.SIZ_u_big=SIZ_u_big;
ArrayBigSma2Duv.SIZ_v_big=SIZ_v_big;
ArrayBigSma2Duv.eta_u_big=eta_u_big;
ArrayBigSma2Duv.xi_u_big=xi_u_big;
ArrayBigSma2Duv.eta_v_big=eta_v_big;
ArrayBigSma2Duv.xi_v_big=xi_v_big;
