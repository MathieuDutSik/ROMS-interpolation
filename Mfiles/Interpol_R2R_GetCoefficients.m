function [ListETA, ListXI, ListCoeff]=Interpol_R2R_GetCoefficients(...
    LON_rho, LAT_rho, ...
    ListLON, ListLAT)

[siz1, siz2]=size(ListLON);
[eta_rho,xi_rho]=size(LON_rho);
[ETAmat, XImat]=Interpol_GetETAXImat(...
    eta_rho, xi_rho);
ETAmat_expand=ETAmat(:);
XImat_expand=XImat(:);
ListMatch=zeros(eta_rho, xi_rho);
nbNodeTot=size(ETAmat_expand,1);
ListNode=zeros(nbNodeTot, 2);
for iNode=1:nbNodeTot
  iEta=ETAmat_expand(iNode,1);
  iXi=XImat_expand(iNode,1);
  ListMatch(iEta, iXi)=iNode;
  ListNode(iNode,1)=LON_rho(iEta, iXi);
  ListNode(iNode,2)=LAT_rho(iEta, iXi);
end;

nbTrig=2*(eta_rho-1)*(xi_rho-1);
ListTrig=zeros(nbTrig,3);
IndexETA=zeros(nbTrig,1);
IndexXI=zeros(nbTrig,1);
ListNature=zeros(nbTrig,1);
iTrig=0;
for iEta=1:eta_rho-1
  for iXi=1:xi_rho-1;
    eCoord11=ListMatch(iEta, iXi);
    eCoord21=ListMatch(iEta+1, iXi);
    eCoord12=ListMatch(iEta, iXi+1);
    eCoord22=ListMatch(iEta+1, iXi+1);
    %
    iTrig=iTrig+1;
    ListTrig(iTrig,1)=eCoord11;
    ListTrig(iTrig,2)=eCoord21;
    ListTrig(iTrig,3)=eCoord12;
    IndexETA(iTrig,1)=iEta;
    IndexXI(iTrig,1)=iXi;
    ListNature(iTrig, 1)=1;
    %
    iTrig=iTrig+1;
    ListTrig(iTrig,1)=eCoord22;
    ListTrig(iTrig,2)=eCoord21;
    ListTrig(iTrig,3)=eCoord12;
    IndexETA(iTrig,1)=iEta;
    IndexXI(iTrig,1)=iXi;
    ListNature(iTrig, 1)=2;
  end;
end;
%ListTrig=delaunay(LON_rho_expand, LAT_rho_expand);
SIZ=size(ListLON);
ListLON_expand=ListLON(:);
ListLAT_expand=ListLAT(:);
ListPT=[ListLON_expand, ListLAT_expand];
ListTrigNrExpand=TRIG_Findelems(ListTrig, ListNode, ListPT);
ListTrigNr=reshape(ListTrigNrExpand, SIZ);

siz=size(ListLON, 1);
ListETA=zeros(siz1, siz2);
ListXI=zeros(siz1, siz2);
ListCoeff=zeros(siz1, siz2, 4);
isparse=0;
A1=zeros(3,3);
A2=zeros(3,3);
for i1=1:siz1
  for i2=1:siz2
    iTrig=ListTrigNr(i1, i2);
    if (iTrig > -1)
      iEta=IndexETA(iTrig,1);
      iXi=IndexXI(iTrig,1);
      iNature=ListNature(iTrig,1);
      B=[1 ListLON(i1, i2)  ListLAT(i1, i2)];
      if (iNature == 1)
	A1(1,:)=[1 LON_rho(iEta, iXi) LAT_rho(iEta,iXi)];
	A1(2,:)=[1 LON_rho(iEta+1, iXi) LAT_rho(iEta+1,iXi)];
	A1(3,:)=[1 LON_rho(iEta, iXi+1) LAT_rho(iEta,iXi+1)];
	eSol1=B*inv(A1);
	lambda1=eSol1(1,2);
	lambda2=eSol1(1,3);
	eCoeff00=(1-lambda1)*(1-lambda2);
	eCoeff10=lambda1*(1-lambda2);
	eCoeff01=(1-lambda1)*lambda2;
	eCoeff11=lambda1*lambda2;
      else
	A2(1,:)=[1 LON_rho(iEta+1, iXi+1) LAT_rho(iEta+1,iXi+1)];
	A2(2,:)=[1 LON_rho(iEta+1, iXi) LAT_rho(iEta+1,iXi)];
	A2(3,:)=[1 LON_rho(iEta, iXi+1) LAT_rho(iEta,iXi+1)];
	eSol2=B*inv(A2);
	lambda1=eSol2(1,3);
	lambda2=eSol2(1,2);
	eCoeff00=lambda1*lambda2;
	eCoeff10=(1-lambda1)*lambda2;
	eCoeff01=lambda1*(1-lambda2);
	eCoeff11=(1-lambda1)*(1-lambda2);
      end;
      if (eCoeff00 < -0.05 || eCoeff00 > 1.05 ...
	  || eCoeff10 < -0.05 || eCoeff10 > 1.05 ...
	  || eCoeff01 < -0.05 || eCoeff01 > 1.05 ...
	  || eCoeff11 < -0.05 || eCoeff11 > 1.05)
	disp(['Please debug from here in Interpol_R2R_GetCoefficients.m']);
	keyboard;
      end;
    else
      iEta=-1;
      iXi=-1;
      eCoeff00=0;
      eCoeff01=0;
      eCoeff10=0;
      eCoeff11=0;
    end;
    ListETA(i1,i2)=iEta;
    ListXI(i1,i2)=iXi;
    ListCoeff(i1, i2, 1)=eCoeff00;
    ListCoeff(i1, i2, 2)=eCoeff10;
    ListCoeff(i1, i2, 3)=eCoeff01;
    ListCoeff(i1, i2, 4)=eCoeff11;
  end;
end;
