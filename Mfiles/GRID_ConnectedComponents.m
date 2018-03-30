function [ListConn, nbConn]=GRID_ConnectedComponents(MSK_rho)

[eta_rho, xi_rho]=size(MSK_rho);

ETAmat=zeros(eta_rho, xi_rho);
XImat=zeros(eta_rho, xi_rho);
for iEta=1:eta_rho
  for iXi=1:xi_rho
    ETAmat(iEta, iXi)=iEta;
    XImat(iEta, iXi)=iXi;
  end;
end;
ETAmatExp=ETAmat(:);
XImatExp=XImat(:);
MSK_rhoExp=MSK_rho(:);
siz=size(ETAmatExp,1);
CorrespMat=zeros(eta_rho, xi_rho);
idx=0;
for i=1:siz
  if (MSK_rhoExp(i, 1) == 1)
    idx=idx+1;
    iEta=ETAmatExp(i,1);
    iXi=XImatExp(i,1);
    CorrespMat(iEta, iXi)=idx;
  end;
end;
nbVert=idx;

UpperBoundNbEdges=eta_rho*xi_rho*2;
ListEdges=zeros(UpperBoundNbEdges, 2);
nbEdge=0;
for iEta1=2:eta_rho
  for iXi1=2:xi_rho
    for i=1:2
      if (i == 1)
	iEta2=iEta1-1;
	iXi2=iXi1;
      else
	iEta2=iEta1;
	iXi2=iXi1-1;
      end;
      if (MSK_rho(iEta1, iXi1) == 1 && MSK_rho(iEta2, iXi2) == 1)
	idx1=CorrespMat(iEta1, iXi1);
	idx2=CorrespMat(iEta2, iXi2);
	nbEdge=nbEdge+1;
	ListEdges(nbEdge,1)=idx1;
	ListEdges(nbEdge,2)=idx2;
      end;
    end;
  end;
end;
ListEdgesRed=ListEdges(1:nbEdge,1:2);
ListVertexStatus=GRAPH_ConnectedComponent(...
    ListEdgesRed, nbVert);

ListConn=zeros(eta_rho, xi_rho);
idx=0;
for i=1:siz
  if (MSK_rhoExp(i, 1) == 1)
    iEta=ETAmatExp(i,1);
    iXi=XImatExp(i,1);
    idx=idx+1;
    ListConn(iEta, iXi)=ListVertexStatus(idx, 1);
  end;
end;
nbConn=max(ListVertexStatus);
