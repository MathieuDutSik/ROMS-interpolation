function TheArray=GRID_GetArray(GridFile)
nc=netcdf(GridFile, 'nowrite');
LON_rho=nc{'lon_rho'}(:);
LAT_rho=nc{'lat_rho'}(:);
LON_psi=nc{'lon_psi'}(:);
LAT_psi=nc{'lat_psi'}(:);
LON_u=nc{'lon_u'}(:);
LAT_u=nc{'lat_u'}(:);
LON_v=nc{'lon_v'}(:);
LAT_v=nc{'lat_v'}(:);
MSK_rho=nc{'mask_rho'}(:);
DEP_rho=nc{'h'}(:);
ANG_rho=nc{'angle'}(:);
pm_rho=nc{'pm'}(:);
pn_rho=nc{'pn'}(:);
fCoriolis=nc{'f'}(:);
dndx=nc{'dndx'}(:);
dmde=nc{'dmde'}(:);
close(nc);
%
TheArray.name=GridFile;
TheArray.LON_rho=LON_rho;
TheArray.LAT_rho=LAT_rho;
TheArray.LON_psi=LON_psi;
TheArray.LAT_psi=LAT_psi;
TheArray.LON_u=LON_u;
TheArray.LAT_u=LAT_u;
TheArray.LON_v=LON_v;
TheArray.LAT_v=LAT_v;
TheArray.fCoriolis=fCoriolis;
%
TheArray.dndx=dndx;
TheArray.dmde=dmde;
%
[eta_rho, xi_rho]=size(LON_rho);
[eta_u, xi_u]=size(LON_u);
[eta_v, xi_v]=size(LON_v);
TheArray.eta_rho=eta_rho;
TheArray.eta_psi=eta_rho-1;
TheArray.eta_u=eta_u;
TheArray.eta_v=eta_v;
TheArray.xi_rho=xi_rho;
TheArray.xi_psi=xi_rho-1;
TheArray.xi_u=xi_u;
TheArray.xi_v=xi_v;
%
SIZ_rho=eta_rho*xi_rho;
SIZ_u=eta_u*xi_u;
SIZ_v=eta_v*xi_v;
TheArray.SIZ_rho=SIZ_rho;
TheArray.SIZ_u=SIZ_u;
TheArray.SIZ_v=SIZ_v;
%
MSK_v=MSK_rho(1:eta_rho-1,:).*MSK_rho(2:eta_rho,:);
MSK_u=MSK_rho(:,1:xi_rho-1).*MSK_rho(:,2:xi_rho);
MSK_psi=MSK_rho(1:eta_rho-1,1:xi_rho-1).*...
        MSK_rho(2:eta_rho,1:xi_rho-1).*...
        MSK_rho(1:eta_rho-1,2:xi_rho).*...
        MSK_rho(2:eta_rho,2:xi_rho);
TheArray.MSK_rho=MSK_rho;
TheArray.MSK_psi=MSK_psi;
TheArray.MSK_u=MSK_u;
TheArray.MSK_v=MSK_v;
MSKbound_rho=zeros(eta_rho, xi_rho);
MSKbound_rho(1, :)=1;
MSKbound_rho(:, 1)=1;
MSKbound_rho(eta_rho, :)=1;
MSKbound_rho(:, xi_rho)=1;
TheArray.MSKbound_rho=MSKbound_rho;
Kbound_rho=find(MSKbound_rho == 1);
TheArray.Kbound_rho=Kbound_rho;
%
RMat=GRID_RoughnessMatrix(DEP_rho, MSK_rho);
rx0=max(RMat(:));
TheArray.RMat=RMat;
TheArray.rx0=rx0;
%
ANG_u=(ANG_rho(:, 1:xi_u)+ANG_rho(:, 2:xi_rho))/2;
ANG_v=(ANG_rho(1:eta_v, :)+ANG_rho(2:eta_rho, :))/2;
ANG_psi=(ANG_rho(1:eta_rho-1,1:xi_rho-1)+...
           ANG_rho(2:eta_rho,1:xi_rho-1)+...
           ANG_rho(1:eta_rho-1,2:xi_rho)+...
           ANG_rho(2:eta_rho,2:xi_rho))/4;
TheArray.ANG_rho=ANG_rho;
TheArray.ANG_psi=ANG_psi;
TheArray.ANG_u=ANG_u;
TheArray.ANG_v=ANG_v;
%
SigSize=[num2str(eta_rho) 'x' num2str(xi_rho)];
TheArray.SigSize=SigSize;
%
mn_rho=pm_rho.*pn_rho;
mn_psi=(mn_rho(1:eta_rho-1,1:xi_rho-1)+mn_rho(1:eta_rho-1,2:xi_rho)+...
	mn_rho(2:eta_rho,2:xi_rho)+mn_rho(2:eta_rho,1:xi_rho-1))/4;
TheArray.pm_rho=pm_rho;
TheArray.pn_rho=pn_rho;
TheArray.pm=pm_rho;
TheArray.pn=pn_rho;
TheArray.mn_rho=mn_rho;
TheArray.mn_psi=mn_psi;
%
DEP_u=(DEP_rho(:, 1:xi_u)+DEP_rho(:, 2:xi_rho))/2;
DEP_v=(DEP_rho(1:eta_v, :)+DEP_rho(2:eta_rho, :))/2;
DEP_psi=(DEP_rho(1:eta_rho-1,1:xi_rho-1)+...
           DEP_rho(2:eta_rho,1:xi_rho-1)+...
           DEP_rho(1:eta_rho-1,2:xi_rho)+...
           DEP_rho(2:eta_rho,2:xi_rho))/4;
TheArray.DEP_rho=DEP_rho;
TheArray.DEP_psi=DEP_psi;
TheArray.DEP_u=DEP_u;
TheArray.DEP_v=DEP_v;
%
[ListConn, nbConn]=GRID_ConnectedComponents(MSK_rho);
TheArray.ListConn=ListConn;
TheArray.nbConn=nbConn;
%
KlandPsi=find(MSK_psi == 0);
KseaPsi=find(MSK_psi == 1);
KlandRho=find(MSK_rho == 0);
KseaRho=find(MSK_rho == 1);
KlandU=find(MSK_u == 0);
KseaU=find(MSK_u == 1);
KlandV=find(MSK_v == 0);
KseaV=find(MSK_v == 1);
TheArray.KlandPsi=KlandPsi;
TheArray.KseaPsi=KseaPsi;
TheArray.KlandRho=KlandRho;
TheArray.KseaRho=KseaRho;
TheArray.KlandU=KlandU;
TheArray.KseaU=KseaU;
TheArray.KlandV=KlandV;
TheArray.KseaV=KseaV;
%
LONnorth_u=LON_u(eta_u,:);
LONsouth_u=LON_u(1,:);
LONwest_u=LON_u(:,1);
LONeast_u=LON_u(:,xi_u);
LATnorth_u=LAT_u(eta_u,:);
LATsouth_u=LAT_u(1,:);
LATwest_u=LAT_u(:,1);
LATeast_u=LAT_u(:,xi_u);
ANGnorth_u=ANG_u(eta_u,:);
ANGsouth_u=ANG_u(1,:);
ANGwest_u=ANG_u(:,1);
ANGeast_u=ANG_u(:,xi_u);
LONnorth_v=LON_v(eta_v,:);
LONsouth_v=LON_v(1,:);
LONwest_v=LON_v(:,1);
LONeast_v=LON_v(:,xi_v);
LATnorth_v=LAT_v(eta_v,:);
LATsouth_v=LAT_v(1,:);
LATwest_v=LAT_v(:,1);
LATeast_v=LAT_v(:,xi_v);
ANGnorth_v=ANG_v(eta_v,:);
ANGsouth_v=ANG_v(1,:);
ANGwest_v=ANG_v(:,1);
ANGeast_v=ANG_v(:,xi_v);
%
LONbound=[LON_rho(1,1) LON_rho(1, xi_rho) ...
	  LON_rho(eta_rho,xi_rho) LON_rho(eta_rho,1) LON_rho(1,1)];
LATbound=[LAT_rho(1,1) LAT_rho(1, xi_rho) ...
	  LAT_rho(eta_rho,xi_rho) LAT_rho(eta_rho,1) LAT_rho(1,1)];
TheArray.LONbound=LONbound;
TheArray.LATbound=LATbound;
%
minDepth=min(DEP_rho(KseaRho));
TheArray.minDepth=minDepth;
%
nbWetRho=size(KseaRho,1);
nbWetPsi=size(KseaPsi,1);
nbWetU=size(KseaU,1);
nbWetV=size(KseaV,1);
TheArray.nbWetRho=nbWetRho;
TheArray.nbWetPsi=nbWetPsi;
TheArray.nbWetU=nbWetU;
TheArray.nbWetV=nbWetV;
TheArray.xy_rho=nbWetRho;
TheArray.xy_u=nbWetU;
TheArray.xy_v=nbWetV;
%
InflMatrixSph_rho=GRID_GetMaxRadiusInfluence_kernel(...
    LON_rho, LAT_rho, MSK_rho, 'spherical');
InflMatrixSph_u=GRID_GetMaxRadiusInfluence_kernel(...
    LON_u, LAT_u, MSK_u, 'spherical');
InflMatrixSph_v=GRID_GetMaxRadiusInfluence_kernel(...
    LON_v, LAT_v, MSK_v, 'spherical');
InflMatrixSph_psi=GRID_GetMaxRadiusInfluence_kernel(...
    LON_psi, LAT_psi, MSK_psi, 'spherical');
InflMatrixEuc_rho=GRID_GetMaxRadiusInfluence_kernel(...
    LON_rho, LAT_rho, MSK_rho, 'euclidean');
InflMatrixEuc_u=GRID_GetMaxRadiusInfluence_kernel(...
    LON_u, LAT_u, MSK_u, 'euclidean');
InflMatrixEuc_v=GRID_GetMaxRadiusInfluence_kernel(...
    LON_v, LAT_v, MSK_v, 'euclidean');
InflMatrixEuc_psi=GRID_GetMaxRadiusInfluence_kernel(...
    LON_psi, LAT_psi, MSK_psi, 'euclidean');
TheArray.InflMatrixSph_rho=InflMatrixSph_rho;
TheArray.InflMatrixSph_u=InflMatrixSph_u;
TheArray.InflMatrixSph_v=InflMatrixSph_v;
TheArray.InflMatrixSph_psi=InflMatrixSph_psi;
TheArray.InflMatrixEuc_rho=InflMatrixEuc_rho;
TheArray.InflMatrixEuc_u=InflMatrixEuc_u;
TheArray.InflMatrixEuc_v=InflMatrixEuc_v;
TheArray.InflMatrixEuc_psi=InflMatrixEuc_psi;
%
[LON_psi2, LAT_psi2, MSK_psi2]=GRID_ExtendedPsiSec(...
    MSK_rho, LON_rho, LAT_rho, LON_psi, LAT_psi, ...
    LON_u, LAT_u, LON_v, LAT_v);
TheArray.LON_psi2=LON_psi2;
TheArray.LAT_psi2=LAT_psi2;
TheArray.MSK_psi2=MSK_psi2;

%
[nbPositive, nbNegative]=GRID_GridOrientationInfo(LON_rho, LAT_rho);
if (nbPositive > 0 && nbNegative == 0)
  TheOrientation='positive';
elseif (nbPositive == 0 && nbNegative > 0)
  TheOrientation='negative';
else
  TheOrientation='unknown';
end;
TheArray.TheOrientation=TheOrientation;
%
result=GRID_SpatialDisc(MSK_rho, LON_rho, LAT_rho);
TheArray.TheOrientation=TheOrientation;
TheArray.MaxDist_ETA=result.MaxDist_ETA;
TheArray.MinDist_ETA=result.MinDist_ETA;
TheArray.AveDist_ETA=result.AveDist_ETA;
TheArray.MaxDist_XI=result.MaxDist_XI;
TheArray.MinDist_XI=result.MinDist_XI;
TheArray.AveDist_XI=result.AveDist_XI;
%
MinLon=min(LON_rho(KseaRho));
MaxLon=max(LON_rho(KseaRho));
MinLat=min(LAT_rho(KseaRho));
MaxLat=max(LAT_rho(KseaRho));
TheArray.MinLon=MinLon;
TheArray.MaxLon=MaxLon;
TheArray.MinLat=MinLat;
TheArray.MaxLat=MaxLat;
TheArray.TheQuad=[MinLon MaxLon, MinLat MaxLat];
%
EarthRadius=6367000;
xl=EarthRadius*GeodesicDistance(...
    [LON_rho(1,1), LAT_rho(1,1)], ...
    [LON_rho(1,xi_rho), LAT_rho(1,xi_rho)]);
el=EarthRadius*GeodesicDistance(...
    [LON_rho(1,1), LAT_rho(1,1)], ...
    [LON_rho(eta_rho,1), LAT_rho(eta_rho,1)]);
TheArray.xl=xl;
TheArray.el=el;
%
[str1, str2, str3, str4, ...
 str5, str6, str7, str8]=GRID_GetEWNSinfo(MSK_rho);
TheArray.str1=str1;
TheArray.str2=str2;
TheArray.str3=str3;
TheArray.str4=str4;
TheArray.str5=str5;
TheArray.str6=str6;
TheArray.str7=str7;
TheArray.str8=str8;
%
