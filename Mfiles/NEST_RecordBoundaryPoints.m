function TheRecord=NEST_RecordBoundaryPoints(gridNEST)

nc=netcdf(gridNEST);
lon_rho=nc{'lon_rho'}(:);
lat_rho=nc{'lat_rho'}(:);
mask_rho=nc{'mask_rho'}(:);
eta_rho=size(lon_rho,1);
xi_rho=size(lon_rho,2);

lon_u=nc{'lon_u'}(:);
lat_u=nc{'lat_u'}(:);
mask_u=nc{'mask_u'}(:);
eta_u=size(lon_u,1);
xi_u=size(lon_u,2);

lon_v=nc{'lon_v'}(:);
lat_v=nc{'lat_v'}(:);
mask_v=nc{'mask_v'}(:);
eta_v=size(lon_v,1);
xi_v=size(lon_v,2);
DEP_rho=nc{'h'}(:);
ANG_rho=nc{'angle'}(:);
close(nc);
DEP_u=(DEP_rho(:, 1:xi_u)+DEP_rho(:, 2:xi_rho))/2;
DEP_v=(DEP_rho(1:eta_v, :)+DEP_rho(2:eta_rho, :))/2;
ANG_u=(ANG_rho(:, 1:xi_u)+ANG_rho(:, 2:xi_rho))/2;
ANG_v=(ANG_rho(1:eta_v, :)+ANG_rho(2:eta_rho, :))/2;

nb=0;

% the rho points.

lon_rho_west = lon_rho(1:eta_rho,1);
lat_rho_west = lat_rho(1:eta_rho,1);
dep_rho_west = DEP_rho(1:eta_rho,1);
ang_rho_west = ANG_rho(1:eta_rho,1);
clip = find(mask_rho(1:eta_rho,1)==0);
lon_rho_west(clip) = [];
lat_rho_west(clip) = [];
dep_rho_west(clip) = [];
ang_rho_west(clip) = [];
nb_rho_west=size(lon_rho_west(:), 1);
STA_rho_west=zeros(nb_rho_west, 1);
ORI_rho_west=ones(nb_rho_west, 1);
rho_west_info=ones(eta_rho,1);
rho_west_info(clip)=0;
disp(['We have ' num2str(nb_rho_west) ' west rho points'])
rho_west_begin=nb+1;
rho_west_end=nb+nb_rho_west;
nb=nb+nb_rho_west;


lon_rho_south = lon_rho(1,1:xi_rho);
lat_rho_south = lat_rho(1,1:xi_rho);
dep_rho_south = DEP_rho(1,1:xi_rho);
ang_rho_south = ANG_rho(1,1:xi_rho);
clip = find(mask_rho(1,1:xi_rho)==0);
lon_rho_south(clip) = [];
lat_rho_south(clip) = [];
dep_rho_south(clip) = [];
ang_rho_south(clip) = [];
nb_rho_south=size(lon_rho_south(:), 1);
STA_rho_south=zeros(nb_rho_south, 1);
ORI_rho_south=2*ones(nb_rho_south, 1);
rho_south_info=ones(xi_rho,1);
rho_south_info(clip)=0;
disp(['We have ' num2str(nb_rho_south) ' south rho points'])
rho_south_begin=nb+1;
rho_south_end=nb+nb_rho_south;
nb=nb+nb_rho_south;

lon_rho_east = lon_rho(1:eta_rho,xi_rho);
lat_rho_east = lat_rho(1:eta_rho,xi_rho);
dep_rho_east = DEP_rho(1:eta_rho,xi_rho);
ang_rho_east = ANG_rho(1:eta_rho,xi_rho);
clip = find(mask_rho(1:eta_rho,xi_rho)==0);
lon_rho_east(clip) = [];
lat_rho_east(clip) = [];
dep_rho_east(clip) = [];
ang_rho_east(clip) = [];
nb_rho_east=size(lon_rho_east(:), 1);
STA_rho_east=zeros(nb_rho_east, 1);
ORI_rho_east=3*ones(nb_rho_east, 1);
rho_east_info=ones(eta_rho,1);
rho_east_info(clip)=0;
disp(['We have ' num2str(nb_rho_east) ' east rho points'])
rho_east_begin=nb+1;
rho_east_end=nb+nb_rho_east;
nb=nb+nb_rho_east;

lon_rho_north = lon_rho(eta_rho,1:xi_rho);
lat_rho_north = lat_rho(eta_rho,1:xi_rho);
dep_rho_north = DEP_rho(eta_rho,1:xi_rho);
ang_rho_north = ANG_rho(eta_rho,1:xi_rho);
clip = find(mask_rho(eta_rho,1:xi_rho)==0);
lon_rho_north(clip) = [];
lat_rho_north(clip) = [];
dep_rho_north(clip) = [];
ang_rho_north(clip) = [];
nb_rho_north=size(lon_rho_north(:), 1);
STA_rho_north=zeros(nb_rho_north, 1);
ORI_rho_north=4*ones(nb_rho_north, 1);
rho_north_info=ones(xi_rho,1);
rho_north_info(clip)=0;
disp(['We have ' num2str(nb_rho_north) ' north rho points'])
rho_north_begin=nb+1;
rho_north_end=nb+nb_rho_north;
nb=nb+nb_rho_north;

% the u-points.

lon_u_west = lon_u(1:eta_u,1);
lat_u_west = lat_u(1:eta_u,1);
dep_u_west = DEP_u(1:eta_u,1);
ang_u_west = ANG_u(1:eta_u,1);
clip = find(mask_u(1:eta_u,1)==0);
lon_u_west(clip) = [];
lat_u_west(clip) = [];
dep_u_west(clip) = [];
ang_u_west(clip) = [];
nb_u_west=size(lon_u_west(:), 1);
STA_u_west=1*ones(nb_u_west, 1);
ORI_u_west=1*ones(nb_u_west, 1);
u_west_info=ones(eta_u,1);
u_west_info(clip)=0;
disp(['We have ' num2str(nb_u_west) ' west u-points'])
u_west_begin=nb+1;
u_west_end=nb+nb_u_west;
nb=nb+nb_u_west;

lon_u_south = lon_u(1,1:xi_u);
lat_u_south = lat_u(1,1:xi_u);
dep_u_south = DEP_u(1,1:xi_u);
ang_u_south = ANG_u(1,1:xi_u);
clip = find(mask_u(1,1:xi_u)==0);
lon_u_south(clip) = [];
lat_u_south(clip) = [];
dep_u_south(clip) = [];
ang_u_south(clip) = [];
nb_u_south=size(lon_u_south(:), 1);
STA_u_south=1*ones(nb_u_south, 1).*1;
ORI_u_south=2*ones(nb_u_south, 1);
u_south_info=ones(xi_u,1);
u_south_info(clip)=0;
disp(['We have ' num2str(nb_u_south) ' south u-points'])
u_south_begin=nb+1;
u_south_end=nb+nb_u_south;
nb=nb+nb_u_south;

lon_u_east = lon_u(1:eta_u,xi_u);
lat_u_east = lat_u(1:eta_u,xi_u);
dep_u_east = DEP_u(1:eta_u,xi_u);
ang_u_east = ANG_u(1:eta_u,xi_u);
clip = find(mask_u(1:eta_u,xi_u)==0);
lon_u_east(clip) = [];
lat_u_east(clip) = [];
dep_u_east(clip) = [];
ang_u_east(clip) = [];
nb_u_east=size(lon_u_east(:), 1);
STA_u_east=1*ones(nb_u_east, 1);
ORI_u_east=3*ones(nb_u_east, 1);
u_east_info=ones(eta_u,1);
u_east_info(clip)=0;
disp(['We have ' num2str(nb_u_east) ' east u-points'])
u_east_begin=nb+1;
u_east_end=nb+nb_u_east;
nb=nb+nb_u_east;

lon_u_north = lon_u(eta_u,1:xi_u);
lat_u_north = lat_u(eta_u,1:xi_u);
dep_u_north = DEP_u(eta_u,1:xi_u);
ang_u_north = ANG_u(eta_u,1:xi_u);
clip = find(mask_u(eta_u,1:xi_u)==0);
lon_u_north(clip) = [];
lat_u_north(clip) = [];
dep_u_north(clip) = [];
ang_u_north(clip) = [];
nb_u_north=size(lon_u_north(:), 1);
STA_u_north=1*ones(nb_u_north, 1);
ORI_u_north=4*ones(nb_u_north, 1);
u_north_info=ones(eta_u,1);
u_north_info(clip)=0;
disp(['We have ' num2str(nb_u_north) ' north u-points'])
u_north_begin=nb+1;
u_north_end=nb+nb_u_north;
nb=nb+nb_u_north;

% the v-points

lon_v_west = lon_v(1:eta_v,1);
lat_v_west = lat_v(1:eta_v,1);
dep_v_west = DEP_v(1:eta_v,1);
ang_v_west = ANG_v(1:eta_v,1);
clip = find(mask_v(1:eta_v,1)==0);
lon_v_west(clip) = [];
lat_v_west(clip) = [];
dep_v_west(clip) = [];
ang_v_west(clip) = [];
nb_v_west=size(lon_v_west(:), 1);
STA_v_west=2*ones(nb_v_west, 1);
ORI_v_west=1*ones(nb_v_west, 1);
v_west_info=ones(eta_v,1);
v_west_info(clip)=0;
disp(['We have ' num2str(nb_v_west) ' west v-points'])
v_west_begin=nb+1;
v_west_end=nb+nb_v_west;
nb=nb+nb_v_west;

lon_v_south = lon_v(1,1:xi_v);
lat_v_south = lat_v(1,1:xi_v);
dep_v_south = DEP_v(1,1:xi_v);
ang_v_south = ANG_v(1,1:xi_v);
clip = find(mask_v(1,1:xi_v)==0);
lon_v_south(clip) = [];
lat_v_south(clip) = [];
dep_v_south(clip) = [];
ang_v_south(clip) = [];
nb_v_south=size(lon_v_south(:), 1);
STA_v_south=2*ones(nb_v_south, 1);
ORI_v_south=2*ones(nb_v_south, 1);
v_south_info=ones(xi_v,1);
v_south_info(clip)=0;
disp(['We have ' num2str(nb_v_south) ' south v-points'])
v_south_begin=nb+1;
v_south_end=nb+nb_v_south;
nb=nb+nb_v_south;

lon_v_east = lon_v(1:eta_v,xi_v);
lat_v_east = lat_v(1:eta_v,xi_v);
dep_v_east = DEP_v(1:eta_v,xi_v);
ang_v_east = ANG_v(1:eta_v,xi_v);
clip = find(mask_v(1:eta_v,xi_v)==0);
lon_v_east(clip) = [];
lat_v_east(clip) = [];
dep_v_east(clip) = [];
ang_v_east(clip) = [];
nb_v_east=size(lon_v_east(:), 1);
STA_v_east=2*ones(nb_v_east, 1);
ORI_v_east=3*ones(nb_v_east, 1);
v_east_info=ones(eta_v,1);
v_east_info(clip)=0;
disp(['We have ' num2str(nb_v_east) ' east v-points'])
v_east_begin=nb+1;
v_east_end=nb+nb_v_east;
nb=nb+nb_v_east;

lon_v_north = lon_v(eta_v,1:xi_v);
lat_v_north = lat_v(eta_v,1:xi_v);
dep_v_north = DEP_v(eta_v,1:xi_v);
ang_v_north = ANG_v(eta_v,1:xi_v);
clip = find(mask_v(eta_v,1:xi_v)==0);
lon_v_north(clip) = [];
lat_v_north(clip) = [];
dep_v_north(clip) = [];
ang_v_north(clip) = [];
nb_v_north=size(lon_v_north(:), 1);
STA_v_north=2*ones(nb_v_north, 1);
ORI_v_north=4*ones(nb_v_north, 1);
v_north_info=ones(xi_v,1);
v_north_info(clip)=0;
disp(['We have ' num2str(nb_v_north) ' north v-points'])
v_north_begin=nb+1;
v_north_end=nb+nb_v_north;
nb=nb+nb_v_north;

lon_nest_rho = [lon_rho_west(:); lon_rho_south(:); ...
		lon_rho_east(:); lon_rho_north(:);];
lon_nest_u = [lon_u_west(:); lon_u_south(:); ...
	      lon_u_east(:); lon_u_north(:);];
lon_nest_v = [lon_v_west(:); lon_v_south(:); ...
	      lon_v_east(:); lon_v_north(:);];
lon_nest=[lon_nest_rho(:); lon_nest_u(:); lon_nest_v(:);];
%
lat_nest_rho = [lat_rho_west(:); lat_rho_south(:); ...
		lat_rho_east(:); lat_rho_north(:);];
lat_nest_u = [lat_u_west(:); lat_u_south(:); ...
	      lat_u_east(:); lat_u_north(:);];
lat_nest_v = [lat_v_west(:); lat_v_south(:); ...
	      lat_v_east(:); lat_v_north(:);];
lat_nest=[lat_nest_rho(:); lat_nest_u(:); lat_nest_v(:);];
%
dep_nest_rho = [dep_rho_west(:); dep_rho_south(:); ...
		dep_rho_east(:); dep_rho_north(:);];
dep_nest_u = [dep_u_west(:); dep_u_south(:); ...
	      dep_u_east(:); dep_u_north(:);];
dep_nest_v = [dep_v_west(:); dep_v_south(:); ...
	      dep_v_east(:); dep_v_north(:);];
dep_nest=[dep_nest_rho(:); dep_nest_u(:); dep_nest_v(:);];
%
ang_nest_rho = [ang_rho_west(:); ang_rho_south(:); ...
		ang_rho_east(:); ang_rho_north(:);];
ang_nest_u = [ang_u_west(:); ang_u_south(:); ...
	      ang_u_east(:); ang_u_north(:);];
ang_nest_v = [ang_v_west(:); ang_v_south(:); ...
	      ang_v_east(:); ang_v_north(:);];
ang_nest=[ang_nest_rho(:); ang_nest_u(:); ang_nest_v(:);];
%
STA_nest_rho = [STA_rho_west(:); STA_rho_south(:); ...
		STA_rho_east(:); STA_rho_north(:);];
STA_nest_u = [STA_u_west(:); STA_u_south(:); ...
	      STA_u_east(:); STA_u_north(:);];
STA_nest_v = [STA_v_west(:); STA_v_south(:); ...
	      STA_v_east(:); STA_v_north(:);];
STA_nest=[STA_nest_rho(:); STA_nest_u(:); STA_nest_v(:);];
%
ORI_nest_rho = [ORI_rho_west(:); ORI_rho_south(:); ...
		ORI_rho_east(:); ORI_rho_north(:);];
ORI_nest_u = [ORI_u_west(:); ORI_u_south(:); ...
	      ORI_u_east(:); ORI_u_north(:);];
ORI_nest_v = [ORI_v_west(:); ORI_v_south(:); ...
	      ORI_v_east(:); ORI_v_north(:);];
ORI_nest=[ORI_nest_rho(:); ORI_nest_u(:); ORI_nest_v(:);];
%
clear('ListType');
clear('ListComment');
nbPoint=size(STA_nest,1);
for iPoint=1:nbPoint
  clear('TheRec');
  TheRec.ang=ang_nest(iPoint,1);
  if (STA_nest(iPoint,1) == 0)
    strType='rho';
  elseif (STA_nest(iPoint,1) == 1)
    strType='u';
  elseif (STA_nest(iPoint,1) == 2)
    strType='v';
  else
    disp('We have a problem');
    keyboard;
  end;
  TheRec.Status=strType;
  if (ORI_nest(iPoint,1) == 1)
    strOri='west';
  elseif (ORI_nest(iPoint,1) == 2)
    strOri='south';
  elseif (ORI_nest(iPoint,1) == 3)
    strOri='east';
  elseif (ORI_nest(iPoint,1) == 4)
    strOri='north';
  else
    disp('We have a problem 2');
    keyboard;
  end;
  TheRec.orientation=strOri;
  TheRec.dep=dep_nest(iPoint,1);
  TheRec.lon=lon_nest(iPoint,1);
  TheRec.lat=lat_nest(iPoint,1);
  ListType{iPoint}=TheRec;
  ListComment{iPoint}=[strType '-point oriented ' strOri];
end;

%
% TheRecord of full informations
%
TheRecord.rho_west_info=rho_west_info;
TheRecord.rho_east_info=rho_east_info;
TheRecord.rho_north_info=rho_north_info;
TheRecord.rho_south_info=rho_south_info;
TheRecord.u_west_info=u_west_info;
TheRecord.u_east_info=u_east_info;
TheRecord.u_north_info=u_north_info;
TheRecord.u_south_info=u_south_info;
TheRecord.v_west_info=v_west_info;
TheRecord.v_east_info=v_east_info;
TheRecord.v_north_info=v_north_info;
TheRecord.v_south_info=v_south_info;
%
TheRecord.rho_west_begin=rho_west_begin;
TheRecord.rho_west_end=rho_west_end;
TheRecord.rho_east_begin=rho_east_begin;
TheRecord.rho_east_end=rho_east_end;
TheRecord.rho_north_begin=rho_north_begin;
TheRecord.rho_north_end=rho_north_end;
TheRecord.rho_south_begin=rho_south_begin;
TheRecord.rho_south_end=rho_south_end;
%
TheRecord.u_west_begin=u_west_begin;
TheRecord.u_west_end=u_west_end;
TheRecord.u_east_begin=u_east_begin;
TheRecord.u_east_end=u_east_end;
TheRecord.u_north_begin=u_north_begin;
TheRecord.u_north_end=u_north_end;
TheRecord.u_south_begin=u_south_begin;
TheRecord.u_south_end=u_south_end;
%
TheRecord.v_west_begin=v_west_begin;
TheRecord.v_west_end=v_west_end;
TheRecord.v_east_begin=v_east_begin;
TheRecord.v_east_end=v_east_end;
TheRecord.v_north_begin=v_north_begin;
TheRecord.v_north_end=v_north_end;
TheRecord.v_south_begin=v_south_begin;
TheRecord.v_south_end=v_south_end;
%
TheRecord.nb_rho_west=nb_rho_west;
TheRecord.nb_rho_east=nb_rho_east;
TheRecord.nb_rho_north=nb_rho_north;
TheRecord.nb_rho_south=nb_rho_south;
TheRecord.nb_u_west=nb_u_west;
TheRecord.nb_u_east=nb_u_east;
TheRecord.nb_u_north=nb_u_north;
TheRecord.nb_u_south=nb_u_south;
TheRecord.nb_v_west=nb_v_west;
TheRecord.nb_v_east=nb_v_east;
TheRecord.nb_v_north=nb_v_north;
TheRecord.nb_v_south=nb_v_south;
%
TheRecord.lon_nest=lon_nest;
TheRecord.lat_nest=lat_nest;
TheRecord.dep_nest=dep_nest;
TheRecord.ang_nest=ang_nest;
TheRecord.STA_nest=STA_nest;
TheRecord.ORI_nest=ORI_nest;
TheRecord.ListType=ListType;
TheRecord.ListComment=ListComment;
%
TheRecord.lon_nest_rho=lon_nest_rho;
TheRecord.lat_nest_rho=lat_nest_rho;
TheRecord.lon_nest_u=lon_nest_u;
TheRecord.lat_nest_u=lat_nest_u;
TheRecord.lon_nest_v=lon_nest_v;
TheRecord.lat_nest_v=lat_nest_v;
