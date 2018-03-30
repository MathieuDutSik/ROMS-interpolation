function dist=GeodesicDistance(Pair1, Pair2)
lon1=pi*Pair1(1)/180;
lat1=pi*Pair1(2)/180;
x1=cos(lon1)*cos(lat1);
y1=sin(lon1)*cos(lat1);
z1=sin(lat1);

lon2=pi*Pair2(1)/180;
lat2=pi*Pair2(2)/180;
x2=cos(lon2)*cos(lat2);
y2=sin(lon2)*cos(lat2);
z2=sin(lat2);
scalprod=x1*x2+y1*y2+z1*z2;
dist=acos(scalprod);
%dist=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
