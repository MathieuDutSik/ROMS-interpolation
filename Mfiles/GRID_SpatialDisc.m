function result=GRID_SpatialDisc(MSK_s, LON_s, LAT_s)

[ETA_s, XI_s]=size(MSK_s);
EarthRadius=6367;

nbDist=0;
ListDist=zeros(0,1);
for i=1:ETA_s-1
  for j=1:XI_s
    i1=i;
    j1=j;
    i2=i+1;
    j2=j;
    if (MSK_s(i1, j1) == 1 && MSK_s(i2, j2) == 1)
      nbDist=nbDist+1;
      ePair1=[LON_s(i1, j1) LAT_s(i1,j1)];
      ePair2=[LON_s(i2, j2) LAT_s(i2,j2)];
      ListDist(nbDist,1)=GeodesicDistance(ePair1, ePair2)*EarthRadius;
    end;
  end;
end;
result.MaxDist_ETA=max(ListDist);
result.MinDist_ETA=min(ListDist);
result.AveDist_ETA=sum(ListDist(:))/nbDist;

nbDistTotal=nbDist;
sumDistTotal=sum(ListDist(:));

nbDist=0;
ListDist=zeros(0,1);
for i=1:ETA_s
  for j=1:XI_s-1
    i1=i;
    j1=j;
    i2=i;
    j2=j+1;
    if (MSK_s(i1, j1) == 1 && MSK_s(i2, j2) == 1)
      nbDist=nbDist+1;
      ePair1=[LON_s(i1, j1) LAT_s(i1,j1)];
      ePair2=[LON_s(i2, j2) LAT_s(i2,j2)];
      ListDist(nbDist,1)=GeodesicDistance(ePair1, ePair2)*EarthRadius;
    end;
  end;
end;
result.MaxDist_XI=max(ListDist);
result.MinDist_XI=min(ListDist);
result.AveDist_XI=sum(ListDist(:))/nbDist;

nbDistTotal=nbDistTotal+nbDist;
sumDistTotal=sumDistTotal+sum(ListDist(:));

result.AveDist=sumDistTotal/nbDistTotal;
