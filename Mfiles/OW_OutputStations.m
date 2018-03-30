function OW_OutputStations(StationFile, ListLon, ListLat, ListComment)


nbStation=size(ListLon, 1);
disp(['writing ' num2str(nbStation) ' points in ' StationFile]);
RemoveFileIfExist(StationFile); % actually not needed, 'w' erase it.
fid = fopen(StationFile,'w');
OW_HeaderStationFile(fid);
str=['    NSTATION ==  ' int2str(nbStation)];
fprintf(fid,'%s\n',str);
str='POS =  GRID  FLAG      X-POS       Y-POS     COMMENT';
fprintf(fid,'%s\n',str);
for iSta=1:nbStation
  str=['         1    1      ' ...
       num2str(ListLon(iSta, 1)) 'd0  ' ...
       num2str(ListLat(iSta,1)) 'd0 ! ' ...
       num2str(iSta) ' : ' ListComment{iSta}];
  fprintf(fid,'%s\n',str);
end
fclose(fid);
