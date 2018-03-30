GridFile='ncom2_376x136_TOP7samp_simple_r0_20.nc';
InitialFile='init_2007_12_09_00_00_00_ncom2.nc';

AddiRecordVertDisc.ThetaS=3;
AddiRecordVertDisc.ThetaB=0.35;
AddiRecordVertDisc.N=20;
AddiRecordVertDisc.hc=3;

%
% the longitude (totally arbitrary...)
% 
ListLonADCP=[14+04/60; 15+10/60];
ListLatADCP=[43+55/60; 42+58/60];

%
% The depth, 4 levels for the first
%            3 levels for the second
ADCPverticalLevels=[ -5 -10 -12 -20;
		    -10 -15 -20 NaN];

%
% creating the total array used later on.
%
TotalArray=ADCP_R2A_GetTotalArray(...
    GridFile, AddiRecordVertDisc, ...
    ListLonADCP, ListLatADCP, ADCPverticalLevels);

%
% Reading data from the initial file
%
nc=netcdf(InitialFile, 'nowrite');
ZETAroms=nc{'zeta'}(1, :, :);
UBARroms=nc{'ubar'}(1, :, :);
VBARroms=nc{'vbar'}(1, :, :);
Uroms=nc{'u'}(1, :, :, :);
Vroms=nc{'v'}(1, :, :, :);
TEMProms=nc{'temp'}(1, :, :, :);
SALTroms=nc{'salt'}(1, :, :, :);
close(nc);

%
% the interpolation itself.
% 

ZETAadcp=ADCP_R2A_Interpolation2Dfield(TotalArray, ZETAroms);
[UVbar_east, UVbar_north]=ADCP_R2A_Interpolation2Duvfield(...
    TotalArray, UBARroms, VBARroms);
TEMPadcp=ADCP_R2A_Interpolation3Dfield(TotalArray, TEMProms);
SALTadcp=ADCP_R2A_Interpolation3Dfield(TotalArray, SALTroms);
[UV_east, UV_north]=ADCP_R2A_Interpolation3Duvfield(...
    TotalArray, Uroms, Vroms);

