%
% the relevant grids.
%
BigGridFile='ncom2_376x136_TOP7samp_simple_r0_20.nc';
SmallGridFile='ROVINJ_700m_grd_V4.nc';

%
% their vertical parametrization
%

AddiRecordInfo.BigThetaS=3;
AddiRecordInfo.BigThetaB=0.35;
AddiRecordInfo.Nbig=20;
AddiRecordInfo.Bighc=3;

AddiRecordInfo.SmaThetaS=3;
AddiRecordInfo.SmaThetaB=0.35;
AddiRecordInfo.Nsma=20;
AddiRecordInfo.Smahc=3;

%
% the method used for interpolation:
% --UseSpMat=1 for using sparse matrix.
% --UseMemEff=1 for using memory efficient method
% If both selected the results are compared.
%

UseSpMat=1;
UseMemEff=1;

%
% the initial file we want to interpolate
%

InitialFileBig='init_2007_12_09_00_00_00_ncom2.nc';

%
% the place where we store matrix informations.
%

PrefixInterpData='ncomTOP7sampR0_20_RovinjV4';

%
% loading the information from the big grid.
%

nc=netcdf(InitialFileBig, 'nowrite');
ZETAbig=nc{'zeta'}(1, :, :);
UBARbig=nc{'ubar'}(1, :, :);
VBARbig=nc{'vbar'}(1, :, :);
Ubig=nc{'u'}(1, :, :, :);
Vbig=nc{'v'}(1, :, :, :);
TEMPbig=nc{'temp'}(1, :, :, :);
SALTbig=nc{'salt'}(1, :, :, :);
close(nc);

%
% Computing the initial informations for interpolating
% between the two grids.
%

FileSave=[PrefixInterpData '_basicArr.mat'];
if (IsExistingFile(FileSave) == 0)
  tic
  TotalArray=InterpolGetTotalArray_R2R_fields(...
      BigGridFile, SmallGridFile);
  t=toc;
  disp(['TotalArray computed, running time ' num2str(t)]);
  save(FileSave, 'TotalArray');
else
  load(FileSave, 'TotalArray');
end;
disp('We have TotalArray');

%
% Computing the sparse matrices for sparse matrix computation
% (Note that in practice you might wish to compute only a subset
%  of them, for examepl ArrayBigSma2D and ArrayBigSma2Duv)
%

if (UseSpMat == 1)
  FileSave=[PrefixInterpData '_addiArr.mat'];
  if (IsExistingFile(FileSave) == 0)
    tic
    ArrayBigSma2D=InterpolGetSpMat_R2R_2Dfield(TotalArray);
    ArrayBigSma2Duv=InterpolGetSpMat_R2R_2Duvfield(TotalArray);
    ArrayBigSma3D=InterpolGetSpMat_R2R_3Dfield(TotalArray, AddiRecordInfo);
    ArrayBigSma3Duv=InterpolGetSpMat_R2R_3Duvfield(...
	TotalArray, AddiRecordInfo);
    t=toc;
    disp(['Sparse matrices computed, running time ' num2str(t)]);
    save(FileSave, 'ArrayBigSma2D', 'ArrayBigSma2Duv', 'ArrayBigSma3D', 'ArrayBigSma3Duv');
  else
    load(FileSave, 'ArrayBigSma2D', 'ArrayBigSma2Duv', 'ArrayBigSma3D', 'ArrayBigSma3Duv');
  end;
end;
disp('We have ArrayBigSma*');

%
% Computation of the interpolation using sparse matrices.
%

if (UseSpMat == 1)
  tic
  ZETAsma_sm=InterpolSpMat_R2R_2Dfield(ArrayBigSma2D, ZETAbig);
  disp('zeta done');
  ZETAsma=ZETAsma_sm;
  %
  TEMPsma_sm=InterpolSpMat_R2R_3Dfield(ArrayBigSma3D, TEMPbig);
  disp('temp done');
  TEMPsma=TEMPsma_sm;
  %
  SALTsma_sm=InterpolSpMat_R2R_3Dfield(ArrayBigSma3D, SALTbig);
  disp('salt done');
  SALTsma=SALTsma_sm;
  %
  [UBARsma_sm, VBARsma_sm]=InterpolSpMat_R2R_2Duvfield(...
      ArrayBigSma2Duv, UBARbig, VBARbig);
  disp('ubar/vbar done');
  UBARsma=UBARsma_sm;
  VBARsma=VBARsma_sm;
  %
  [Usma_sm, Vsma_sm]=InterpolSpMat_R2R_3Duvfield(...
      ArrayBigSma3Duv, Ubig, Vbig);
  disp('u/v done');
  Usma=Usma_sm;
  Vsma=Vsma_sm;
  t=toc;
  disp(['Sparse matrix interpolations done. runtime = ' num2str(t)]);
end;

%
% Computation of the interpolation using memory efficient 
% procedures. (Normally slower, but run several times to really
% get an idea)
%

if (UseMemEff == 1)
  tic
  ZETAsma_me=InterpolMemEff_R2R_2Dfield(TotalArray, ZETAbig);
  disp('zeta done');
  ZETAsma=ZETAsma_me;
  %
  TEMPsma_me=InterpolMemEff_R2R_3Dfield(...
      TotalArray, AddiRecordInfo, TEMPbig);
  disp('temp done');
  TEMPsma=TEMPsma_me;
  %
  SALTsma_me=InterpolMemEff_R2R_3Dfield(...
      TotalArray, AddiRecordInfo, SALTbig);
  SALTsma=SALTsma_me;
  %
  [UBARsma_me, VBARsma_me]=InterpolMemEff_R2R_2Duvfield(...
      TotalArray, UBARbig, VBARbig);
  disp('ubar/vbar done');
  UBARsma=UBARsma_me;
  VBARsma=VBARsma_me;
  %
  [Usma_me, Vsma_me]=InterpolMemEff_R2R_3Duvfield(...
      TotalArray, AddiRecordInfo, Ubig, Vbig);
  disp('u/v done');
  Usma=Usma_me;
  Vsma=Vsma_me;
  t=toc;
  disp(['Memory Efficient interpolations done. runtime = ' num2str(t)]);
end;

%
% Comparing the results (Equal if no bug but this does not 
% guarantee correctness in dynamical sense, which is another story)
%

if (UseMemEff == 1 && UseSpMat == 1)
  deltaZETA=sum(abs(ZETAsma_me(:)-ZETAsma_sm(:)));
  disp(['deltaZETA=' num2str(deltaZETA)]);
  %
  deltaTEMP=sum(abs(TEMPsma_me(:)-TEMPsma_sm(:)));
  disp(['deltaTEMP=' num2str(deltaTEMP)]);
  %
  deltaSALT=sum(abs(SALTsma_me(:)-SALTsma_sm(:)));
  disp(['deltaSALT=' num2str(deltaSALT)]);
  %
  deltaUVBAR=sum(abs(UBARsma_me(:)-UBARsma_sm(:)))+...
      sum(abs(VBARsma_me(:)-VBARsma_sm(:)));
  disp(['deltaUVBAR=' num2str(deltaUVBAR)]);
  %
  deltaUV=sum(abs(Usma_me(:)-Usma_sm(:)))+...
      sum(abs(Vsma_me(:)-Vsma_sm(:)));
  disp(['deltaUV=' num2str(deltaUV)]);
  %
  disp('Comparison finished');
end;
