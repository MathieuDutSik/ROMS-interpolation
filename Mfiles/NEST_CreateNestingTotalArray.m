function NestingTotalArray=NEST_CreateNestingTotalArray(...
    BigGridFile, SmallGridFile, AddiRecordInfo, UseSpMat)

GrdArrBig=GRID_GetArray(BigGridFile);
GrdArrSma=GRID_GetArray(SmallGridFile);

eta_rho_sma=GrdArrSma.eta_rho;
xi_rho_sma=GrdArrSma.xi_rho;
eta_u_sma=GrdArrSma.eta_u;
xi_u_sma=GrdArrSma.xi_u;
eta_v_sma=GrdArrSma.eta_v;
xi_v_sma=GrdArrSma.xi_v;

MSKsma_rho=GrdArrSma.MSK_rho;
MSKsma_rho(2:eta_rho_sma-1,2:xi_rho_sma-1)=0;
GrdArrSma.MSK_rho=MSKsma_rho;
%
MSKsma_u=GrdArrSma.MSK_u;
MSKsma_u(2:eta_u_sma-1,2:xi_u_sma-1)=0;
GrdArrSma.MSK_u=MSKsma_u;
%
MSKsma_v=GrdArrSma.MSK_v;
MSKsma_v(2:eta_v_sma-1,2:xi_v_sma-1)=0;
GrdArrSma.MSK_v=MSKsma_v;
%
NestingTotalArray=InterpolGetTotalArray_R2R_fields_kernel(...
    GrdArrBig, GrdArrSma);
%
NestingTotalArray.GrdArrBig=GrdArrBig;
%
if (UseSpMat == 1)
  ArrayBigSma2D=InterpolGetSpMat_R2R_2Dfield(NestingTotalArray);
  ArrayBigSma2Duv=InterpolGetSpMat_R2R_2Duvfield(NestingTotalArray);
  ArrayBigSma3D=InterpolGetSpMat_R2R_3Dfield(...
      NestingTotalArray, AddiRecordInfo);
  ArrayBigSma3Duv=InterpolGetSpMat_R2R_3Duvfield(...
      NestingTotalArray, AddiRecordInfo);
  NestingTotalArray.ArrayBigSma2D=ArrayBigSma2D;
  NestingTotalArray.ArrayBigSma2Duv=ArrayBigSma2Duv;
  NestingTotalArray.ArrayBigSma3D=ArrayBigSma3D;
  NestingTotalArray.ArrayBigSma3Duv=ArrayBigSma3Duv;
end;
%
NestingTotalArray.AddiRecordInfo=AddiRecordInfo;
NestingTotalArray.UseSpMat=UseSpMat;
