function TotalArray=InterpolGetTotalArray_R2R_fields(...
    BigGridFile, SmallGridFile)
%
GrdArrBig=GRID_GetArray(BigGridFile);
GrdArrSma=GRID_GetArray(SmallGridFile);
%
TotalArray=InterpolGetTotalArray_R2R_fields_kernel(...
    GrdArrBig, GrdArrSma);
