function ZETAsma=InterpolSpMat_R2R_2Dfield(ZETAbig, TotalArray)



ZETAbigExp=ZETAbig(:);
ZETAsmaExp=TotalArray.ArrayBigSma2D.SPmat*ZETAbigExp;
ZETAsma=reshape(...
    ZETAsmaExp, ...
    TotalArray.ArrayBigSma2D.eta_rho_sma, ...
    TotalArray.ArrayBigSma2D.xi_rho_sma);
