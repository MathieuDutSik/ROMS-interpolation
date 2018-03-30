function TEMPusual=Field3D_wp2usual(MSK_rho, N, TEMPwp)

[eta_rho, xi_rho]=size(MSK_rho);
TEMPusual=zeros(N, eta_rho, xi_rho);
idx=0;
for iN=1:N
  for iEta=1:eta_rho
    for iXi=1:xi_rho
      if (MSK_rho(iEta, iXi) == 1)
	idx=idx+1;
	TEMPusual(iN, iEta,iXi)=TEMPwp(1,idx);
      end;
    end;
  end;
end;
