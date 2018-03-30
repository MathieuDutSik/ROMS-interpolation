function NewField=Field3D_PutToFillVal(InputMSK, TEMP_rho, FillVal)

Kland=find(InputMSK == 0);
[N, eta_rho, xi_rho]=size(TEMP_rho);
NewField=zeros(N, eta_rho, xi_rho);
for i=1:N
  Ksel=squeeze(TEMP_rho(i, :, :));
  Ksel(Kland)=FillVal;
  NewField(i, :, :)=Ksel;
end;
