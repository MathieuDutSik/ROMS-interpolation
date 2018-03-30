function qlev=interp_vertical(zref,qref,zlev)

% zref, qref : input
%       zlev : output


% zadaje se kao [qlev]=intepr_vertical(zref,qref,zlev);
% koristi se basis1d1; ako je zlev ispod zadnjeg nivoa zref tada se
% jednostavno kopira vrijednost sa zadnjeg nivoa qref(nnv), isto tako ako je zlev
% iznad prvog nivoa (cisto sumljam) tada se kopira qref(1) sa prvog nivoa
nnv=length(zlev);
H=find(isnan(qref) == 0);
if (sum(H) == 1)
  idx=H(1,1);
  qlev=qref(idx)*ones(nnv,1);
else
  qlev=NaN*zeros(nnv,1);
  %moram imati pozitivne vrijednosti za dubine tj sortirane 
  for i=1:nnv
    [b2,n2,b1,n1]=basis1d1(zref,zlev(i));
    qlev(i)=qref(n2)*b2+qref(n1)*b1;
  end; 
end;

