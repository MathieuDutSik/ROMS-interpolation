function TheStr=StringNumber(nb, nbDigit)


if (nb > 10^(nbDigit))
  disp('the number is greater than the allowed size');
  disp('please correct');
  keyboard;
end;

idx=1;
while(1)
  if (nb < 10^(idx))
    TheStr='';
    for U=1:nbDigit-idx
      TheStr=[TheStr '0'];
    end
    TheStr=[TheStr num2str(nb)];
    break;
  end;
  idx=idx+1;
end;


