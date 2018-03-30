function CorrespMat=GetCorrespMatrix(nbdim, ListDim)

if (nbdim == 2)
  dim1=ListDim(1,1);
  dim2=ListDim(2,1);
  ListCoord1=zeros(dim1, dim2);
  ListCoord2=zeros(dim1, dim2);
  for iC1=1:dim1
    for iC2=1:dim2
      ListCoord1(iC1, iC2)=iC1;
      ListCoord2(iC1, iC2)=iC2;
    end;
  end;
  CorrespMat=zeros(dim1, dim2);
  ListCoord1exp=ListCoord1(:);
  ListCoord2exp=ListCoord2(:);
  TheSiz=size(ListCoord2exp,1);
  for i=1:TheSiz
    iC1=ListCoord1exp(i,1);
    iC2=ListCoord2exp(i,1);
    CorrespMat(iC1, iC2)=i;
  end;
elseif (nbdim == 3)
  dim1=ListDim(1,1);
  dim2=ListDim(2,1);
  dim3=ListDim(3,1);
  ListCoord1=zeros(dim1, dim2, dim3);
  ListCoord2=zeros(dim1, dim2, dim3);
  ListCoord3=zeros(dim1, dim2, dim3);
  for iC1=1:dim1
    for iC2=1:dim2
      for iC3=1:dim3
	ListCoord1(iC1, iC2, iC3)=iC1;
	ListCoord2(iC1, iC2, iC3)=iC2;
	ListCoord3(iC1, iC2, iC3)=iC3;
      end;
    end;
  end;
  CorrespMat=zeros(dim1, dim2, dim3);
  ListCoord1exp=ListCoord1(:);
  ListCoord2exp=ListCoord2(:);
  ListCoord3exp=ListCoord3(:);
  TheSiz=size(ListCoord3exp,1);
  for i=1:TheSiz
    iC1=ListCoord1exp(i,1);
    iC2=ListCoord2exp(i,1);
    iC3=ListCoord3exp(i,1);
    CorrespMat(iC1, iC2, iC3)=i;
  end;
else
  disp('We need to program more');
  keyboard;
end;