function [Sc_w, Cs_w, Sc_r, Cs_r, TheMult]=...
    GRID_GetSc_Cs(N, theta_s, theta_b)

cff=1/N;

Sc_w=zeros(1,N);
Cs_w=zeros(1,N);
Sc_r=zeros(1,N);
Cs_r=zeros(1,N);

for k=1:N
  Sc_w(1,k)=(k-N)/N;
  Sc_r(1,k)=(k-N-0.5)/N;
  cff1=1/sinh(theta_s);
  cff2=1/(2*tanh(theta_s/2));
  if (theta_s>0)
    Cs_w(1,k)=(1-theta_b)*cff1*sinh(theta_s*Sc_w(1,k))+...
	    theta_b*(cff2*tanh(theta_s*(Sc_w(1,k)+0.5))-0.5);
    Cs_r(1,k)=(1-theta_b)*cff1*sinh(theta_s*Sc_r(1,k))+...
	      theta_b*(cff2*tanh(theta_s*(Sc_r(1,k)+0.5))-0.5);
  else
    Cs_w(1,k)=Sc_w(1,k);
    Cs_r(1,k)=Sc_r(1,k);
  end;
end;

TheMax=0;
phi=zeros(1,N+1);
phi(1,1)=-1;
for k=1:N
  phi(1,k+1)=Cs_w(1,k);
end;
for k=2:N+1
  TheAlpha=abs(phi(1,k)+phi(k-1))/abs(phi(1,k)-phi(1,k-1));
  TheMax=max(TheMax, TheAlpha);
end;
TheMult=TheMax;

TheMult=0;
for i=1:N-1
  alpha=abs((Cs_w(1,i+1)+Cs_w(1,i))/(Cs_w(1,i+1)-Cs_w(1,i)));
  TheMult=max(TheMult, alpha);
end;

