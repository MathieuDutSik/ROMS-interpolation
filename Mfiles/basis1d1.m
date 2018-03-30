% function [b2,n2,b1,n1]=basis1d1(zref,zdes)
%-----------------------------------------------------------------------
% purpose: This subroutine evaluates the value of the basis functions 
%            alive on a 1-D linear element at a point
%
% restrictions: Applicable only for  1-D linear elements
%
% inputs:  nn is the number of entries in the reference z array
%          zref(nn) is the 1-D array containing the reference z array
%               NOTE: zref must increase from zref(1) to zref(nn)
%          zdes is the z value at which the basis functions are desired
%
% outputs: n1 is the reference node immediately below zdes
%          b1 is the value of the n1's basis function at zdes
%          n2 is the reference node immediately above zdes
%          b2 is the value of the n2's basis function at zdes
%
% notes:   zdes<zref(1)  => n1=nn,b1=0.0,n2=1,b2=1.0 
%          zdes>zref(nn) => n1=nn,b1=1.0,n2=1,b2=0.0 
%
% history:  Written by Christopher E. Naimie
%           Dartmouth College
%           26 AUGUST 1992
%-----------------------------------------------------------------------
function [b2,n2,b1,n1]=basis1d1(zref,zdes)
%
zdiff=zref-zdes;
nn=length(zref);
if zdiff(1) >= 0;
   n1=nn;
   b1=0.0;
   n2=1;
   b2=1.0; 
elseif zdiff(length(zdiff)) <= 0;
   n1=nn;
   b1=1.0;
   n2=1;
   b2=0.0;
else
   zpast=find(zdiff<=0);
   n1=zpast(length(zpast));
   n2=n1+1;
   dz=zref(n2)-zref(n1);
   b1=(zref(n2)-zdes)/dz;
   b2=(zdes-zref(n1))/dz;
end
