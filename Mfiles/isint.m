%
% ISINT determine if input matrix is integers
%
% ISINT a Matlab function to determine which input matrix positions 
%       contain an integer or one that can be represented as an integer.
%       ISINT returns 1 where matrix entry is integer and a 0 otherwise.
%
%       If x is a vector or matrix, ISINT returns a vector or matrix 
%       the same size as x with each element of the input being 
%       evaluated as above.
%
%       Ex. 1;  if x=3.2, ISINT returns 0
%               if x=3.0, ISINT returns 1
%               if x=3, ISINT returns 1
%               if x=0, ISINT returns 1
%
%       Ex. 2;  >> x=[3   3.5
%                     4.5 6  ];
%               >> isint(x)
%                 
%               ans = 1 0 
%                     0 1
%
%
% Call as: >> retval=isint(x)
%
% Written by : Brian O. Blanton
%
function retval=isint(x)
%ints=x-floor(x)==0;
%retval=ints.*x;
retval=x-floor(x)==0;

%
%        Brian O. Blanton
%        Curr. in Marine Science
%        15-1A Venable Hall
%        CB# 3300
%        Uni. of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        919-962-4466
%        blanton@marine.unc.edu
%
