% invpsi Inverse of the digamma (psi) function (uses Newton's method) 
% USAGE 
%   y=invpsi(x);
% INPUT
%   x : n-vector of numbers
% OUTPUT
%   y : n-vector of value of the inverse of the psi function
%
% The digamma (psi) function is the derivative of the log gamma function
function y=invpsi(x)
tol=5e-15;
maxit=100;
% initial estimate (the maximal error occurs at -1.6605 and equals 0.1216)
y = ifthenelse(x >= -1.6605, exp(x) + 0.5, -1./(x+0.5772));
for i=1:maxit
  r=psi(y)-x;
  if norm(r,inf)<tol
    break;
  end
  y = y - r./psi(1,y);
end
%disp(i)
return   
 