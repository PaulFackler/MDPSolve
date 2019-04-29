% Burr12mom2param Computes the Burr-12 parameters fromthe mean, CV and skewness  
% USAGE
%   [a,b,c]=Burr3mom2param(mu,cv,gamma);
%
% This is the inverse of Burr12mom
% see: Burr12mom
function [a,b,c]=Burr12mom2param(mu,cv,gamma)
if nargin<3, gamma = 0; end
eta=log( [(1 + cv.^2);  (gamma.*cv.^3 + 3*cv.^2 + 1)] );
bc = log([1/cv;1/cv]); 
bc = log([5;5]); 
usenewton = true;   % uses CompEcon procedure implementing Newton's method
if usenewton
  [bc,~,~,it] = newton(@resfunc,bc,eta); 
  bc1=exp(bc);
  if any(isinf(bc1) | isnan(bc1))
    [a,b,c]=Burr12mom2param(mu,cv/2,2*gamma);
    [bc,~,~,it] = newton(@resfunc,log(bc1),eta);
  end
else
  res0 = inf;
  tol = 1e-7;
  for it=1:100
    [res,J]=resfunc(bc,eta);
    if norm(res,inf)<tol || norm(res./res0-1)<tol
      break
    end
    bc = bc - J\res;
  end
end
%fprintf('iteration count: %1.0f\n',it)
bc=exp(bc);
if bc(1)*bc(2)<=3
  error('could not find parameter values')
end
ex=bc(1)*beta( 1+1/bc(2), bc(1)-1/bc(2) );
a=mu/ex;
if nargout<=1
  a=[a;bc(1);bc(2)];
else
  b=bc(1);
  c=bc(2);
end

function [res,J]=resfunc(bc,eta)
% solve for log of parameters to avoid negative value issues 
logb=bc(1);
bc=exp(bc);
b=bc(1);
c=bc(2);
%b=ifthenelse(b<0,rand(1,1),b);
c=ifthenelse(b*c<=3, (3+5*rand(1,1))/b, c);
%c=max(c, 3+5*rand(1,1));

gb1  = gammaln(b+1);
gc1  = gammaln(1+1./c);
gbc1 = gammaln(b-1./c);
res = [gammaln(b-2./c) + gammaln(1+2./c) +    gb1 - logb  - 2*(gbc1 + gc1);
       gammaln(b-3./c) + gammaln(1+3./c) + 2*(gb1 - logb) - 3*(gbc1 + gc1)] - eta;

if nargout>1
  pb1  = psi(b+1);
  pc1  = psi(1+1./c);
  pbc1 = psi(b-1./c);
  pc2  = psi(1+2./c);
  pbc2 = psi(b-2./c);
  pc3  = psi(1+3./c);
  pbc3 = psi(b-3./c);    
  J = [(pbc2 +   pb1 - 2*pbc1)*b-1  2*((pbc2 - pc2) - (pbc1 - pc1))/c;
       (pbc3 + 2*pb1 - 3*pbc1)*b-2  3*((pbc3 - pc3) - (pbc1 - pc1))/c];
   
end
        