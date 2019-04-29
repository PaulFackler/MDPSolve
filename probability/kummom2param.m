% kummom2param Computes the Kumaraswamy parameters from the mean and CV
% USAGE
%   [a,b] = kummom2param(mu,cv);

function [a,b]=kummom2param(mu,cv)
c=log([mu; mu.^2.*(1+cv.^2)]);

% not clear how to pick good starting values
ab0=[mu*cv; cv/mu]; 
ab0=[1;1];
%ab0=[1/(mu*cv);1/(mu*cv)];
%kapprox = kumarapprox;
%ab0=kapprox(mu,cv)';
if 0
  ab=log(ab0);
  res0=inf;
  for it=1:500
    [res,J]=resfunc(ab,c);
    ab = ab - J\res;
    if norm(res./res0-1,inf)<1e-8
      break
    end
    res0 = res;
  end
else
  maxit = optget('newton','maxit');
  optset('newton','maxit',1000)
  optset('newton','tol',5e-7)
  optset('newton','showiters',0);
  [ab,~,~,it] = newton(@resfunc,ab0,c); 
  optset('newton','maxit', maxit)
  optset('newton','tol',1e-8)
  optset('newton','showiters',0);
end
disp(it)
ab=exp(ab);
if nargout==1
  a=ab;
else
  a = ab(1);
  b = ab(2);
end

function [res,J]=resfunc(ab,c)
% solve for log of parameters to avoid negative value issues 
lb=ab(2);
ab=exp(ab);
a=ab(1);
b=ab(2);
glb = gammaln(b);
res = [lb + gammaln(1+1./a) + glb - gammaln(1+1./a+b);
       lb + gammaln(1+2./a) + glb - gammaln(1+2./a+b)] - c;
if nargout>1
  pb = psi(b);
  pab1 = psi(1+1./a+b);
  pab2 = psi(1+2./a+b);
  J = [  (pab1-psi(1+1./a))/a  1+(pb-pab1)*b;
       2*(pab2-psi(1+2./a))/a  1+(pb-pab2)*b];
end
        