function [a,b]=gamma2Burr(alpha,beta)
if nargin==1
  beta=alpha(2);
  alpha=alpha(1);
end
[x,w]=qnwgamma(101,alpha,beta);
c=w'*log(x);
f=@(a) -score(a,x,w,c);
a=fminbnd(f,alpha/2,alpha*2);
b=1./(w'*log(1+x.^a));
if nargout<=1
  a=[a;b];
end
return

function s=score(a,x,w,c)
z=w'*log(1+x.^a);
b=1./z;
s=log(a*b)+(a-1)*c-(b+1)*z;