function x=icdfraisedcos(u)
maxit=28;
c=(2*u-1)*pi;
x=zeros(size(u));
for i=1:maxit
  res=c-(x+sin(x));
  x=x+res./(1+cos(x));
end
x=x/pi;