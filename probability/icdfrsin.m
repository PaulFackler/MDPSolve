function x=icdfraisedcos(u)
tol=1e-11;
maxit=100;
c=(2*u-1)*pi;
x=zeros(size(u));
for i=1:maxit
  res=x+sin(x)-c;
  x=x-res./(1+cos(x));
  if all(abs(res)<tol)
    break
  end
end