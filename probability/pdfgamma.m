function f=pdfgamma(x,a,lnG)
if nargin<3 || isempty(lnG), lnG=gammaln(a); end
f=exp((a-1)*log(x)-x-lnG);
f(a<=0)=NaN;
f(x<0)=0;
f(x==0 & a>1)=0;
f(x==0 & a<1)=inf;
f(x==0 & a==1)=1;