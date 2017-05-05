function P=ProbMofN(p,n)
z=p==0;
p(z)=[];
if nargin>1
  n(z)=[];
  p=p(n);
end
alpha=prod(p);
s=-(1-p)./p;
S=poly(s);
P=alpha*S(:);
