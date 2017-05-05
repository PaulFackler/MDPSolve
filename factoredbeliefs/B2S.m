% B2S converts from a factored belief grid to a simplex grid
function S=B2S(x,q,p)

d=length(q);
if d==1
  d=length(p);
  q=q+zeros(size(p));
elseif length(p)==1
  p=p+zeros(size(q));
elseif length(p)~=d
  error('q and p must be the same size or scalar')
end

if size(x,2)~=sum(q)
  error(['x must have ' num2str(sum(q)) ' columns'])
end

k=0;
S=x(:,k+1:k+q(1));
k=k+q(1);
for i=2:d
  S=kronrow(S,x(:,k+1:k+q(i)));
  k=k+q(i);
end