% S2B converts from a simplex grid to factored belief grid
%  (with possible information loss)
function x=S2B(S,q,p)

d=length(q);
if d==1
  d=length(p);
  q=q+zeros(size(p));
elseif length(p)==1
  p=p+zeros(size(q));
elseif length(p)~=d
  error('q and p must be the same size or scalar')
end

if size(S,2)~=prod(q)
  error(['S must have ' num2str(prod(q)) ' columns'])
end

k=0;
x=zeros(size(S,1),sum(q));
k=0;
for i=1:d
  x(:,k+1:k+q(i))=S*marginalOp(q,i);
  k=k+q(i);
end
end

function M=marginalOp(q,i)
  M=kron(kron(ones(sum(q(1:i-1)),1),eye(q(i))),ones(sum(q(i+1:end)),1));
end