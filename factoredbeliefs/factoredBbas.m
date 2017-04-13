% factoredBgrid Creates an interpolation matrix for a factored belief space
% USAGE
%   B=factoredBbas(x,q,p);
% INPUTS
%   x : m x sum(q) matrix of evaluation points 
%   q : d-vector containing the number of values for each variable
%   p : number of mesh intervals (scalar or d-vector)
% OUTPUT
%   B : sum(q) x m interpolation matrix
%
% see: factoredBgrid
function B=factoredBbas(x,q,p)
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
B=simplexbas(x(:,k+1:k+q(1)),q(1),p(1),1);
k=k+q(1);
for i=2:d
  B=kroncol(B,simplexbas(x(:,k+1:k+q(i)),q(i),p(i),1));
  k=k+q(i);
end