% factoredBgrid Creates a factored belief grid
% USAGE
%   x=factoredBgrid(q,p);
% INPUTS
%   q : d-vector containing the number of values for each variable
%   p : number of mesh intervals (scalar or d-vector)
% OUTPUT
%   x : sum(q)-column grid of belief points
function x=factoredBgrid(q,p)
d=length(q);
if d==1
  d=length(p);
  q=q+zeros(size(p));
elseif length(p)==1
  p=p+zeros(size(q));
elseif length(p)~=d
  error('q and p must be the same size or scalar')
end

x=simplexgrid(q(1),p(1),1);
for i=2:d
  x=rectgrid(x,simplexgrid(q(i),p(i),1));
end