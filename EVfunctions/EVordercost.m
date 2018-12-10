% EVordercost Cost of using a particular order for an EV function
% USAGE
%   w=EVordercost(I,n,q,mx,me);
% INPUTS
%   I  : d permuation vector of values {1,...,d}
%   n  : d-vector of state variable sizes
%   q  : d-element cell array with parent vectors for each state variable
%   mx : vector of conditioning variable sizes
%   me : vector of noise variable sizes
% OUTPUT
%   cost : number of flops used to compute the EV function
function cost=EVordercost(I,n,q,mx,me)
% not implemented yet for noise terms
if nargin < 5, me = []; end
d=length(I);
q=q(I);
n=n(I);
combined=q{1};
N=prod(n);
cost=N*prod(mx(combined));
for i=2:d
  combined=union(combined,q{i});
  N = N/n(i-1);
  cost = cost + N*prod(mx(combined));
end