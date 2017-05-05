% raisedcosgen Random number generator for Raised Cosine distribution
% USAGE
%   x=raisedcosgen(n,rv,parents,z);
% INPUTS
%   n       : number of variates to generate
%   rv      : rv structure
%   parents : not used
%   z       : n-vector of random uniform values
% OUTPUT
%   x       : n-vector of random Raised Cosine variates
function x=raisedcosgen(n,rv,parents,z)
if nargin<4 || isempty(z)
  z=rand(n,1);
end
x=icdfraisedcos(z);
x=x*(rv.parameters(2)/pi)+rv.parameters(1);

% generates a standard (0,pi) variate
function x=icdfraisedcos(u)
maxit=28;
c=(2*u-1)*pi;
x=zeros(size(u));
for i=1:maxit
  res=c-(x+sin(x));
  x=x+res./(1+cos(x));
end