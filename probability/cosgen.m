% cosgen Random number generator for Cosine distribution
% USAGE
%   x=cosgen(n,rv,parents,z);
% INPUTS
%   n       : number of variates to generate
%   rv      : rv structure
%   parents : not used
%   z       : n-vector of random uniform values
% OUTPUT
%   x       : n-vector of random Cosine distributed variates
function x=cosgen(n,rv,parents,z)
if nargin<4 || isempty(z)
  z=rand(n,1);
end
ab=rv.parameters;
if isnumeric(ab)
  x=asin(2*z-1).*((2/pi)*rv.ab(2))+ab(1);
else
  ab=(ab(parents{:}))';
  x=asin(2*z-1).*((2/pi)*ab(:,2))+ab(:,1);
end