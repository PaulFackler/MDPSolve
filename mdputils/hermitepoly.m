% Hermitepoly Constructs a complete set of Hermite polynomials up to order k
% USAGE
%   Z=hermitepoly(X,k,m,s);
% INPUTS
%   X : n x d matrix
%   k : integer scalar
%   m : 1 x d vector of variable means
%   s : 1 x d vector of variable standard deviations
% OUTPUT
%   Z : n x p matrix of polynomial functions of X w/ p=(k+d)!/(k!d!)
function Z=hermitepoly(X,k,m,s)
if nargin>=4
  Z=bsxfun(@rdivide,bsxfun(@minus,X,m(:)'),s(:)');
else
  Z = X;
end
Z=crossorthmom(Z,k);