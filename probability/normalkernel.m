% normalkernel Normal (Gaussian) density function without normalizing constant
% USAGE
%   f=normalkernel(x,mu,sigma2);
% INPUT
%   x      : an n-vector of values of a normal random variable
%   mu     : mean
%   sigma2 : variance
% OUTPUT
%   f : an n-vector of non-normalized density values
function f=normalkernel(x,mu,sigma2)
if nargin<=1
  f=exp(-0.5*x.*x);
else
  z=(x-mu)/sqrt(sigma2);
  f=exp(-0.5*z.*z);
end