% lognkernel Lognormal density function without normalizing constant
% USAGE
%   f=lognkernel(x,mu,sigma2);
% INPUT
%   x      : an n-vector of values of a log normal random variable
%   mu     : mean of ln(x)
%   sigma2 : variance of ln(x)
% OUTPUT
%   f      : an n-vector of non-normalized density values
function f=lognkernel(x,mu,sigma2)
if nargin>1
  z=(log(x)-mu)./sqrt(sigma2);
else
  z=log(x);
end
f=exp(-0.5*z.*z)./x;
f(x==0)=0;