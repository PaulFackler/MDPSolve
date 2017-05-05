% gammakernel Gamma density function without normalizing constant
% USAGE
%   f=gammakernel(x,lambda,theta);
% INPUTS
%   x : an n-vector of values of a Gamma random variable
%   lambda : a scalar shape parameter 
%   theta  : a scalar scale parameter
% OUTPUT
%   f      : an n-vector of non-normalized density values
% 
% f(x) = (x/theta)^(lambda-1)*exp(-x/theta)
% E[x]   = lambda*theta
% Var[x] = lambda*theta^2
%
% values of f for x<0 or x=NaN are set to 0
function f=gammakernel(x,lambda,theta) 
x(isnan(x))=0;
x(x==0)=realmin;
f=x/theta;
f=exp((lambda-1)*log(f) - f );
f(x<0)=0;