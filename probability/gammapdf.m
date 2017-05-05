% gammapdf Gamma density function
% USAGE
%   f=gammapdf(x,lambda,theta);
% INPUTS
%   x      : an n-vector of values of a Gamma random variable
%   lambda : a scalar shape parameter 
%   theta  : a scalar scale parameter
% OUTPUT
%   f      : an n-vector of density values
% 
% f(x) = (x/theta)^(lambda-1)*exp(-x/theta)
% E[x]   = lambda*theta
% Var[x] = lambda*theta^2
%
% values of f for x<0 or x=NaN are set to 0
function f=gammapdf(x,lambda,theta) 
x(isnan(x))=0;
x(x==0)=realmin;
f=x/theta;
f=exp((lambda-1)*log(f) - f - gammaln(lambda))/theta;
f(x<0)=0;