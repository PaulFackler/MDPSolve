% betakernel Beta density function without normalizing constant
% USAGE
%   f=betakernel(x,a,b);
% INPUTS
%   x   : an n-vector of values of a Beta random variable
%   a,b : scalar Beta distribution parameters
% OUTPUT
%   f      : an n-vector of non-normalized density values
% 
% f(x) = x^(a-1)*(1-x)^(b-1)
% E[x]   = a/(a+b)
% Var[x] = ab/(a+b)^2/(a+b+1)
%
% values of f for x<0, x>1 or x=NaN are set to 0
function f=betakernel(x,a,b) 
x(isnan(x))=0;
x(x==0)=realmin;
x(x==1)=1-eps;
f=(a-1)*log(x) + (b-1)*log(1-x);
f=exp( f );
f(x<0 | x>1)=0;