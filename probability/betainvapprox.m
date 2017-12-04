% betainvapprox Computes Beta variates using approximate inverse CDF
% USAGE
%   [x,v]=betainvapprox(u,a,b);
% or
%   x=betainvapprox(u,v);
% INPUTS
%   u   : n-vector of uniformly distributed random values
%   a,b : scalar values of the Beta distribution parameters
%   v   : interpolating function information
% OUTPUTS
%   x   : n-vector of Beta(a,b) random values
%   v   : interpolating function information
%
% For repeated use with alternative u vectors use
%   [x,v]=betainvapprox(u,a,b);
% the first time calling and
%   x=betainvapprox(u,v);
% on subsequent calls (this eliminates the setup on subsequent calls)
%
% This function uses Hermite interpolation to approximate the inverse CDF
% USES: hermiteinv, hermiteinterp
function [x,v]=betainvapprox(u,a,b)
if nargin==2;
  v=a;
else
  % problem setup
  lnB=gammaln(a)+gammaln(b)-gammaln(a+b);
  f =@(x) betainc(x,a,b);
  df=@(x) exp((a-1).*log(x) + (b-1).*log(1-x) - lnB);

  % get approximation to the inverse CDF
  opts=struct('tol',1e-10,'mindy',1e-10,'chunk',5000);
  v=hermiteinv(f,df,[0;0.5;1],opts);
end
% get approximate values of x
if isempty(u)
  x=[];
else
  x=hermiteinterp(u,v{:});
end