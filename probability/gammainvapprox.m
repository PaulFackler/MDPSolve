% gammainvapprox Computes Gamma variates using approximate inverse CDF
% USAGE
%   [x,v]=gammainvapprox(u,a,S);
% or
%   x=gammainvapprox(u,v);
% INPUTS
%   u   : n-vector of uniformly distributed random values
%   a   : scalar value of the Gamma distribution shape parameter
%   s   : scalar value of the Gamma distribution scale parameter [default: 1]
%   v   : interpolating function information
% OUTPUTS
%   x   : n-vector of Gamma(a) random values
%   v   : interpolating function information
%
% For repeated use with alternative u vectors use
%   [x,v]=gammainvapprox(u,a,s);
% the first time calling and
%   x=gammainvapprox(u,v);
% on subsequent calls (this eliminates the setup on subsequent calls)
%
% If a 2-parameter Gamma distribution is used, multiply x by the scale parameter
%
% This function uses Hermite interpolation to approximate the inverse CDF
% USES: hermiteinv, hermiteinterp
function [x,v]=gammainvapprox(u,a,s)
if ~isscalar(a);
  v=a;
else
  % problem setup
  if nargin<3, s=1; end
  lna=gammaln(a);
  f =@(x) gammainc(x/s,a);
  df=@(x) exp((a-1).*log(x/s) - x - lna)/s;
  
  hi=2*a;
  while f(hi)<1-2*eps, hi=hi*2; end

  % get approximation to the inverse CDF
  opts=struct('tol',1e-12,'mindy',1e-12,'chunk',5000);
  v=hermiteinv(f,df,[0;hi],opts);
end
% get approximate values of x
if isempty(u)
  x=[];
else
  x=hermiteinterp(u,v);
end