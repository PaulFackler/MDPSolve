% HERMITEINTERP Interpolates using piecewise cubic Hermite interpolation
% USAGE
%    [yi,dyi]=hermiteinterp(xi,v);
%   or
%    [yi,dyi]=hermiteinterp(xi,x,y,dy);
% INPUTS
%  xi       : an m-vector of interpolation points
%  v        : an 3-element cell array containing p-vectors
%                x  : values of input variable (sorted and unique values)
%                y  : values of output variable
%                dy : values of derivative
% OUTPUTS
%  yi       : m-vector of interpolated function values
%  dyi      : m-vector of interpolated derivative values
%
% Based on:
% Wolfgang Hörmann and Josef Leydold (2003)
% "Continuous Random Variate Generation by Fast Numerical Inversion"
% ACM Transactions on Modeling and Computer Simulation, 13: 347-362.
%
% Note: v can be obtained from hermiteinv
function [y,dy]=hermiteinterp(x,xi,dxi,yi,a1,a2,a3)
if nargin<3
  v=xi;
  xi=v{1};
  dxi=v{2};
  yi=v{3};
  a1=v{4};
  a2=v{5};
  a3=v{6};
end


% make endpoints -inf and inf to get index values
xx=xi; xx(1)=-inf; xx(end)=inf;  
[~,ind]=histc(x,xx);

lambda = (x-xi(ind))./dxi(ind);

if nargout>1
  a1=a1(ind);
  a2=a2(ind);
  a3=a3(ind);
  y=((a3.*lambda+a2).*lambda+a1).*lambda+yi(ind);
  dy=((3*a3.*lambda+2*a2).*lambda+a1)./dxi(ind);
else
  y=((a3(ind).*lambda+a2(ind)).*lambda+a1(ind)).*lambda+yi(ind);
end
