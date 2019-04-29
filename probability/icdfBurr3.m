% icdfBurr12 Computes the inverse CDF of the Burr-12 distribution 
%   F(x) = (1+(x/a)^-c)^-b
% USAGE
%  x = icdfBurr12(p,a,b,c);
% INPUTS
%   p   : values on [0,1)
%   a   : scale parameter
%   b,c : shape parameters
% OUTPUT
%   x   : inverse CDF values
function x=icdfBurr12(p,a,b,c)
x = a.*( p.^(-1./b) - 1 ).^(-1./c);
