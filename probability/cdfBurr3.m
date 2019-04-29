% cdfBurr3 Computes the CDF of the Burr-3 distribution F(x) = (1+(x/a)^-c)^-b
% USAGE
%   p = cdfBurr3(x,a,b,c);
% INPUTS
%   x   : values of the random variable
%   a   : scale parameter
%   b,c : shape parameters
% OUTPUT
%   p   : CDF values
function p=cdfBurr3(x,a,b,c)
z=(x./a).^c;
p=(z./(1+z)).^b;