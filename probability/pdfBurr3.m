% pdfBurr3 Computes the PDF of the Burr-3 distribution F(x) = (1+(x/a)^-c)^-b
% USAGE
%   p = pdfBurr3(x,a,b,c);
% INPUTS
%   x   : values of the random variable
%   a   : scale parameter
%   b,c : shape parameters
% OUTPUT
%   p   : PDF values
function p=pdfBurr3(x,a,b,c)
z=(x./a);
p=(b*c/a)*z.^-(c+1)./(1+z.^-c).^(b+1);
