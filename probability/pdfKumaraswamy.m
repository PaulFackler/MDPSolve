% pdfKumaraswamy Computes the PDF of the Kumaraswamy distribution 
% USAGE
%   p = pdfKumaraswamy(x,a,b);
% INPUTS
%   x   : values of the random variable
%   a,b : parameters
% OUTPUT
%   p   : PDF values
function p = pdfKumaraswamy(x,a,b)
p = a*b*x.^(a-1).*(1-x.^a).^(b-1);