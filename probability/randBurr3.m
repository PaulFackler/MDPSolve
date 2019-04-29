% randBurr3 Generates random values from the Burr-3 distribution 
%             F(x) = (1+(x/a)^-c)^-b
% USAGE
%   x = randBurr3(n,a,b,c);
% INPUTS
%   n   : number of values to generate 
%   a   : scale parameter
%   b,c : shape parameters
% OUTPUT
%   x   : n-vector of random values
function x=randBurr3(n,a,b,c)
x = a.*( rand(n,1).^(-1./b) - 1 ).^(-1./c);


%x = rand(n,1).^(1./b);
%x = a.*(x./(1-x)).^(1./c);