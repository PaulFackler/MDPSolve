% randepanechnikov Generates random variables from the Epanechnikov density
% USAGE
%   x = randepanechnikov(n);
% Returns an n x 1 vector of random values. 
% Uses Newton's method to implement the inverse CDF method.
% The Epanechnikov density is f(x) = 3/4(1-x^2) for x in [-1,1]
% The CDF is F(x) = (3*x-x.^3+2)/4
function x = randepanechnikov(n)
c=2*rand(n,1)-1;
x=c;
d=2/3;
for k=1:16
  x = d*(x+(x-c)./(x.^2-1));
end
