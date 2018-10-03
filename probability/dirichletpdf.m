% dirichletpdf Probability density for the Dirichlet distribution
% USAGE 
%   p=dirichletmom(x,a);
% INPUT
%   x : n x d matrix of non-negative values with row sums equal to 1
%   a : d-vector of parameters
% OUTPUT
%   p : n x 1 vector of probability values
function p=dirichletpdf(x,a)
p=gammaln(sum(a))-sum(gammaln(a))+log(x)*(a(:)-1);

