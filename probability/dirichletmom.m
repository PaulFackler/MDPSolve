% dirichletmom First 2 moments of the Dirichlet distribution
% USAGE 
%   [mean,covariance,correlation]=dirichletmom(a);
% INPUT
%   a : d-vector of parameters
% OUTPUTS
%   mean         : d-vector of mean values
%   covariance   : d x d covariance matrix
%   correlation  : d x d correlation matrix
%
% The domain of the Dirichlet distribution is a d-vector of 
% non-negative values that sum to 1. Note that this implies that
% the covariance and correlation matrices are rank d-1
function [mean,covariance,correlation]=dirichletmom(a)
mean = a/sum(a);
if nargout>1
  variance = a.*(sum(a)-a)/(sum(a).^2*(sum(a)+1));
  C = -a(:)*a(:)'/(sum(a).^2*(sum(a)+1));
  covariance=C-diag(diag(C))+diag(variance);
  if nargout>2
    correlation=covariance./(sqrt(variance)'*sqrt(variance));
  end
end
