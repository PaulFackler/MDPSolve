% bbpdf Probability distribution of the beta-binomial distribution
% USAGE
%   p=bbpdf(x,n,a,b);
% or
%   p=bbpdf(x,n,mu,theta,1);
% INPUTS
%   x   : a vector of integer values between 0 and n
%   n   : a scalar positive integer
%   a,b : positive parameters of the beta distribution
%   mu, theta : alternative parameterization with mu=a/(a+b) and theta=1/(a+b)
%                 mu is the expected value of a single trial
% OUTPUT
%   p   : a vector of probability values
%
% The beta-binomial distribution applies to a binomial random variable for which
%   the success probability has a Beta(a,b) distribution.
% Note that the entire distribution can be obtained be setting x=0:n;
function p=bbpdf(x,n,a,b,flag)
if nargin<5, flag=0; end
if flag
  mu=a; theta=b;
  if theta==0
    p=gammaln(n+1) + x.*log(mu) + (n-x)*log(1-mu) - (gammaln(x+1)+gammaln(n-x+1));
    p=exp(p);
    return
  else
    a=mu/theta; b=(1-mu)/theta;
  end
end
  p=(gammaln(n+1)+gammaln(a+b))-(gammaln(a+b+n)+gammaln(a)+gammaln(b))+ ...
    (gammaln(a+x)+gammaln(b+n-x)) - (gammaln(x+1)+gammaln(n-x+1));
  p=exp(p);
  