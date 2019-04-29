% randbeta generates a vector of ranom Beta variates
% USAGE
%   x=randbeta(n,a,b);
function x=randbeta(n,a,b)
x=betaincinv(rand(n,1),a,b);