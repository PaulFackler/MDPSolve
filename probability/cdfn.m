% cdfn Computes the CDF of the standard normal distribution
% USAGE
%   p=cdfn(x);

function p=cdfn(x)
p = 0.5 * erfc(-0.7071067811865475*x);   % x/-sqrt(2)
