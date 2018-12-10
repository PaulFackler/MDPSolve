% randdir Generate random Dirichlet vectors
% USAGE
%    x = randdir(n,a);
% INPUTS
%   n  : number of replicates
%   a  : q-vector of positive values
% OUTPUT
%   x  : n x q matrix of random Direchlet values
function x = randdir(n,a)
q = length(a);
x = gamrnd(repmat(a(:)',n,1),1,n,q);
%x = x ./ repmat(sum(x,2),1,q);
x = bsxfun(@rdivide,x,sum(x,2));