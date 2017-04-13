% KLdiv Kullback-Leibler divergence
% USAGE
%   KL=KLdiv(p,q);
% INPUTS
%   p,q : nxm matrices with non-negative elements and rows that sum to 1
% OUTPUT
%  KL   : n-vector of values the the KL divergence from p to q
function KL=KLdiv(p,q)
ind=p>0 | q>0;
pind=p(ind);
KL = sum(pind.*log2((pind+realmin)./q(ind)),2);