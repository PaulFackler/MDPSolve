% qnwrvs Creates quadrature nodes and weights for a set of rv structures
% USAGE
%   [e,w]=qnwrvs(rvs);
% INPUT
%   rvs  a p-element cell array of rv strctures (defined by rvdef)
% OUTPUTS
%   e    an nxp matrix of nodal values
%   w    an nx1 vector of probability weights
function [e,w]=qnwrvs(rvs)

e=cellfun(@(x)x.values,rvs,'UniformOutput',0);
w=cellfun(@(x)x.cpt,rvs,'UniformOutput',0);

e=rectgrid(e);
w=prod(rectgrid(w),2);
