% matchindices Finds the values of ind1 that match those of ind0
% USAGE
%   [ind,indx]=matchindices(ind0,ind1);
% INPUTS
%   ind0 : n-vector of numbers
%   ind1 : m-vector of numbers
% OUTPUTS
%   ind  : p-vector with ind(i) = place in ind1 that matches a value in ind0
%   indx : q-vector of palces in ind0 with no meching values in ind1
%
% Note that p+q = n 
% Values of ind are in {1,...,m}
% Values of indx are in {1,...,n}
% Example:
%   [ind,indx] = matchindices([1 4 6 2],[4 3 2 1]) 
% yields:
%    ind =
%         4     1     3
%    indx =
%         3
function [ind,indx]=matchindices(ind0,ind1)
  [~,ind]=ismember(ind0,ind1);
  indx=find(ind==0);
  ind=ind(ind>0);