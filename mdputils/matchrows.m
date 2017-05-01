% matchrows Finds the rows in B that match those in A
% USAGE
%   [I,Q]=matchrows(A,B);
% INPUTS
%   A  : m x n matrix
%   B  : p x n matrix
% OUTPUTS
%   I  : m x 1 vector given the row numbers of B corresponding to each row in A
%   Q  : m x p selection matrix such that A=Q'*B
%
% Note: if A has rows that do not match any row in B, I is set equal to 0 and
%   Q=0 in the row; a warning is issued.

function [I,Q]=matchrows(A,B)
[~,I]=ismember(A,B,'rows'); 
ii=find(I>0);
if length(ii)~=size(A,1)
  s=sprintf('A contains %1.0f row(s) out of %1.0f that do not match any rows in B\n',...
    size(A,1)-length(ii),size(A,1));
  warning(s)
end
if nargout>1
  Q=sparse(I(ii),ii,1,size(B,1),size(A,1));
end