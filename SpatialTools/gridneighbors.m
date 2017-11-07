% gridneighbors Creates neighbor pairs for a rectangular grid
% USAGE
%   [net,G,A]=gridneighbors(m,n);
% INPUTS
%   m,n  : grid size mn nodes on an m x n grid
% OUTPUT :
%   net  : 2mn-m-n x 2 matrix of neighbor pairs
%   G    : m x n matrix of node indices
%   A    : mn x mn adjacency matrix
%
% Example:
% [net,G,A]=gridneighbors(3,2)
% net =
%      1     2
%      1     4
%      2     3
%      2     5
%      3     6
%      4     5
%      5     6
% G =
%      1     4
%      2     5
%      3     6

function [net,G,A]=gridneighbors(m,n)
if     m==1, net=[(1:n-1)' (2:n)'];
elseif n==1, net=[(1:m-1)' (2:m)'];
else
  G=reshape(1:m*n,m,n);
  net=[];
  for j=1:n-1
    net=[net; G(1:m-1,j) G(2:m,j); G(1:m,j) G(1:m,j+1)];
  end
  net=[net; G(1:m-1,end) G(2:m,end)];
end
net=sortrows(net);

if nargout>2
  A=sparse([net(:,1);net(:,2)],[net(:,2);net(:,1)],true,m*n,m*n);
end