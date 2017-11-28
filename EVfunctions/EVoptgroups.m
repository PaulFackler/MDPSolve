% EVoptgroups Determines the optimal way to group factors in an EV function
% USAGE
%   [order,c]=EVoptgroups(p,m,options);
% INPUTS
%   p        : d-vector of cumulative factor sizes from i-d
%   m       : d-vector of sizes of common conditioning variables 1-i
%   options : options structure
%               penalty - penalizes more factors
%               indexed - use the operation count for an indexed evaluation
% OUTPUTS
%    order : a vector of length <= d indicting the number of factors in each group
%              note: order sums to d
%    c     : cost of optimal order
%
% Example: for a 6 variable problem order=[3 1 2] indicates that the 
%          first 3 factors should be merged as should the last 2
function [order,c]=EVoptgroups(p,m,options)
penalty=0;
indexed=false;
if exist('options','var') && ~isempty(options)
  if isfield(options,'penalty'),   penalty = options.penalty; end
  if isfield(options,'indexed'),   indexed = options.indexed; end
end
if indexed, m=min(m,p(1)); end
d=length(p);
C=triu(p'*m);
Q=cell(d,d);
for i=1:d, for j=i:d, Q{i,j}={[i j]}; end; end

for j=2:d
  for i=j-1:-1:1
    for k=1:j-i
      temp=(C(i,j-k)+C(j-k+1,j))*(1+penalty);
      if temp < C(i,j)
        C(i,j)=temp;
        Q{i,j}=[Q{i,j-k} Q{j-k+1,j}];
      end
    end
  end
end

order=Q{1,d};
order=reshape([order{:}],2,length(order));
order=order(2,:)-order(1,:)+1;
c=C(1,d);
