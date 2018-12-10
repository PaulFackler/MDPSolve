% EVoptgroups Determines the optimal way to group factors in an EV function
% USAGE
%   [mergevec,c]=EVoptgroups(p,parents,X,e,options);
% INPUTS
%   p       : m-element cell array of conditional probability matrices
%   parents : m-element cell array of conditioning variable (parent) indices
%   X       : matrix or cell array of state/action variables
%   e       : cell array of rv structures (discrete or w/ discrete approximations)
%   options : options structure
%               penalty - penalizes more factors
%               indexed - use the operation count for an indexed evaluation
% OUTPUTS
%    mergevec : a vector of length <= d indicting the number of factors in each group
%                 note: sizes sums to d
%    c        : cost of optimal grouping
%
% Example: for a 6 variable problem mergevec=[3 1 2] indicates that the 
%          first 3 factors should be merged as should the last 2
function [mergevec,c]=EVoptgroups(p,parents,X,e,options)
penalty=0;
indexed=false;
if exist('options','var') && ~isempty(options)
  if isfield(options,'penalty'),   penalty = options.penalty; end
  if isfield(options,'indexed'),   indexed = options.indexed; end
end

% get row sizes
N=fliplr(cumprod(fliplr(cellfun(@(x)size(x,1),p))));
%get column sizes of X variables
xparents=cellfun(@(x)x(x>0),parents,'UniformOutput',false);
if iscell(X)
  ds = length(parents);
  nx = cellfun(@length,X);
  m=zeros(1,ds);
  common = xparents{1};
  m(1)=prod(nx(common));
  for i=2:ds
    common=union(common, xparents{i});
    m(i)=prod(nx(common));
  end
else
  [Ip,Iy] = EVindices(xparents,X); 
  m=cellfun(@(x)numel(x),Ip);
end
  
if nargin<4 || isempty(e)
  de = 0;
else
  if ~iscell(e), e={e}; end
  de = length(e);
  [first,last] = getenterexit(parents,e);
  ne = cellfun(@(x)length(x.values),e);
end
  
if any(diff(m)<0)
  error('m must be non-decreasing')
end
if any(diff(N)>=0)
  error('p must be increasing')
end
if indexed, m=min(m,N(1)); end
ds=length(N);
C=triu(N'*m);
Q=cell(ds,ds);
for i = 1:ds, 
  for j = i:ds, 
    for k = 1:de
      if i>first(k) || j<last(k)
        C(i,j) = C(i,j)*ne(k);  % increase by size of included shocks
      end
    end
    Q{i,j} = {[i j]};  % used for storing optimal splits
  end; 
end

for j = 2:ds
  for i = j-1:-1:1
    for k = 1:j-i
      temp = ( C(i,j-k) + C(j-k+1,j) )*(1+penalty);
      if temp < C(i,j)
        C(i,j) = temp;
        Q{i,j} = [ Q{i,j-k}  Q{j-k+1,j} ];
      end
    end
  end
end

mergevec = Q{1,ds};
mergevec = reshape([mergevec{:}],2,length(mergevec));
mergevec = mergevec(2,:)-mergevec(1,:)+1;
c=C(1,ds);
