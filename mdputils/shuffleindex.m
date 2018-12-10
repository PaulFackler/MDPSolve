% shuffleindex Creates an index vector that reorders a nd-array
% USAGE 
% ind=shuffleindex(n,order);
% INPUTS
%   n             : m-vector of dimension sizes
%   order         : m-vector permutation of 1:m
%   lexicographic : 0/1 set to 0 if reverse lexicographic [default: 1]
% OUTPUT
%   ind           : prod(n) vector that reshuffles an nd-array
%                     returns [] if order=1:m (no reordering)
% 
% Suppose that X is an m-dimensional array with n=size(X) (reverse lexicographic order)
% Let 
%   Y=permute(X,order);
% Y can also be obtained using 
%   ind=shuffleindex(n,order,0); Y=reshape(X(ind),n(order));
%
% Suppose that X is an m-dimensional array with n=fliplr(size(X)) (lexicographic order)
% Let 
%   Y=permute(X,m+1-fliplr(order));
% Y can also be obtained using 
%   ind=shuffleindex(n,order,1); Y=reshape(X(ind),fliplr(n(order))); 
%
% Example: rearrange columns of a lexicographically ordered matrix and then
%   use index vector from shuffleindex to rearrange the rows so they are
%   lexicographically ordered
% n=[2 3 4]; 
% order=[2 1 3]; 
% X=rectgrid(arrayfun(@(x)(1:x)',n,'UniformOutput', false)); 
% ind1=shuffleindex(n,order); 
% X(ind1,order)

function ind=shuffleindex(n,order,lexicographic)
  if nargin<3 || isempty(lexicographic), lexicographic=true; end
  m=length(n);
  if length(order)~=m
    error('inputs must be the same length')
  end
  if isequal(order(:)',1:m)
    ind=[]; return
  end
  if lexicographic
    n=fliplr(n);
    order=m+1-fliplr(order);
  end
  ind=reshape(1:prod(n),n);
  ind=permute(ind,order);
  ind=ind(:);