% CKRONX The product of repeated Kronecker products and a matrix.
% USAGE      
%   z=ckronx(B,c,ind,transpose);
% INPUTS
%   B   : a d-element cell array with element i an m(i) x n(i) matrix
%   c   : a compatible matrix prod(n) x p (prod(m) x p if transpose=1)
%   ind : a selection vector of numbers on 1..d that selects elements of B
%            [default: 1:d]
%   transpose: 1 to use the transpose of the elements of B [default: 0]
% OUTPUT  
%   z :  prod(m) x p matrix (prod(n) x p if transpose=1)
% Solves (B1xB2x...xBd)*c
% where x denotes Kronecker (tensor) product.
% The Bi are passed as a cell array B. 
% B must be a vector cell array containing 2-D numerical arrays.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function z=ckronx(b,c,ind,transpose)

if nargin<2, error('At least two parameters must be passed'), end
if nargin<3 || isempty(ind), ind=1:numel(b); end
if nargin<4 || isempty(transpose), transpose=false; end

if ~iscell(b)                         % b is a matrix: return b*c
  if size(b,2)~=size(c,1)
    error('b and c are not conformable')
  end
  z=b*c;
else                                  % b is a cell array
  d=length(ind);
  n=zeros(d,1);
  if transpose
    for i=1:d, n(i)=size(b{ind(i)},1); end
  else
    for i=1:d, n(i)=size(b{ind(i)},2); end
  end
  if prod(n)~=size(c,1)
    error('b and c are not conformable')
  end
  z=c';
  mm=1;
  if transpose
    for i=1:d
      m=numel(z)/n(i);
      z=reshape(z,m,n(i))';
      z=b{ind(i)}'*z;
      mm=mm*size(z,1);
    end
  else
    for i=1:d
      m=numel(z)/n(i);
      z=reshape(z,m,n(i))';
      z=b{ind(i)}*z;
      mm=mm*size(z,1);
    end
  end
  z=reshape(z,mm,size(c,2));
end
