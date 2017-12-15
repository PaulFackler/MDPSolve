% mergecpts Merges the CPTs for a set of variables
% USAGE
%   [P1,parents1]=mergecpts(P0,X,parents0);
% INPUTS
%   P0       : p-element cell array of CPT matrices for individual variables
%   X        : nx x dx matrix of values of conditioning variables
%   parents0 : p-element cell array of parent vectors (these must be row vectors)
% OUTPUTS
%   P1       : CPT matrix for combined variables
%   parents1 : index vector for the parents of the combined variables
%
% The parent vectors contain a set of unique integers on {1,...,dx} that
%   indicate the columns of X associated with the conditioning variables for
%   each variable being combined. These indices must correspond to the way that
%   the conditioning variables are combined in the associated CPT (hence 
%   they need not be in  sorted order).
% P1 and parents1 are always output in sorted order.
%
% This procedure expands each CPT to the full output size and uses
% kroncol to combine CPTs

function [P1,parents1]=mergecpts(P0,X,parents0,options)
if nargin<4, options=[]; end
parents1=unique([parents0{:}]);
if isempty(X)
  if isempty(options)
    error('an options structure must be passed is X is empty')
  end
  if isfield(options,'nx')
    d=length(options.nx);
    options.nx=options.nx(parents1);
  elseif isfield(options,'maxX')
    d=length(options.maxX);
    options.maxX=options.maxX(parents1);
  else
    error('cannot determine the size of X')
  end
else
  d=size(X,2);
end

if length(parents1)==d
  X1=X;
elseif ~isempty(X)
  X1=unique(X(:,parents1),'rows');
else
  X1=[];
end
% create function to list the columns of X1 associated with a parent vector
var=zeros(1,d);  var(parents1)=1:length(parents1); 
pind=@(i) var(parents0{i});
p=length(P0);
ind=getI(X1,pind(1),options);
P1=P0{1}(:,ind);
for i=2:p
  ind=getI(X1,pind(i),options);
  P1=kroncol(P1,P0{i}(:,ind));
end