% mergecpts Merges the CPTs for a set of variables
% USAGE
%   [P,parents]=mergecpts(p,parents,X,e,eremove,options);
% INPUTS
%   p        : m-element cell array of CPT matrices for individual variables
%   parents  : m-element cell array of parent vectors (these must be row vectors)
%   X        : nx x dx matrix of values of conditioning variables
%   e        : de-element cell array of rv structures (discrete or 
%                w/ discrete approximations)
%   eremove  : de-vector specifying which of the noise variables are summed out
%   options  : structure variable (fields described below)
% OUTPUTS
%   P        : CPT matrix for combined variables
%   parentsc : index vector for the parents of the combined variables
%
% The parent vectors contain a set of unique integers on {1,...,dx} that
%   indicate the columns of X associated with the conditioning variables for
%   each variable being combined. These indices must correspond to the way that
%   the conditioning variables are combined in the associated CPT (hence 
%   they need not be in  sorted order).
% The presense of noise variables is indicated be negative values in the
%   parent vectors. Any noise variables that are not summed out will also
%   enter the P matrix in sorted order.
% P and parentsc are always output in sorted order.

function  [P,parentsc] = mergecpts(p,parents,X,e,eremove)
% get information about noise variables (e)
if nargin<4
  e={};
  de=0;
else
  de = length(e);
end
if nargin<5 || isempty(eremove), eremove=false(de,1); end
if length(eremove) == 1, eremove = repmat(eremove,de,1); end
[eenter,eexit] = getenterexit(parents,e);
[Ip,Iy] = EVindices(parents,X,e,eremove);

m=length(p);
if m>1
  parentsc=unique([parents{:}]);
  parentsc(parentsc<0) = sort(parentsc(parentsc<0),'descend');
else
  parentsc=parents{1};
  parentsc=[sort(parentsc(parentsc<0),'descend') parentsc(parentsc>0)];
end
w=cellfun(@(x)x.cpt,e,'UniformOutput',false);
for i=1:m
  if i==1
    if isempty(Ip{1})
      P=p{1};
    else
      P=p{1}(:,Ip{1});
    end
  else
    P = kroncol(P(:,Iy{i}),p{i}(:,Ip{i}));
  end
  % sum out noise variables as appropriate
  if ~isempty(e)
    if any(eremove & eexit==i)
      ii=find(eremove & eexit==i);
      wi = w{ii(1)};
      for j=2:length(ii)
        wi = w{ii(j)}*wi';  
        wi=wi(:);
      end
      P=reshape(reshape(P,[],length(wi))*wi,size(P,1),[]); 
      % remove summed noise terms from parent list
      parentsc(matchindices(-ii,parentsc)) = [];
    end
  end
end
return

