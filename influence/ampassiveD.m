% ampassiveD Creates a model structure for passive adaptive management
% USAGE
%   model=ampassiveD(D,w,options);
% INPUTS
%   D       : an influence diagram structure
%   w       : weights on unknown parameters
%   options : an options structure passed to d2model 
%               (see doumentation for d2model)
% OUTPUTS
%   model   : model structure with weights attached to parameters
%
% If there is more than one parameter variable the weights input
%   can be defined as a single probability vector or as a cell array of
%   probability vectors (vectors of non-negative numbers that sum to 1)
% If a single vector is used it must have length equal to the number of
%   parameter value combinations, with the combinations in lexicographic
%   order with variables ordered as they are entered into the diagram.
% If a cell array is used each vector must have the same number of values
%   as the number of values for the associated variable.
%
% Example:
%   Suppose p1 hase values 1, 2 and 3 and p2 has values 0 and 1. The weight
%   variable could be w=[0.1;0.2;0.2;0.1;0.2;0.2] or it could be
%   w={[0.2;0.4;0.4],[0.5;0.5]}; these are equivalent. The second approach
%   allows specification of a prior in which the joint does not equal the
%   product of the marginals.
function model=ampassiveD(D,w,options)
if nargin<3, options=[]; end
ii=find(ismember(D.types,'p'));
if length(ii)==1 
  D.types{ii}='c';
  D.cpds{ii}=rvdef('d',w,D.values{ii});
elseif iscell(w)
  if length(ii) ~= length(w)
    error('if w is a cell array it must have the same number of elements as the number of parameters')
  end
  for i=1:length(ii)
     k=ii(i);
     if length(w{i})~=length(D.values{k})
       error(['weights for ' D.names{k} ' not the same length as the values'])
     end
     D.types{k}='c';
     D.cpds{k}=rvdef('d',w{i},D.values{k});
  end
else
  pvals=dvalues(D,ii,'m');
  for i=1:length(ii)
     k=ii(i);
     D.types{k}='c';
     D.parents{k}={'parameters'};
     D.cpds{k}=rvdef('f',@(p)pvals(p,i),D.values{k});
  end
  D.names=[{'parameters'} D.names];
  D.types=[{'c'} D.types];
  D.obs=[1 D.obs];
  D.parents={{} D.parents{:}};
  v=(1:size(pvals,1))';
  D.cpds=[{rvdef('d',w,v)} D.cpds];
  D.values=[{v} D.values];
  D.sizes=[length(v) D.sizes];
end

model=d2model(D,options);