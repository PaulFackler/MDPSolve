% EVindices Creates index vectors for EV functions (called by EVcreate)
% USAGE
%   [Ip,Iy]=EVindices(parents,X,e,eremove);
% INPUTS
%   parents : m-element cell array of conditioning variables (parents)
%               positive integers associate with X variables
%               negative integers associate with e variables            
%   X       : matrix of state/action combinations or 
%               cell array of state/action variables
%   e       : de-element cell array of rv structures (discrete or 
%                w/ discrete approximations)
%   eremove : de-vector specifying which of the noise variables are summed out
%   options : options structure
% OUTPUTS
%   Ip,Iy : index vectors to guide EV function operations 
%
% parents{i} contains the X and e varaibles that condition the ith target variable
% Thus Xpi=unique(X(:,parents{i}),'rows') contains the unique elements of the
%   conditioning variables that condition output variable i.

function [Ip,Iy] = EVindices(parents,X,e,eremove)
m=length(parents);
ne = [];
if nargin<3 || isempty(e)
  e={};
  de = 0;
  eremove=false(0,1);
elseif nargin<4 || isempty(eremove)
  de = length(e);
  ne = cellfun(@(x)length(x.values),e);
  eremove=false(de,1); 
else
  de = length(e);
  ne = cellfun(@(x)length(x.values),e);
end
if length(eremove) == 1, 
  eremove = repmat(eremove,de,1); 
end
if isempty(e)
  eenter = 0; eexit = 0;
else
  [eenter,eexit] = getenterexit(parents,e);
  eexit(~eremove)=0;
end
% redefine the parent vectors to conform to Xe
combined = cell(1,m);
cp=[];
for i=1:m
  parentsi = parents{i};
  parentsi(parentsi>0) = parentsi(parentsi>0) + de;
  parentsi(parentsi<0) = -parentsi(parentsi<0);
  parents{i} = parentsi;
  cp = union(cp,parentsi);
  cp = cp(:)';
  combined{i} = cp; 
  remove = find(eexit==i);
  cp(matchindices(remove,cp)) = [];
end
% merge the X and e variables
evals=cellfun(@(x)x.values,e,'UniformOutput',false);
% get the indices
Ip=cell(1,m);
Iy=cell(1,m);
Xe = combineXe(X,evals,combined{m});
for i = m:-1:1
  if 0
  Xee = combineXe(X,evals,combined{i});
  if i<m
    if ~isequal(Xee,Xe)
      disp(' ')
    end
  end
  end
  ii = matchindices(parents{i},combined{i});
  Ip{i}=getI(Xe,ii);
  if i>1
    [ii,ix] = matchindices(combined{i-1},combined{i});
    [Iy{i},Xe]=getI(Xe,ii);
    % add in the shocks that were removed in the previous stage
    if ~isempty(ix)
      if iscell(Xe)
        Xe=[evals(combined{i-1}(ix)) Xe];
      else
        Xe = rectgrid(evals(combined{i-1}(ix)),Xe);
      end
    end
  else
    Iy{1} = ones(numel(Ip{i}),1);
  end
  Ip{i} = uint64(Ip{i});
  Iy{i} = uint64(Iy{i});
  if 0
    nei = prod(ne(i>=eenter & i<=eexit));
    Ip{i} = reshape(uint64(Ip{i}),[],nei);
    Iy{i} = reshape(uint64(Iy{i}),[],nei);
  end
end
return

function Xe=combineXe(X,evals,vars)
de=length(evals);
if iscell(X),
  Xe = [evals(vars(vars<=de)),X(vars(vars>de)-de)];
else
  Xe = rectgrid(evals(vars(vars<=de)),unique(X(:,vars(vars>de)-de),'rows'));
end

