% EVsetup Creates index vectors used to create an EV function
% USAGE
%   [Is,workspace]=EVsetup(X,parents,options);
% INPUTS
%   X       : matrix of state/action combinations
%   parents : m-element cell array of conditioning variables (parents)
%   options : structure variable, fields listed below
% OUTPUTS
%   Is        : m-element cell array of index vectors to define CPTs
%                 values in columns of ith CPT are X(Is{i},parents{i})
%   workspace : a structure variable containing the index vectors
%                 used by an EV function
% parents{i} contains the columns of X that condition the ith variable
% Thus Xi=X(Is{i},parents{i}) contains the unique elements of the
%   conditioning variables that condition output variable i.
% The conditional Probability Table (CPT) associated with variable i should
%   have columns that match the variables in Xi.
%
% options fields:
%   order : order to process the variables (empty for 1:m)
%
% To create an EV function use the following
%   [Is,ws]=EVsetup(X,parents);
%   for i=1:m
%     Xi=X(Is{i},parents{i});
%     p{i} = f(Xi);  % define the CPT here
%   end
%   EV = EVcreate(p,ws);
function [Is,ws]=EVsetup(X,parents,options)
m=length(parents);
order=[];
if exist('options','var') && ~isempty(options)
  if isfield(options,'order'),     order = options.order;     end
end
if ~isempty(order)
  if length(order)~=m
    error(['options.order must be an ' num2str(m) '-vector']);
  end
  parents=parents(order);
end
Is=cell(1,m); Ip=cell(1,m); Iy=cell(1,m); Jp=cell(1,m); Jy=cell(1,m);
parentscombined=cell(1,m);
parents{1}=parents{1}(:);
parentscombined{1}=sort(parents{1});
for i=2:m
  parents{i}=parents{i}(:);
  parentscombined{i}=union(parentscombined{i-1},parents{i});
end
Iy{m}=[]; Jy{m}=[];
Xcombined=X; Ic=[]; Jc=[];
for i=m:-1:1
  %%%%%%% get expansion indices Ip{i} and Jp{i} for the p{i} %%%%'
  % Ip expands relative to the current combination of parents
  % Jp expands relative to the full X
  % also get the extraction index Is{i} used to form Xi{i}
  % but first check is this has already been obtained
  if isequal(parents{i},parentscombined{i}) 
     Is{i}=Ic;
     Ip{i}=[];
     Jp{i}=Jc;
  else
    % see if parent list has already been used
    skip=false;
    for j=i+1:m
      if isequal(parents{i},parents{j})
        Is{i}=Is{j};
        Jp{i}=Jp{j};
        if isequal(parentscombined{i},parentscombined{j}) 
          Ip{i}=Ip{j};
        else
          ii=associate(parents{i},parentscombined{i});
          [~,Is{i},Ip{i}]=unique(Xcombined(:,ii),'rows');
          if ~isempty(Ic), Is{i}=Ic(Is{i}); end
        end
        skip=true;
        break
      end
    end
    if ~skip
      ii=associate(parents{i},parentscombined{i});
      [~,Is{i},Ip{i}]=unique(Xcombined(:,ii),'rows');
      if ~isempty(Ic), Is{i}=Ic(Is{i}); end
      if isempty(Jc)
        Jp{i}=Ip{i};
      else 
        Jp{i}=Ip{i}(Jc);
      end
    end
  end
   
  %%%%%%% get expansion indices Iy{i} and Jy{i} for y %%%%
  % Iy expands relative to the current combination of parents
  if i>1
    if isequal(parentscombined{i-1},parentscombined{i}) 
      Iy{i}=[];
      Jy{i}=Jc;
    else
      ii=associate(parentscombined{i-1},parentscombined{i});
      [Xcombined,Ici,Iy{i}]=unique(Xcombined(:,ii),'rows');
      if isempty(Jc), Jy{i}=Iy{i};
      else            Jy{i}=Iy{i}(Jc);
      end
      Jc=Jy{i};
      if isempty(Ic),  Ic=Ici;
      else             Ic=Ic(Ici);
      end
    end
  end
  if ~isempty(Is{i}), Is{i}=convert2int(Is{i}); end
  if ~isempty(Ip{i}), Ip{i}=convert2int(Ip{i}); end
  if ~isempty(Iy{i}), Iy{i}=convert2int(Iy{i}); end
  if ~isempty(Jp{i}), Jp{i}=convert2int(Jp{i}); end
  if ~isempty(Jy{i}), Jy{i}=convert2int(Jy{i}); end
end
if ~isempty(order), Is(order)=Is; end
Iy{1}=[];
Jy{1}=[];
ws.Ip=Ip;
ws.Iy=Iy;
ws.Jp=Jp;
ws.Jy=Jy;
ws.order=order;

% Finds the values of ind2 associated with each value of ind1.
% In this application ii(j) is always a value ind2.
% If an element is in ind1 and not ind2 a 0 would be returned.
function ii=associate(ind1,ind2)
m=length(ind1);
ii=zeros(1,m);
for i=1:m
  ii(i)=find(ind1(i)==ind2);
end

function x=convert2int(x,n)
if nargin<2, n=max(x); end
if     n>=2^32, x=cast(x,'uint64');
elseif n>=2^16, x=cast(x,'uint32');
elseif n>=2^8,  x=cast(x,'uint16');
else            x=cast(x,'uint8');
end


