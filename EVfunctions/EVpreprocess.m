% EVpreprocess
% INPUTS
%   xelist  : m-element cell array with lists of integers for each state variable: 
%             positive integers represent X values 
%             negative integers represent e values
%   e       : cell array of rv structures (discrete or w/ discrete approximations)
%   options : structure variable 
%               options.order a permutation of 1:m
% OUTPUTS
%   type     : 0 indicates all states variables have unique conditioning variables
%              1 indicates that all state variables have unique conditioning shocks
%              2 indicates a general model
%   Ax       : logical matrix with A(i,j)=true if X(j) is associated with S(i)
%   Ae       : logical matrix with A(i,j)=true if e(j) is associated with S(i)
%   shocksum : shock j can be summed out after state shocksum(j) is processed
%                shocksum=0 unless type=2
%   reordere : shocks should be reordered so the first to be summed out are lowest in order
%
%   ne       : sizes of the random noise terms
%   ww       : m-element cell array of weight vectors used for summing out in the ith 
%                step of the processing sequence
function [modeltype,p1,parents1,e1,eenter,eexit,ne,ww,reordere,Ae]=EVpreprocess(p,parents,X,e,options)
if isnumeric(parents), parents={parents}; end
if ~iscell(e), e={e}; end
m=length(parents);
order=[];
print=0;
if exist('options','var') && ~isempty(options)
  if isfield(options,'order'),     order = options.order;     end
  if isfield(options,'print'),     print = options.print;     end
end
if ~isempty(order)
  if length(order)~=m
    error(['options.order must be an ' num2str(m) '-vector'])
  end
  p=p(order);
  parents=parents(order);
end

maxx=max(0,max(cellfun(@max,parents)));  % largest X variable index
maxe=max(0,max(-cellfun(@min,parents))); % largest e variable index

Ax = false(m,maxx);
Ae = false(m,maxe);
for i=1:m
  Ax(i, parents{i}(parents{i}>0))=true;   % Ax(i,j)=1 if X variable j conditions target variable i
  Ae(i,-parents{i}(parents{i}<0))=true;   % Ae(i,j)=1 if e variable j conditions target variable i
end
  
p1=p; parents1=parents;
e1=[]; eenter=[]; eexit=[]; ne=[]; ww=repmat({[]},1,m); reordere=[]; 
modeltype=0;          % disjoint conditioning sets            
if any(sum(Ax,1)>1) || any(sum(Ae,1)>1)
  if all(sum(Ae,1)<=1)
    modeltype=1;      % non-disjoint conditioning sets; no noise variables 
  else
    modeltype=2;      % non-disjoint conditioning sets; with noise variables 
  end
end

if print>0
  switch modeltype
    case 0
      disp('Model type = 0: target variables have non-overlapping conditioning sets')
    case 1
      disp('Model type = 1: target variables have overlapping conditioning sets with no conditioning noise variables')
    case 2
      disp('Model type = 2: target variables have overlapping conditioning sets with conditioning noise variables')
  end
end

% e variables with eenter=0 are in no conditioning sets and are ignored
% e variables with eenter=eexit are in a single conditioning set and 
% should be summed out of the associated p and removed.
if maxe>0
  eenter = zeros(1,maxe);  % iteration where e(i) enters combined parents
  eexit  = zeros(1,maxe);   % iteration where e(i) exits combined parents
  for j=1:maxe
    eenter(j) = find(Ae(:,j),1,'first'); 
    eexit(j)  = find(Ae(:,j),1,'last'); 
  end
  
  % reorder so lowest numbered noise variables are summed out first
  sumnow = eenter==eexit & eenter>0; 
  [~,reordere]=sort(eexit.*double(~sumnow)); 
  if print && ~isequal(reordere,1:maxe)
     fprintf('Noise variables have been renumbered\n',i)
     fprintf('  Old numbering: '); fprintf('%1.0f ',1:maxe); fprintf('\n')
     fprintf('  New numbering: '); fprintf('%1.0f ',reordere); fprintf('\n')
  end
  % move all noise variables to be summed out immediately to the front
  eexit  = eexit(reordere);
  eenter = eenter(reordere);
  sumnow = eenter==eexit & eenter>0; 
  Ae=Ae(:,reordere);
  e1=e(reordere);
  w=cellfun(@(x)x.cpt,e1,'UniformOutput', false);
  ne=cellfun(@length,w);
  parents1=parents;
  % renumber the noise terms in the xelists
  for i=1:m
    for j=1:length(parents1{i}), 
      if parents1{i}(j)<0, parents1{i}(j)=-reordere(-parents1{i}(j)); end
    end
  end
  
  
  % rearrange p to conform to new order and sum out noise variables associated with a single CPT
  for i=1:m
    sumevars=find(Ae(i,:) &  sumnow);
    if isempty(sumevars)
      [p1{i},parents1{i}]=reorderp(p{i},parents1{i},X,e);
    else
      [p1{i},parents1{i}]=reorderp(p{i},parents1{i},X,e,sumevars);
      % remove any variables that have been summed out
      eenter(abs(parents{i}(~ismember(parents{i},parents1{i}))))=0;
      eexit(abs(parents{i}(~ismember(parents{i},parents1{i}))))=0;
      if print>0
          if length(sumevars)>1, ss='s'; else ss=''; end
          fprintf('Target variable %1.0f is associated uniquely with the following noise variable%1c: ',i,ss)
          fprintf('%1.0f  ',sumevars); fprintf('\n')
      end
    end
    if print>0 && ~isequal(parents{i},parents1{i})
       fprintf('Conditioning set for target variable %1.0f has been rearranged\n',i)
       fprintf('  Old parent list: '); fprintf('%1.0f ',parents{i}); fprintf('\n')
       fprintf('  New parent list: '); fprintf('%1.0f ',parents1{i}); fprintf('\n')
    end
  end

  % get the product of any probability vectors of the summed noise terms 
  % at each iteration
  ww=cell(1,m);
  for i=2:m
    ind = find(eexit==i);
    if ~isempty(ind)
      ww{i}=w{ind(1)};
      for j=2:length(ind)
        ww{i}=ww{i}*w{ind(j)}';
        ww{i}=ww{i}(:);
      end
    end
  end
end


return


function [p1,parents1]=reorderp(p,parents,X,e,sumevars)
if nargin<5, sumevars=[]; end
xvars = sort(parents(parents>0));
evars = sort(abs(parents(parents<0)));
parents1=[-sort(evars) sort(xvars)];
ii=matchindices(parents,parents1);
% check if any rearrangement is necessary
if any(diff(ii)<=0)
  [~,Xi]=getI(X,xvars);
  evals=cellfun(@(e)e.values,e,'UniformOutput', false);
  Xe=rectgrid(evals(evars),Xi);
  I=getI(Xe,ii);
  p1=p(:,I);
else
  p1=p;
end
if ~isempty(sumevars)
  ii=matchindices(sumevars,evars);
  if length(sumevars)>1, 
    w=cellfun(@(e)e.cpt,e(ii),'UniformOutput', false);
    w=ckron(w); 
  else
    w=e{ii}.cpt;
  end
  p1=reshape( reshape(p,[],length(w))*w, size(p,1), [] );
  ii=matchindices(-sumevars,parents1);
  parents1(ii)=[];
end


% matchindices Finds the values of ind1 that match those of ind0
function ind=matchindices(ind0,ind1)
  [~,ind]=ismember(ind0,ind1);
  ind=ind(ind>0);



