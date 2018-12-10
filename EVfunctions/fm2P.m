% fm2P Creates a transition matrix from a factored model
% USAGE
%   P=fm2P(p,parents,X,e);
% INPUTS
%   p       : m-element cell array of conditional probability matrices
%   parents : m-element cell array of conditioning variable (parent) indices
%   X       : matrix or cell array of state/action variables
%   e       : cell array of rv structures (discrete or w/ discrete approximations)
% OUTPUTS
%   P       : a transition probability matrix
%
% This procedure calls EVcreate with mergevec=m and then extracts the
%   transition matrix from the EV function's workspace 
function [P,ws]=fm2P(p,parents,X,e,options)
if nargin < 4, e={}; end
if nargin<5, options=struct; end
options.mergevec = length(p);
[EV,ws] = EVcreate(p,parents,X,e,options);
P=ws.p{1};
if ~isempty(ws.vreorder)
  %P=P(ii,:);
end
return

% a simple but inefficient alternative
W=[];
evals = {};
Xe = X;
% get information on random noise (e) variables, if any
if nargin>=4 && ~isempty(e)
  % merge the X and e variables
  if ~iscell(e), e={e}; end
  evals=cellfun(@(x)x.values,e,'UniformOutput',false);
  if iscell(X)
    Xe = [evals,X];
  else
    Xe = rectgrid(evals,X);
  end
  %get the joint noise probability vector
  w=cellfun(@(x)x.cpt,e,'UniformOutput',false);
  W=w{1};
  for i=2:length(w)
    W = w{i}*W';
    W = W(:);
  end
end  

% merge the CPTs using kroncol
m=length(p);
for i=1:m
  parentsi = parents{i};
  parentsi(parentsi>0) = parentsi(parentsi>0)+length(evals);
  parentsi(parentsi<0) = -parentsi(parentsi<0);
  ind=getI(Xe,parentsi);
  pi = p{i}(:,ind);
  ind = []; % free memory
  if i==1; P = pi; 
  else     P = kroncol(P,pi);
  end
  pi = [];  % free memory
end

% sum out the noise variables, if any
if ~isempty(W), 
  P=reshape(reshape(P,[],length(W))*W,size(P,1),[]); 
end
