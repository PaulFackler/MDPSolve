% evalstrategy Value function associated with a specified strategy
% USAGE 
%   v=evalstrategy(model,strategy);
% INPUTS
%   model    : an mdpsolve model structure (see mdpsolve)
%   strategy : an ns-element index vector specifying the rows of X
%                 defining a strategy
% OUTPUT
%   v : an ns-vector specifying the value function associated with
%         the specified strategy
% Note: if a results structure is obtained from mdpsolve
%   then results.v is identical to evalmodel(model,results.Ixopt)
function [v]=evalstrategy(model,strategy)

d=model.d;
if isfield(model,'T') && model.T<inf % finite horizon
  T=model.T;
  if size(strategy,2)==1
    warning('This is a finite time horizon problem but the strategy is not a function of time')
    P=model.P(:,strategy);
    R=model.R(strategy);
    if isfield(model,'vterm')
      v(:,T)=R+delta*(P'*vterm);
    else
      v(:,T)=R;
    end
    for t=T-1:-1:1
      v(:,t)=R+d*(P'*v(:,t+1));
    end
  else
    P=model.P;
    R=model.R;
    v=zeros(size(P,1),T);
    st=strategy(:,T);
    if isfield(model,'vterm')
      v(:,T)=R(st)+delta*(P(:,st)'*vterm);
    else
      v(:,T)=R(st);
    end
    for t=T-1:-1:1
      st=strategy(:,t);
      v(:,t)=R(st)+d*(P(:,st)'*v(:,t+1));
    end
  end
else % infinite horizon
  P=model.P(:,strategy);
  R=model.R(strategy);
  v=(R'/(speye(size(P,1))-d*P))';
end