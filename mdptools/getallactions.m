% getallactions Obtains all actions to lead to the same value function
% USAGE
%   [A,v]=getallactions(model,results);
% INPUTS
%   model   : mdpsolve model structure 
%   results : mdpsolve results structure
% OUTPUT
%   A     : ns-element cell array - ith element contains index of rows of X
%             associated with an optimal action
%   v     : the ns element value function vector (v is refined to ensure that
%             it is accurately computed)
%
% In some cases DP models have degenerate solutions in which alternative actions 
%   lead to the same value function. This procedure finds all equivalent actions
%   for each state that lead to the same value.

function [A,v]=getallactions(model,results)
Ia=results.Ixopt;
v=results.v;
R=model.R;
P=model.P;

if isnumeric(P)
  bellop=@(v) R + model.d*(P'*v);
else
  bellop=@(v) R + model.d*P(v);
end

% refine the value function
tol=1e-13;
ev=bellop(v);
for i=1:500, 
  if max(abs(v-ev(Ia)))<tol, break; end 
  v=ev(Ia);
  ev=bellop(v); 
end
%disp(i)


if isfield(model,'Ix')
  Ix=model.Ix;
else
  Ix=getI(model.X,model.svars);
end

ns=length(v);
A=cell(ns,1);
nx=length(ev);
for i=1:nx
  si=Ix(i);
  if abs(ev(i)-v(si))<1e-11
    A{si}=[A{si} i];
  end
end

for i=1:ns
  if ~ismember(Ia(i),A{i})
    error('A does not contain the optimal action')
  end
end
