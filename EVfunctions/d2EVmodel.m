% d2EVmodel Converts a diagram to a model with an EV function
% USAGE
%   model=d2EVmodel(D,delta);
% INPUTS
%   D     : a diagram structure
%   delta : a discount factor, i.e., a scalar on (0,1]
% OUTPUT
%   model : a model structure that can be passed to mdpsolve
%
% Currently this only works for diagrams in which the future states
%   have CPTs that are functions of the states and actions only and
%   the reward is a set of values or a function of the states and actions only.
function model=d2EVmodel(D,delta)
 Is=find(ismember(D.types,'s'));
 Ia=find(ismember(D.types,'a'));
 If=find(ismember(D.types,'f'));
 
 model.X=rectgrid(D.values{Is},D.values{Ia});
 if nargin>1, model.d=delta; end
 model.svars=1:length(Is);
 model.EV=true;
 parents=getparents(D);
 parents=parents(If);
 for i=1:length(If), parents{i}=matchindices(parents{i},[Is Ia]); end
 cpts=cellfun(@(x)x.cpt,D.cpds(If),'UniformOutput',false);
 model.P=EVcreate(cpts,model.X,parents);
 Iu=find(ismember(D.types,'u') | ismember(D.types,'r'));
 if D.cpds{Iu}.type=='f'
   xi=dvalues(D,D.parents{Iu});
   model.R = D.cpds{Iu}.valfunc(xi{:});
 else
   model.R = D.values{Iu};
 end
 
% matchindices Finds the values of ind1 that match those of ind0
function ind=matchindices(ind0,ind1)
  [~,ind]=ismember(ind0,ind1);
  ind=ind(ind>0);