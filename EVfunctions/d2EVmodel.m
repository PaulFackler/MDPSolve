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
 
 model.X=rectgrid(D.values{Ia},D.values{Is});
 if nargin>1, model.d=delta; end
 model.svars=length(Ia)+(1:length(Is));
 model.EV=true;
 parents=getparents(D);
 fparents=parents(If);
 for i=1:length(If), fparents{i}=matchindices(fparents{i},[Ia Is]); end
 cpts=cellfun(@(x)x.cpt,D.cpds(If),'UniformOutput',false);
 model.P=EVcreate(cpts,model.X,fparents);
 Iu=find(ismember(D.types,'u') | ismember(D.types,'r'));
 if D.cpds{Iu}.type=='f'
   ii=matchindices(parents{Iu},[Ia Is]);
   xi=dvalues(D,ii);
   model.R = D.cpds{Iu}.valfunc(xi{:});
 else
   model.R = D.values{Iu};
 end
 
% matchindices Finds the values of ind1 that match those of ind0
function ind=matchindices(ind0,ind1)
  [~,ind]=ismember(ind0,ind1);
  ind=ind(ind>0);