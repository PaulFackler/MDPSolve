% raisedcosdef Defines an rv structure for the Raised Cosine distribution
% USAGE
%   rv=raisedcosdef(parameters,values,cpt);
% INPUTS
%   parameters : 2x1 vector [mu;s]
%   values     : positive scalar integer - # of discrete values
%                vector of values on (mu-s,mu+s)
%   cpt        : vector of probaiblity weights
% OUTPUT
%   rv         : an rv structure
function rv=raisedcosdef(parameters,values,cpt)
rv=struct('type','raisedcos');
rv.parameters=parameters;
if nargin<2
  rv.values=[];
  rv.cpt=[];
  rv.size=0;
else
  mu=parameters(1);
  s=parameters(2);
  if isscalar(values)
    n=values;
    rv.size=n;
    [x,w]=gausslegendre(n,-1,1);
    %x=((1:n)'./(n+1)-0.5)*2; w=ones(n,1);
    rv.values=x*s+mu;
    w=w.*(1+cos(pi*x));
    rv.cpt=w/sum(w);
  else
    rv.size=length(values);
    rv.values=values;
    if nargin<3
      x=(values-mu)/s;
      w=1+cos(pi*x);
      rv.cpt=w/sum(w);
    else
      rv.cpt=cpt;
    end
  end
end