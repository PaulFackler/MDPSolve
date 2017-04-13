% cosdef Defines an rv structure for the Cosine distribution
% USAGE
%   rv=cosdef(parameters);
% or
%   rv=cosdef(parameters,n);
% or
%   rv=cosdef(parameters,values,cpt);
% INPUTS
%   parameters : 2x1 vector [mu;s]
%   n          : positive scalar integer - # of discrete values
%   values     : n-vector of values on (a-b,a+b)
%   cpt        : n-vector of probability weights
% OUTPUT
%   rv         : an rv structure
function rv=cosdef(parameters,values,cpt)
rv=struct('type','cos');
rv.parameters=parameters;
if nargin<2
  rv.values=[];
  rv.cpt=[];
  rv.size=0;
else
  a=parameters(1);
  b=parameters(2);
  if isscalar(values)
    n=values;
    rv.size=n;
    if 1  % gaussian quadrature
      [x,w]=gausslegendre(n,-1,1);
      rv.values=x*b+a;
      w=w.*cos((pi/2)*x);
      rv.cpt=w/sum(w);
    else  % equal weights
      x=(2*(1:n)-1)'/(2*n);
      rv.values=asin(2*x-1)*(2*b/pi)+a;
      rv.cpt=ones(n,1)/n;
    end
  else
    rv.size=length(values);
    rv.values=values;
    if nargin<3
      x=(values-a)/b;
      w=cos((pi/2)*x);
      rv.cpt=w/sum(w);
    else
      rv.cpt=cpt;
    end
  end
end