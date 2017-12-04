% HERMITEINV Generates a Hermite interpolating function for an inverse
% USAGE
%    [v,stats]=hermiteinv(f,df,x0,opts);
% INPUTS
%  f        : function name or handle for function the inverse of which is approximated
%  df       : function name or handle for derivative df/dx
%  x0       : an initial sorted vector of input values (contains the desired endpoints)
%  opts     : an options structure (described below)
% OUPUT
%  v        : 3-cell array of results containing {y f^-1(y) df^-1(y/dy)}
%               sorted by y
%  stats    : number of failed steps (steps less than min size allowed)
%
% The function can be interpolated using
%   x=hermiteinterp(v,y);
% 
% opts can contain the following fields  
%   tol   : the absolute error tolerance for the fitted values
%   mindy : the minimimum distance between two adjacent values of y
%   chunk : the number of elements used to expand solution vectors, when needed
% 
% This algorithm bisects each interval and compares the fitted interpolate
% with the function and its derivative. It will continue to bisect and
% interval until the absolute differences in both the function and its
% derivative are less than a given tolerance.
%
% use hermiteinterp to evaluate the interpolating function.
%
% Based on:
% Wolfgang Hörmann and Josef Leydold (2003)
% "Continuous Random Variate Generation by Fast Numerical Inversion"
% ACM Transactions on Modeling and Computer Simulation, 13: 347-362.
function [v,stats]=hermiteinv(f,df,x0,opts)

% set optional parmameters
tol   = 1.5e-8;
mindy = [];
maxdy = [];
chunk = 1000;
if exist('opts','var') && ~isempty(opts)
  if isfield(opts,'tol'),   tol   = opts.tol;   end
  if isfield(opts,'mindy'), mindy = opts.mindy; end
  if isfield(opts,'maxdy'), maxdy = opts.maxdy; end
  if isfield(opts,'chunk'), chunk = opts.chunk; end
end

n=length(x0);

x=x0;
% get initial values
y=feval(f,x); 
d=feval(df,x);
d=1./d;
 
% set minimum and maximum step sizes, as needed
if isempty(mindy),
  mindy=(y(end)-y(1))*1.5e-8;
end
mindy=max(mindy,0);
if isempty(maxdy),
  maxdy=(y(end)-y(1))/25;
end
stats=0;

% initialize loop data 
ind=[(2:n)';0]; % initial links for linked lists
i=1;
j=2;
yi = y(1);
xi = x(1);
di = d(1);
newpoint=0;
% main loop
while j~=0
  % distance between adjacent points
  dy=y(j)-yi;
  xmid=(xi+x(j))/2;
  % distance between evaluation points is too small or
  % mid point evaluates to one of the interval points
  if dy<=mindy || xmid==x(j) || xmid==xi
    newpoint=1;
    stats=stats+1;
  else
    n=n+1;
    % expand vectors if necessary
    if n>length(ind)
      ind = [ind;zeros(chunk,1)];
      y   = [y  ;zeros(chunk,1)];
      x   = [x  ;zeros(chunk,1)];
      d   = [d  ;zeros(chunk,1)];
    end
    % update links
    ind(i)=n;
    ind(n)=j;
    % new point midway between two existing adjacent points
    x(n)=xmid;
    y(n)=feval(f,x(n));
    dn=feval(df,x(n));
    if isinf(dn), d(n)=0;
    elseif dn==0, d(n)=realmax;
    else          d(n)=1/dn;
    end
    if dy<maxdy
      % evaluate approximate function at new point
      dxi    = x(j)-xi;
      % monotonicity check
      mono=3*dxi/dy;
      if di<mono && d(j)<mono
        lambda = (y(n)-yi)./dy;
        a1     = di.*dy;
        a2     = 3*dxi-d(j).*dy-2*a1;
        a3     = dxi-a1-a2;
        xnfit=((a3.*lambda+a2).*lambda+a1).*lambda+xi;
        % evaluate function and derivative at the approximate
        yfit=feval(f,xnfit);
        % check if approximate is close enough
        % if so, make right endpoint the new left endpoint
        if abs(yfit-y(n))<tol
          newpoint=1;
        end
      end
    end
  end
  if newpoint
    i  = j;
    yi = y(i);
    xi = x(i);
    di = d(i);
    newpoint=0;
  end
  % set the new right endpoint
  j=ind(i);
end

% put results into a 3-element cell array {x fx dfx}
v1=zeros(n,1);
v2=zeros(n,1);
v3=zeros(n,1);
j=1;
for i=1:n
  v1(i)=y(j);
  v2(i)=x(j);
  v3(i)=d(j);
  j=ind(j);
end

xi     = v1(1:n-1);
dxi    = v1(2:n) - xi;
yi     = v2(1:n-1);
dyi    = v2(2:n) - yi;
a1     = v3(1:n-1).*dxi;
a2     = 3*dyi - v3(2:n).*dxi - 2*a1;
a3     = dyi - a1 - a2;

v={xi,dxi,yi,a1,a2,a3};
