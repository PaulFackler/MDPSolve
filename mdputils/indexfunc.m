% indexfunc Creates a function index a set of points
% USAGE
%   [Ic,Im]=indexfunc(x);
% INPUT
%   x : n x d matrix or 
%       d-element call array composed of column vectors (x <- rectgrid(x))
% OUTPUTS
%   Ic : function handle of the form I(X1,X2,...,Xd)
%   Im : function handle of the form I(X)
%
% When passed a d column vectors, each of length m, returns a vector of m integer
%   values on {1,...,n} giving the row of x associated with each of the values of X

%g={linspace(0,1,7)',linspace(1,10,10)'}; x=rectgrid(g); II=indexfunc(x); i1=randi(size(X,1),20,1); 
%X=x(i1,:); XX=mat2cell(X,size(X,1),ones(1,size(X,2))); i2=II(XX{:}); isequal(i1,i2)
function [Ic,Im]=indexfunc(x)
if iscell(x)
  x=rectgrid(x);
end
n=size(x,1);
y=(1:n)';
beta=[ones(n,1) x]\y;
alpha=beta(1);
beta=beta(2:end);
e=y-round(alpha+x*beta);
if any(e~=0)
  error('I cannot be found');
end

clear n y e x
Ic=@(varargin) Ifunc(alpha,beta,varargin{:});
Im=@(X) round(alpha+X*beta);


function ival=Ifunc(alpha,beta,varargin)
ival=alpha;
for i=1:length(beta)
  ival = ival + varargin{i}*beta(i);
end
ival=round(ival);