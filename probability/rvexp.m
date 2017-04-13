% rvexp Generates expected values of functions of a specified distribution
% USAGE 
%   ef=rvexp(rv,f,n,parents);
% INPUTS
%   rv       : an rv structure (define using rvdef)
%   f        : a function handle that accepts and returns an n-vector
%   n        : number of discretization values to use
%   parents  : a cell array composed of n-vectors
%                if parameters is a function, parents are passed to it get
%                conditional parameter values
% OUTPUTS
%    ef : an approximation of E[f(x)]
%
% See rvdef for available distributions
%
function ef=rvexp(rv,f,parents)
if ~isempty(rv.values) && ~isempty(rv.weights)
  ef=f(rv.values)'*rv.weights;
else
  if isempty(n)
    error('n (# of nodes) must be specified')
  end
  if nargin<4 || isempty(parents)
    p=rv.parameters();
  else
    p=rv.parameters(parents);
  end
  rv2=rvdef(rv.type,p,n);
  ef=rvexp(rv2,f);
end
