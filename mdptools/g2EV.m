% g2EV Creates a expected value function
% USAGE
%   [EV,p,Is,ws,warnings]=g2EV(g,s,X,e,xelist,options);
% INPUTS
%   g       : d-element cell array of function handles
%   s       : d-element cell array of column vectors of state variable values 
%   X       : k-row matrix of state/action combinations
%   e       : q-element cell array of rv structures [optional]
%   xelist  : d-element cell array of lists of variable indices
%               element i should contain a list specifying which variables
%               in X and e should be passed to g{i}. Specify X variables
%               with positive integers and e variables with negative integers.
%               Variables must be listed in the order defined by g{i}.
%   options : a structure variable to control the procedure
%               cleanup    : 0/1/2 - Determines how extrapolation is handled
%                            0) no adjustments
%                            1) negative values are set to 0 and values are 
%                                 adjusted so columns sum to 1.
%                            2) values of any variable beyond bounds are set to
%                                 the nearby boundary value
%               rectinterp : 0) use Freudenthal interpolation 
%                            1) use linear interpolation [default: 1]
% OUTPUTS
%   EV       : e function handle that can be used with mdpsolve models
%   p        : d-element cell array of individual variable transition matrices (column stochastic)
%   Is       : d-element cell array of selection indices for X from EVsetup
%   ws       : a cell array of selection indices for the EV function from EVsetup
%   warnings : a cell array of warning messages (if not requested any messages
%                will be displayed)
%
% The functions in g should accept a k-row column vectors of states, action and 
%   noise variables and return a k-row column vector of next period state values

% Programming note:
% the function implements interpolation for simplex and scattered grids as well 
% but these features have not be adequately tested or documented

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2017, Paul L. Fackler (paul_fackler@ncsu.edu)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without  
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%        this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice, 
%        this list of conditions and the following disclaimer in the 
%        documentation and/or other materials provided with the distribution.
%    * Neither the name of the North Carolina State University nor of Paul L. 
%        Fackler may be used to endorse or promote products derived from this 
%        software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% For more information, see the Open Source Initiative OSI site:
%   http://www.opensource.org/licenses/bsd-license.php

function [EV,p,Is,ws,warnings]=g2EV(g,s,X,e,xelist,options)
warnings={};
if nargin<6, options=0; end
geometry=0;
if isnumeric(options)
  cleanup=options;
  rectinterp=1;
else
  getopts(options, ...
         'cleanup',     0, ... % ensures that P is a probability matrix
         'rectinterp',  1, ...
         'order',       []); 
end
if rectinterp==0, geometry=3; end

g=inputcheck(g,'g','function_handle');
s=inputcheck(s,'s','numeric');
e=inputcheck(e,'e','struct');
xelist=inputcheck(xelist,'xlist','numeric');
if ~isnumeric(X),
  error('X must be a numeric matrix')
end

ds=length(s);
dx=size(X,2);
de=length(e);

xlist=cell(1,ds);
elist=cell(1,ds);
for i=1:ds
  if length(xelist{i})~=nargin(g{i})
    error(['xelist has a different number of inputs than g for variable ' num2str(i)])
  end
  elist{i}=-xelist{i}(xelist{i}<0);
  if any(elist{i}>de)
    error(['xelist{' num2str(i) '} contains a value greater than the # of noise variables'])
  end
  for j=1:i-1
    if ~isempty(intersect(elist{i},elist{j}))
      error('g2EV is not implemented for models with overlapping noise variables')
    end
  end
  xlist{i}=xelist{i}(xelist{i}>0);
  if any(xlist{i}>dx)
    error(['xelist{' num2str(i) '} contains a value greater than the # of X variables'])
  end
end

% convert rv structures to cell arrays of value and probability vectors
w=cell(1,de);
for i=1:de
  w{i}=e{i}.cpt;
  e{i}=e{i}.values;
end
[Is,ws]=EVsetup(X,xlist,options);
p=cell(1,ds);
for i=1:ds
  if isempty(Is{i})
    Xi=X(:,xlist{i});
  else
    Xi=X(Is{i},xlist{i});
  end
  nx=size(Xi,1);
  wi=prod(rectgrid(w{elist{i}}),2);
  nw=length(wi);
  Xi=mat2cell(Xi,size(Xi,1),ones(1,size(Xi,2)));
  Xei=cell(1,length(xelist{i}));
  try % expand X and e together
    [Xei(xelist{i}>0)]=cellfun(@(X)repmat(X,nw,1),Xi,'UniformOutput',0);
    [Xei(xelist{i}<0)]=cellfun(@(X)vec(repmat(X',nx,1)),e(elist{i}),'UniformOutput',0);
    sval=g{i}(Xei{:});
    clear Xi Xei
    pi=getPk(sval,s{i},geometry,cleanup);
    p{i}=pi*kron(wi,speye(nx)); 
  catch % loop over values of e
    disp('using low memory approach in g2EV')
    ei=cell(1,length(elist{i}));
    [ei{:}]=rectgrid(e{elist{i}});
    Xei(xelist{i}>0)=Xi;
    clear Xi
    pi=sparse(1,1);
    eik=cell(1,length(elist{i}));
    for j=1:size(ei{1},1)
      for k=1:length(elist{i}), eik{k}=ei{k}(j); end
      Xei(xelist{i}<0)=eik;
      sval=g{i}(Xei{:});
      pi=pi+getPk(sval,s{i},geometry,cleanup)*wi(j);
    end
    p{i}=pi;
  end
end

% eliminate inappropriate values arising from extrapolation
if cleanup==1 && geometry~=1
  for i=1:ds
    p{i}=normalizeP(p{i});
  end
end
EV=EVcreate(p,ws);

if nargout<5
  for i=1:numel(warnings)
    disp(warnings{i})
  end
end


function Pk=getPk(gval,s,geometry,cleanup)
  switch geometry
    case 0
      Pk = rectbas(gval,s,[],cleanup);
    case 1
      Pk = scatterbas(gval,s);
    case 2
      Pk = simplexbas(gval,s.params{:});
    case 3
      Pk = freudenthal(gval,s,cleanup);
  end
  
function x=inputcheck(x,xname,class)
if ~iscell(x)
  if isa(x,'class')
    x={x};
  else
    error([xname ' must be a single element or a cell array of elements of class ' class])
  end
end