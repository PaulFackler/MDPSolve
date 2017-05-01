% geneticcomb Genetic optimization solver for combinatoric problems
% Solves optimization problems over permutations on length n
% USAGE
%   [xopt,fbest,pop,fhistory] = geneticcomb(f,n,options);
% INPUTS
%   f       : a function handle with syntax f(x) where x is an n x m that 
%             returns an m-vector of costs associated with the m columns of x.
%             Each column of x is a permutation of the integers 1 through n
%   n       : size of the inputs
%   options : optional structure variable (described below)
% OUTPUTS
%   xopt        : best permutation found
%   fbest       : best (lowest) cost associated with xopt
%   pop         : final generation population (can be placed in
%                   options.results to copntinue search)
%   fhistory    : best cost at each iteration
%
% Options fields:
%  'results',      []      results structure or initial set of points
%  'popSize',     100      size of each generation (must be divisible by 4)
%  'fopt',       -inf      termination level
%  'tol',        1e-8      termination tolerance
%  'numiter',   10000      number of iterations to perform
%  'maxnochange', 500      maximum time w/o improvement
%  'nowrap',        0      set to 1 to prevent wrapping of recombination operator
%  'showiter',      0      = 1 (0) if you do (not) want to see each iter
%  'maxflag',       0      set to 1 for max problems, 0 for min problems 
%  'vecflag',       0      = 1 (0) if f does (not) accept point matrices

% Based on tsp_ga (release 2.3,11/07/11) written by Joseph Kirk (jdkirk630@gmail.com)
% Available at: 
%   http://www.mathworks.com/matlabcentral/fileexchange/13680-traveling-salesman-problem-genetic-algorithm
% Copyright (c) 2007, Joseph Kirk
% All rights reserved.

% Copyright (c) 2013, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function [xopt,fbest,pop,fhistory] = geneticcomb(f,n,options)
if nargin<3, options=[]; end
getopts(options, ...
 'results',     [],...       % results structure or initial set of points
 'popSize',     100,...      % size of each generation
 'fopt',        -inf,...     % termination level
 'tol',         1e-8,...     % termination tolerance
 'numiter',     10000,...    % number of iterations to perform
 'maxnochange', 500,...      % maximum time w/o improvement
 'nowrap',      false,...    % set to 1 to prevent wrapping of recombination operator
 'showiter',    false,...    % = 1 (0) if you do (not) want to see each iter 
 'maxflag',     false,...    % set to 1 for max problems, 0 for min problems
 'vecflag',     false);      % = 1 (0) if f does (not) accept point matrices

% adjust option values
popSize  = 4*ceil(popSize/4);
numiter  = max(1,round(real(numiter(1))));
showiter = logical(showiter(1));
vecflag  = logical(vecflag(1));
maxflag  = logical(maxflag(1));

% Initialize the Population
if isempty(results)
  pop = zeros(n,popSize);
  pop(:,1) = (1:n)';
  for k = 2:popSize
    pop(:,k) = randperm(n)';
  end
else
  if isnumeric(results)
    pop=results;
  else
    pop=results.pop;  % not currently implemented
  end
  [n,popSize]=size(pop);
end
  
% Run the GA
fbest = inf;
totalCost = zeros(1,popSize);
fhistory = NaN+zeros(1,numiter);
if showiter
  disp('iteration   iterations with no change   lowest cost')
end
nochange=0;
for iter = 1:numiter
    % cost of each population member 
    if vecflag
      totalCost(:) = f(pop);
    else
      for i=1:popSize, totalCost(i) = f(pop(:,i)); end
    end
    if maxflag, totalCost=-totalCost; end
    % best permutation in the population
    [minDist,index] = min(totalCost);
    fhistory(iter) = minDist;
    if minDist < fbest
      nochange=0;
      fbest = minDist;
      xopt = pop(:,index);
    else
      nochange=nochange+1;
    end
    if showiter, 
      fprintf('%d   %d   %1.4f\n',iter,nochange,fbest); 
    end
    % produce next generation
    randomOrder = randperm(popSize);
    %routeInsertionPoints = sort(ceil(n*rand(2,popSize/4)),1);
    routeInsertionPoints = ceil(n*rand(2,popSize/4));
    if nowrap, routeInsertionPoints=sort(routeInsertionPoints,1); end
    for p = 4:4:popSize
        ii=randomOrder(p-3:p);
        [ignore,idx] = min(totalCost(ii)); 
        bestOf4 = pop(:,ii(idx));
        I = routeInsertionPoints(1,p/4);
        J = routeInsertionPoints(2,p/4);
        if I<=J, ind=I:J; 
        else     ind=[I:n 1:J]; 
        end
        pop(:,ii)        = bestOf4*ones(1,4);
        pop(ind,ii(2))   = bestOf4(ind(end:-1:1));        % Flip
        pop([I J],ii(3)) = bestOf4([J I]);                % Swap
        pop(ind,ii(4))   = bestOf4([ind(2:end) ind(1)]);  % Slide
    end
    % check for convergence
    if nochange>maxnochange, break; end   % end if no improvement after specified number of iterations
    if fbest<fopt+tol,   break; end       % end if close enough to known min
end
fhistory=fhistory(1:iter);
if maxflag, fbest=-fbest; fhistory=-fhistory; end


% GETOPTS Returns options values in an options structure
% USAGE
%   [value1,value2,...]=getopts(options,field1,default1,field2,default2,...)
% INPUTS
%   options  : a structure variable
%   field    : a field name
%   default  : a default value 
% OUTPUTS
%   value    : value in the options field (if it exists) or the default value
%
% Variables with the field names will be created in the caller's workspace
% and set to the value in the option variables field (if it exists) or to the 
% default value.
% 
% Example called from a function:
%   getopts(options,'tol',1e-8,'maxits',100);
% where options contains the single field 'tol' with value equal to 1
% The function have two variable defined in the local workspace, tol with a
% value of 1 and maxits with a value of 100.
%
% If options contains a field name not in the list passed to getopts, a
% warning is issued.

function varargout=getopts(options,varargin)
K=fix(nargin/2);
if nargin/2==K
  error('fields and default values must come in pairs')
end
if isa(options,'struct'), optstruct=1; else optstruct=0; end
varargout=cell(K,1);
k=0;
ii=1;
for i=1:K
  if optstruct && isfield(options,varargin{ii})
    assignin('caller',varargin{ii},getfield(options,varargin{ii}));
    k=k+1;
  else
    assignin('caller',varargin{ii},varargin{ii+1});
  end
  ii=ii+2;
end
  
if optstruct && k~=size(fieldnames(options),1)
  warning('options variable contains improper fields')
end
