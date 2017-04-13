% mdpreport Generates a report using the results structure produced by MDPSOLVE
% USAGE
%   mdpreport(results)

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function mdpreport(results)
if ~isstruct(results)
  if iscell(results)
    for i=1:length(results)
      dispmess(results{i})
    end
  elseif isnumeric(results)
    for i=1:length(results)
      dispmess(results(i))
    end
  else
    disp('Unrecognized variable passed to mdpreport')
  end
  return
end

disp(' ')

if isfield(results,'name') && ~isempty(results.name)
  disp(['Problem name: ' results.name])
end

if isfield(results,'warnings') && ~isempty(results.warnings)
  disp('The following warnings were generated:')
  for i=1:length(results.warnings)
    dispmess(results.warnings{i})
  end
end

if isfield(results,'errors') && ~isempty(results.errors)
  disp('The following errors were generated:')
  for i=1:length(results.errors)
    dispmess(results.errors{i})
  end
end

% display convergence results
if isfield(results,'algorithm') && ~isempty(results.algorithm)
  switch results.algorithm(1)
    case 'p'
      if length(results.algorithm)>1
        if results.algorithm(2)=='v'
          disp('MDP solved using policy iteration with vanishing discount approach')
        else
          disp('MDP solved using policy iteration with relative value approach')
        end
      else
        disp('MDP solved using policy iteration')
      end
    case 'f'
      if length(results.algorithm)>1
        if results.algorithm(2)=='v'
          disp('MDP solved using function iteration with vanishing discount approach')
        else
          disp('MDP solved using function iteration with relative value approach')
        end
      else
        disp('MDP solved using function iteration')
      end
    case 'b'
      disp('MDP solved using finite horizon backwards iteration')
  end
end
if isfield(results,'time') && ~isempty(results.time)
  disp(['MDPSOLVE ran for ' num2str(results.time) ' seconds'])
end
if isfield(results,'iter') && ~isempty(results.iter)
  disp(['MDPSOLVE ran for ' num2str(results.iter) ' iterations'])
end
if isfield(results,'change') && ~isempty(results.change)
  disp(['Maximum change in the value function on the last iteration: ' num2str(results.change)])
end
if isfield(results,'numnochange') && ~isempty(results.numnochange)
  disp(['Number of iterations since the last change in the policy:   ' num2str(results.numnochange)])
end
disp(' ')

function dispmess(i)
if iscell(i)
  messnum=i{1};
else
  messnum=i;
end
switch messnum
  case  1, disp(['When Ix is not defined R should be an ns (' num2str(i{2}) ') by na (' num2str(i{3}) ') matrix'])
  case  2, disp(['R should have nx (' num2str(i{2}) ') elements'])
  case  3, disp(['Ix should have nx (' num2str(i{2}) ') elements'])
  case  4, disp(['Ix should be composed of integers between 1 and ns (' num2str(i{2}) ')'])
  case  5, disp(['Iexpand should have nx (' num2str(i{2}) ') elements'])
  case  6, disp('Iexpand should be composed of positive integers')
  case  7, disp('discount has an improper size')
  case  8, disp('discount factor must be defined')
  case 11, disp('If EV option is used P is defined as a function that should accept 1 or 2 arguments')
  case 12, disp('P function is not correctly specified')
  case 13, disp(['If EV option is used the P function should return an nx (' num2str(i{2}) ') by 1 vector'])
  case 14, disp('Expectation function does not interpolate a constant function - checks results carefully') 
  case 15, disp(['Expecting that P should be ns (' num2str(i{2}) ') by nx (' ...
      num2str(i{3}) ') rather than (' num2str(i{4}) ') by nx (' num2str(i{5}) ')'])
  case 16, disp(['Expecting that P should be nx (' num2str(i{2}) ') by ns (' ...
      num2str(i{3}) ') rather than (' num2str(i{4}) ') by nx (' num2str(i{5}) ')'])
  case 17, disp('P and Iexpand are not compatible')
  case 18, disp(['P has columns that do not sum to 1: maximum deviation = ' num2str(i{2})])
  case 19, disp(['P has rows that do not sum to 1: maximum deviation = ' num2str(i{2})])
  case 20, disp(['Either the number of rows or the number of columns of P must equal ns (' num2str(i{2}) ')'])
  case 21, disp('Cannot determine if P is column or row stochastic - options.colstoch must be set')
  case 22, disp(['Cannot determine if P in stage ' num2str(i{2}) ' is column or row stochastic - options.colstoch must be set'])
  case 23, disp('P does not appear compatible with ns and nx')
  case 24, disp(['P does not appear compatible with ns and nx in stage ' num2str(i{2})])
  case 25, disp('nx field is not compatible with Ix, X and/or R fields')
  case 26, disp('P contains NaNs')
  case 27, disp('Cannot link state/actions to states - set Ix')
  case 31, disp('Incorrect specification for algorithm option - when T=inf it must be either ''p'' or ''f''')
  case 32, disp('model.vterm has improper size')
  case 33, disp('options.v has improper size')
  case 35, disp(['NaNs or infinities encountered in updating value function - iterations stopped after ' num2str(i{2}) ' iterations'])
  case 50, disp(['The following warnings were generated in stage ' num2str(i{2})])
  case 51, disp('Policy iteration not implemented with EV option')
  case 52, disp('EV option not allowed with policy iteration - changed to function iteration')
  case 53, disp(['Failure to converge in ' num2str(i{2}) ' iterations'])
  case 61, disp('R is improperly specified')
  case 62, disp('nx is improperly specified or cannot be determined')
  case 63, disp('ns is improperly specified or cannot be determined')
  case 65, disp('model must be a structure variable')
  case 66, disp('cannot link state/action combinations to states')
  case 67, disp(['cannot link state/action combinations to states in stage ' num2str(i{2})])
  case 71, disp('T and horizon fields cannot both be specified')
  case 72, disp('T must be a positive integer or inf')
  case 73, disp('model cannot contain both R and reward fields')
  case 74, disp('Reward (R) must be a function handle or a numeric vector or matrix')
  case 75, disp('model cannot contain both discount and d fields')
  case 76, disp('Discount rate (d) must be a function handle or a numeric scalar or vector')
  case 77, disp('model cannot contain both P and transprob fields')
  case 78, disp('Transition matrix (P) must be a function handle or a numeric matrix')
  case 79, disp('ns and nx must contain positive integers')
  case 80, disp('nrep must contain positive integers')
  case 81, disp('nstage must be a positive integer')
  case 82, disp([i{2} ' is incompatible with nstage'])
  case 83,  
    if i{2}==1
      disp('If R is a matrix it must have ns rows')
    else
      disp(['If R is a matrix it must have ns rows - check stage ' num2str(i{2})])
    end
  case 84, disp('X, svars or Ix has an improper number of stages');
  case 85, 
    if i{2}==1
      disp('Either X and svars or Ix must be specified')
    else
      disp(['Either X and svars or Ix must be specified - check stage ' num2str(i{2})])
    end
  case 86, 
    if i{2}==1
      disp('Ix and R are not compatible')
    else
      disp(['Ix and R are not compatible - check stage ' num2str(i{2})])
    end
  case 87,
    if i{2}==1
      disp('X and R are not compatible')
    else
      disp(['X and R are not compatible - check stage ' num2str(i{2})])
    end
  case 88, 
    if i{2}==1
      disp('svars must contain positive integers')
    else
      disp(['svars must contain positive integers - check stage ' num2str(i{2})])
    end
  case 89, disp('There seems to be no way to determine the number of states - set model.ns')
  case 90, disp('ns must contain positive integers - set model.ns')
  case 91, 
    if i{2}==1
      disp(['Iexpand should have ' num2str(i{3}) ' values'])
    else
      disp(['Element ' num2str(i{2}) ' of Iexpand should have ' num2str(i{3}) ' values'])
    end
  case 92, 
    if i{2}==1
      disp('EV is true but P is not a function handle')
    else
      disp(['EV is true but P is not a function handle - check stage ' num2str(i{2})])
    end
  case 93, disp([i{2} ' should be numeric or have ' num2str(i{3}) ' stages'])
  case 94, disp('An out of memory error was generated in computing Aopt - use Ixopt instead')
  case 95, disp('Aopt can only be obtained for nrep=1, nstage<=12 and keepall=0')
  case 96, disp('Aopt includes all of the columns of X - to get only actions define model.svars')
  case 97, disp('Unable to obtain Xopt - use getA')
  case 99, 
    if length(i)>1
      disp(['An unrecognized error occurred: ' i{2}.message])
    else
      disp(['An unrecognized error occurred: ' lasterr]) %#ok<LERR>
    end
      
  otherwise
    disp(['Message ' num2str(messnum) ' is not recognized'])
end