% amdppassive Solves a passive Adaptive Management Markov Decision Problem
% USAGE
%   [v,a,B]=amdppassive(p,P,R,model,options);
% INPUTS
%   p       : belief mesh size (# of belief intervals for each model)
%   P       : an nmodel-element cell array of model specific transition matrices
%   R       : an nmodel-element cell array of model specific reward matrices
%               if empty model.R will be used for all models
%   model   : an MDPSOLVE model structure containing other model information
%   options : an MDPSOLVE options structure
% OUTPUTS
%   v       : ns x Q matrix of value functions
%   a       : ns x Q matrix of optimal actions
%   B       : Q x nmodel matrix of belief weights
%
%  Q = (p+q-1)!/p!/(q-1)! is the number of belief values used to discretize 
%    the belief space (where q=nmodel).
%  ns is the number of states in the no-uncertainty case.

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

function [v,a,B]=amdppassive(p,P,R,model,options)
  if nargin<5, options=struct([]); end
  nmodel=numel(P);
  if all(size(p)==1)
    B=simplexgrid(nmodel,p,1,1);
  else
    B=p;
  end
  nb=size(B,1);
  ns=size(P{1},1);
  v=zeros(ns,nb);
  a=zeros(ns,nb);
  % Only a single R is specified
  if isempty(R)
    % loop over the belief state values
    for i=1:nb
      if nb>1, 
        nchar = fprintf('Processing belief %1i out of %1i\n',i,nb);
      end
      % get the expected transition matrix and rewards
      Pi=B(i,1)*P{1};
      for j=2:nmodel
        Pi=Pi+B(i,j)*P{j};
      end
      model.P=Pi;
      results=mdpsolve(model,options);
      if isempty(results.errors)
        v(:,i)=results.v; a(:,i)=results.Ixopt;
      else
        mdpreport(results);
        error('Unable to solve one of the subproblems')
      end
    end
  % Both P and R have multiple specifications
  else
    % loop over the belief state values
    for i=1:nb
      if nb>1, 
        nchar = fprintf('Processing belief: %5i out of %5i\n',i,nb); 
      end
      % get the expected transition matrix and rewards
      Pi=B(i,1)*P{1};
      Ri=B(i,1)*R{1};
      for j=2:nmodel
        Pi=Pi+B(i,j)*P{j};
        Ri=Ri+B(i,j)*R{j};
      end
      model.P=Pi;
      model.R=Ri;
      results=mdpsolve(model,options);
      if isempty(results.errors)
        v(:,i)=results.v; a(:,i)=results.Ixopt;
      else
        mdpreport(results);
        error('Unable to solve one of the subproblems')
      end
    end
  end
  fprintf('\n')