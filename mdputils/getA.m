% getA Get values of variables 
% USAGE
%   A=getA(X,Ia,vars,stages,reps);
% INPUTS
%   X      : nx x k matrix or nstage cell array
%   Ia     : ns vector of variable values on {1,...,k} (e.g. results.Ixopt field)
%   stages : vector of values on {1,...,nstage}
%   reps   : vector of values
% OUTPUT
%   A      : ns row matrix of values or length(stages) cell array
%
% For single stage models the last two imputs can be omitted.
% For any of inputs 2-5 that are omitted or empty, all possible values are returned

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

  function A=getA(X,Ia,vars,stages,reps)
    if nargin<2, Ia     = []; end
    if nargin<3, vars   = []; end
    if nargin<4, stages = []; end
    if nargin<5, reps   = []; end
    if isnumeric(X),                  X  = {X};  end
    if isnumeric(Ia) && ~isempty(Ia), Ia = {Ia}; end
    if isempty(stages)
      nstage=max(numel(X),numel(Ia));
      stages=1:nstage;
    else
      nstage=length(stages);
    end
    A=cell(1,nstage);
    for i=1:nstage
      Xi=X{min(stages(i),numel(X))};
      if isempty(Ia)
        Ii=(1:size(Xi,1))';
      else
        Ii=Ia{min(stages(i),numel(Ia))};
      end
      if isempty(vars)
        if isempty(reps)
          Ai=Xi(Ii,:);
        else
          Ai=Xi(Ii(:,reps),:);
        end
      else
        if isempty(reps)
          Ai=Xi(Ii,vars);
        else
          Ai=Xi(Ii(:,reps),vars);
        end
      end
      if length(reps)>1
        Ai=reshape(Ai,[size(Ii,1),size(Ai,1)/size(Ii,1),size(Ai,2)]);
        Ai=permute(Ai,[1 3 2]);
      end
      A{i}=Ai;
      clear Ai;
    end
    if nstage==1
      A=A{1};
    end