% f2cpt Converts a function to a CPT
% USAGE
%   [cpt,pvar]=f2cpt(D,list,cleanup,passforward);
% or
%   D=f2cpt(D,[],[],cleanup);
% INPUTS
%   D           : an influence diagram structure
%   list        : list of variables for replacement
%   cleanup     : cleanup variable passed to rectbas to handle extrapolation
%                   0) no adjustment 1) adjust at end 2) adjust for each shock
%   passforward : a 0/1 scalar or vector indicating which variables defined
%                   by functions should be passed forward (1) or
%                   discretized (0)
% OUTPUTS
%   D           : diagram structure with functions replaced

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2014, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function D=f2cpt(D,list,cleanup,passforward,print)
  d=length(D.sizes);
  if nargin<2 || isempty(list)
    list=1:d; 
  end
  if nargin<3 || isempty(cleanup)
    cleanup=0; 
  end
  if nargin<4 || isempty(passforward)
    passforward=true; 
  end
  if nargin<5 || isempty(print)
    print=0; 
  end
  if all(size(passforward)==[1 1])
    passforward=passforward(1,ones(1,d));
  end
  parents=getparents(D); 
  D.parents=parents;
  options=struct('passforward',passforward);
  list=sort(list);
  for i=length(list):-1:1
    var=list(i);
    cpd=D.cpds{var};
    if ~isfield(cpd,'cpt') || isempty(cpd.cpt)
      if ~isfield(cpd,'values')
         error('cannot create a cpt without specifying values')
      end
      switch cpd.type
        case 'f'
          if ~isfield(cpd,'expvalues')
            [X,pvars,D]=dvalues(D,D.parents{var},options);
            v=double(cpd.valfunc(X{:}));
          else
            v=cpd.expvalues;
            pvars=cpd.expparents;
          end
          cpd.cpt=rectbas(v,D.values{var},[],cleanup);
          D.cpds{var}=cpd;
          D.parents{var}=pvars;
        case 'd'
          [X,pvars,D]=dvalues(D,D.parents{var},options);
          cpd.cpt=cpd.parameters(X{:});
          D.cpds{var}=cpd;
          D.parents{var}=pvars;
        case 'logit'
          if isfield(options,'matrix')
            temp=options.matrix;
          else
            temp=[];
          end
          options.matrix=true;
          [X,pvars,D]=dvalues(D,D.parents{var},options);
          cpd.cpt=cpd.parameters*X';
          D.cpds{var}=cpd;
          D.parents{var}=pvars;
          if isempty(temp), options=rmfield(options,'matrix');
          else              options.matrix=temp;
          end
        otherwise
          error('rv type not recognized')
      end
    end
  end
  D.parents=cellfun(@(x)D.names(x),D.parents,'UniformOutput',false);
  return