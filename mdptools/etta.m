% etta Expected time to absorption for a recurrent absorbing Markov chain
% USAGE
%   ET=etta(P);
% INPUT
%  P  : ns x ns transition probabiity matrix with the first state absorbing
% OUTPUT
%  ET : ns x 1 expected time to absorption beginning in each state

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

function ET=etta(P)
ns=size(P,1);
if P(1,1)~=1
  error('The first state must be absorbing');
end
if any(P(1,2:end)~=0) && any(P(2:end,1)~=0)
    error('not a valid probability matrix')
end
if all(P(1,2:end)==0)      % row stochastic form
  ET=(sparse(1:ns,1:ns,[2 ones(1,ns-1)])-P)\[0;ones(ns-1,1)];
elseif all(P(2:end,1)==0)  % column stochastic form
  ET=([0 ones(1,ns-1)]/(sparse(1:ns,1:ns,[2 ones(1,ns-1)])-P))';
end