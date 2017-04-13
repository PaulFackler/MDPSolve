% mdp_getparams Used by mdpsolve and otehr MDP solvers
% gets stage i parameters values and performs checks if needed
% called once for non-staged models
% called at each stage change for staged models

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

function [errors,warnings,R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP] = ...
  mdp_getparams(i,checks,R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP,debug)

if isempty(EV),       EV=false;      end

if iscell(R), R=R{min(numel(R),i)}; end
if iscell(P), P=P{min(numel(P),i)}; end
if iscell(d), d=d{min(numel(d),i)}; end
if iscell(Ix), Ix=Ix{min(numel(Ix),i)}; end
if iscell(Iexpand), Iexpand=Iexpand{min(numel(Iexpand),i)}; end
if iscell(EV), EV=EV{min(numel(EV),i)}; end
ns=ns(min(numel(ns),i));
nx=nx(min(numel(nx),i));
Xindexed=Xindexed(min(numel(Xindexed),i));
expandP=expandP(min(numel(expandP),i));
EV=EV(min(numel(EV),i));

if i>1, nsnext=ns(min(numel(ns),i-1));
else    nsnext=ns(min(numel(ns),length(ns)));
end
if ~isnumeric(P) && ~EV, P=P(); end
if ~isnumeric(R), R=R(); end
if ~isnumeric(d), d=d(); end
if ~isnumeric(Ix), Ix=Ix(); end
if ~isnumeric(Iexpand), Iexpand=Iexpand(); end

% make sure these are columns vectors
if ~isempty(d), d=d(:); end
if ~isempty(Ix), Ix=Ix(:); end
if ~isempty(Iexpand), Iexpand=Iexpand(:); end

if ~isempty(colstoch)
  colstoch=colstoch(min(numel(colstoch),i));
end

if checks
  [errors,warnings]=mdp_checks(R,P,d,ns,nx,Ix,Iexpand,colstoch,Xindexed,expandP,EV,debug);
else
  errors={};
  warnings={};
end