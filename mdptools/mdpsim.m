% mdpsim Simulation of a controlled Markov process
% USAGE
%   [SI,XI,err] = mdpsim(P,s0,T,Ix,colstoch)
% INPUTS
%   P          : ns x nx state transition matrix (or nx x ns if colstoch=0)
%   s0         : k by 1 vector of initial states
%   T          : number of simulated time periods
%   Ix         : ns-vector defining the control - this is an index vector; the ith element
%                  determines which column of P is associated with the ith state
%   colstoch   : 0/1 variable - 1 if P is in colstoch form
% OUTPUTS
%   SI   : k x T+1 matrix of simulated state indices
%   XI   : k x T+1 matrix of simulated state/action indices
%   err  : error code (if non-zero there is a problem)
%            1: P and colstoch appear to be incompatible
%            2: P appears to be invalid
%            3: Can't determine whether P is row or column stochastic - set colstoch

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2013, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function [SI,XI,err] = mdpsim(P,s0,T,Ix,colstoch)
err=0;
if nargin<4 || isempty(Ix)
  if size(P,1)~=size(P,2)
    error('P must be square if Ix is not defined')
  end
  noIx=true;
  Ix=(1:size(P,1))'; 
else
  if length(Ix)~=size(P,1)
    error('P iand Ix are not compatible')
  end
  noIx=false;
end
% check if P is valid
%Pcode=checkP(P,numel(Ix),max(Ix));
Pcode=checkP(P,numel(Ix),size(P,2));
if nargin>=5
  if (colstoch && ~Pcode==1) || (~colstoch && ~Pcode==0)
    err=1;
    if nargout<3
      error('P and colstoch appear to be incompatible')
    end
  elseif Pcode==-1
    err=2;
    if nargout<3
      error('P appears to be invalid')
    end
  end
else
  switch Pcode
    case -2
      err=3;
      if nargout<2
        error('Can''t determine whether P is row or column stochastic - set colstoch')
      end
    case -1
      err=2;
      if nargout<3
        error('P appears to be invalid')        
      end
    case 0
      colstoch=false;
    case 1
      colstoch=true;
  end
end

if err>0
  SI=[]; XI=[]; return
end

k = length(s0);
SI=zeros(k,T+1);
if colstoch
  ns = size(P,1);
  u  = ones(ns,1);
  if isempty(Ix), cp=cumsum(P,1);  
  else            cp=cumsum(P(:,Ix),1);  
  end
  s=s0(:);
  SI(:,1) = s;
  for t=2:T+1
    r = rand(1,k); 
    s = 1+sum(r(u,:)>cp(:,s),1);
    SI(:,t)=s';
  end
else
  ns = size(P,2);
  u  = ones(ns,1);
  if isempty(Ix), cp=cumsum(P,2);  
  else            cp=cumsum(P(Ix,:),2);  
  end
  s=s0(:);
  SI(:,1) = s;
  for t=2:T+1
    r = rand(k,1); 
    s = 1+sum(r(:,u)>cp(s,:),2); 
    SI(:,t)=s';
  end
end

if nargout>1
  if noIx
    XI=[];
  else
    XI=zeros(k,T+1);
    XI(:)=Ix(SI(:));
  end
end