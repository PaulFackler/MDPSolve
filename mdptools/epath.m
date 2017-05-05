% epath Time paths for expectations of functions of Markov chains
% USAGE
%   Ef=epath(fs,P,Ix,T,colstoch);
% INPUTS
%   fs : ns-vector of values of Ix function evaluated at each of the states: f(S)
%   P  : ns x nx transition probability matrix (or nx x ns if colstoch=0)
%   Ix : ns-vector defining the decision strategy
%   T  : time horizon (positive integer)
%   colstoch : 0/1 variable; 1 if P is in column stochastic form
% OUTPUT
%   Ef : ns x (T+1) matrix of expected time paths
%
% P(:,Ix) is an ns x ns matrix defining the transition matrix associated with a
% specific decision strategy.
% If P is in row stochastic form P(Ix,:) is used.

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

function Ef=epath(fs,P,Ix,T,colstoch)
% if colstoch is not passed attempt to determine if P is column stochastic
if nargin<5 || isempty(colstoch)
  sizeP=size(P);
  if sizeP(1)>sizeP(2),     colstoch=false; ns=sizeP(2);
  elseif sizeP(1)<sizeP(2), colstoch=true;  ns=sizeP(1);
  else
    error('Cannot tell if matrix is in column or row stochastic form - pass colstoch')
  end
elseif colstoch, ns=size(P,1);
else             ns=size(P,2);
end
% check the size of the function value vector
if numel(fs)~=ns
  error('fs vector is not compatible with P')
else
  fs=fs(:);
end
    
Ef=zeros(ns,T+1);
% P in column stochastic form
if colstoch
  fs=fs(:)';
  % get the ns x ns transition matrix
  if ~isempty(Ix)
    Pa=P(:,Ix);
  else
    Pa=P;
  end
  % loop over time values
  Ef(:,1)=fs';
  for t=2:T+1
    fs=fs*Pa;
    Ef(:,t)=fs';
  end
% P in row stochastic form
else
  fs=fs(:);
  % get the ns x ns transition matri
  if ~isempty(Ix)
    Pa=P(Ix,:);
  else
    Pa=P;
  end
  % loop over time values
  Ef(:,1)=fs;
  for t=2:T+1
    fs=Pa*fs;
    Ef(:,t)=fs;
  end
end