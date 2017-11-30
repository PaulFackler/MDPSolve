% getmatchfunc Gets a index function for discrete random variables with fixed parameters
% INPUTS
%   rv      : an rv structure with a matchfunc field
%   parents : a cell array of the values of the parents of rv
% OUTPUT
%   rv      : rv with matchfunc set with an appropriate function handle

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

function rv=getmatchfunc(rv,parents)
  if isnumeric(parents)
    parents={parents};
  end
  if any(cellfun(@(x)size(x,2)~=1,parents))
    error('parents must be column vectors')
  end
  n=cellfun(@(x)size(x,1),parents);
  np=prod(n);
  if np~=size(rv.cpt,2)
    error('CPT size does not match size of parents')
  end
  dp=length(parents);
  c1=ones(dp,1); c0=1-parents{dp}(1);
  usematch=false;
  usematch=true;
  for i=1:dp
    if isempty(parents{i}) || any(diff(parents{i})~=1)
      usematch=true;
      break
    end
    if i<dp
      c1(i)=c1(i+1)*length(parents{i+1});
      c0=c0-c1(i)*parents{i}(1);
    end
  end
  if usematch
    try
      rv.parameters=indexfunc(parents);
    catch
      evenspacing=false(1,dp);
      for i=1:dp, 
        diffp=diff(parents{i}); 
        if max(abs(diffp-diffp(1)))<1e-14, evenspacing(i)=true; end
      end
      rv.parameters=getcptm(parents,evenspacing);
    end
  else
    rv.parameters=getcptu(c0,c1);
  end

% use these to avoid memory bug in Matlab anonymous functions
function matchfunc=getcptm(s,evenspacing)
  matchfunc=@(varargin)gridmatch(varargin,s,evenspacing);

function matchfunc=getcptu(c0,c1)
  if all(size(c1)==1) % a single parent
    if c0==0
      matchfunc=@(parent)parent;
    else
      matchfunc=@(parent)parent+c0;
    end
  else
    matchfunc=@(varargin)cellprodsum(c0,c1,varargin{:});
  end
  
  function x=cellprodsum(c0,c1,varargin)
    x=c0;
    for i=1:length(c1)
      x=x+varargin{i}*c1(i);
    end
    