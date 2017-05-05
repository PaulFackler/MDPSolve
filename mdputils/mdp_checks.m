% mdpchecks Checks input data for MDP problems
% USAGE
%   [errors,warnings]=mdp_checks(R,P,d,ns,nx,Ix,Iexpand,colstoch,Xindexed,expandP,EV);
% OUTPUTS
%   errors, warnings : cell arrays with message numbers and supporting information
%
% use mdpreport(errors) or mdpreport(warnings) to display messages

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

function [errors,warnings]=mdp_checks(R,P,d,ns,nx,Ix,Iexpand,colstoch,Xindexed,expandP,EV,debug)
 errors={};
 warnings={};
 if ~Xindexed
    na=nx/ns;
    if any(size(R)~=[ns na])
      if debug, error(' '); end
      errors{end+1}={1,ns,na};
    end
 else
   if numel(R)~=nx    
     if debug, error(' '); end 
     errors{end+1}={2,nx};
   end
   if numel(Ix)~=nx
     if debug, error(' '); end
     errors{end+1}={3,nx};
   end
   if any(Ix~=round(Ix)) || any(Ix<1) || any(Ix>ns)
     if debug, error(' '); end
     errors{end+1}={4,ns};
   end
 end
 
 if expandP
   if numel(Iexpand)~=nx
     if debug, error(' '); end
     errors{end+1}={5,nx};
   end
   if any(Iexpand~=round(Iexpand)) || any(Iexpand<1)
     if debug, error(' '); end
     errors{end+1}=6;
   end
 end
 
 
 if numel(d)~=1 && numel(d)~=ns
     if debug, error(' '); end
   errors{end+1}=7;
 end
      
 
 tol=1e-16*ns; 
 if EV
   if nargin(P)>2
     if debug, error(' '); end
     errors{end+1}=11;
   end
   if nargout(P)>1
     if debug, error(' '); end
     errors{end+1}=11;
   end
   try
     if nargin(P)==1
       ev=P(ones(ns,1));
       if ndims(ev)>2 || any(size(ev)~=[nx 1])
         if debug, error(' '); end
         errors{end+1}={13,nx};
       end
     else
       ev=P(ones(ns,1),(1:ns)');
       if ndims(ev)>2 || any(size(ev)~=[ns 1])
         if debug, error(' '); end
         errors{end+1}={13,ns};
       end
     end
     if any(abs(ev-1)>tol)
       warnings{end+1}=14;
     end
   catch %#ok<CTCH>
     if debug, error(' '); end
     errors{end+1}=12;
   end
 else
    if any(any(isnan(P)))
      if debug, error(' '); end
      errors{end+1}={26};
    end
    sizeP=size(P); 
    if colstoch
      if expandP
        if sizeP(2)<max(Iexpand)
          if debug, error(' '); end
          errors{end+1}=17;
        end
      else
        if sizeP(1)~=ns || sizeP(2)~=nx
          if debug, error(' '); end
          errors{end+1}={15,ns,nx,size(P,1),size(P,2)};
        end
      end
      sP=max(abs(sum(P)-1));
    else
      if expandP
        if sizeP(1)<max(Iexpand)
        if debug, error(' '); end
        errors{end+1}=17;
        end
      else
        if sizeP(1)~=nx || sizeP(2)~=ns
          if debug, error(' '); end
          errors{end+1}={16,nx,ns,size(P,1),size(P,2)};
        end
      end
      sP=max(abs(sum(P,2)-1));
    end
    if sP>tol
      if colstoch        
        warnings{end+1}={18,sP};
      else
        warnings{end+1}={19,sP};
      end
    end
 end
  
