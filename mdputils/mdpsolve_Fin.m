% mdpsolve_Fin  Solves discrete-state/action finite horizon dynamic program
% USAGE
%   results = mdpsolve_Fin(R, P, d, ns, nx, Ix, Iexpand, colstoch, EV, ...
%                           Xindexed,expandP, T, v, keepall, print);
%
% Called by mdpsolve

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

function results = mdpsolve_Fin( ...
   R, P, d, ns, nx, Ix, Iexpand, colstoch, EV,Xindexed,expandP, T, v, keepall, print)

  if keepall
    vv=zeros(ns,T);
    xx=zeros(ns,T);
  end
 
  if print>1
    disp('Solve Bellman equation using backward iteration'); 
  end
 
  na=nx/ns;
  R=R(:);   % normalize so R is a column vector
  iter=T+1;
%%%%%%%% MAIN ITERATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while iter>1
    iter=iter-1;  
    [v,x] = valmax(v); 
    % this should never happen if P is proper and R is bounded
    if ~all(abs(v)<inf) % NaNs or infinities in value function
      results.errors={{35,iter}};
      return
    end
    if keepall
      vv(:,iter)=v;
      xx(:,iter)=x;
    end
  end
  % end of iteration loop
  
  if keepall
    v=vv; x=xx;
  end
    
  % for problems with R ns x na need to transform a to x
  if ~Xindexed
    x=ns*x + (1-ns:0)'*ones(1,size(x,2));
  end
  % collect problem information into results structure
  results=struct('v',v,'Ixopt',x);
    

% gets the maximized value function and associated strategy   
function [vnew,xnew] = valmax(v)
  if EV
    vnew=P(v);
    vnew=vnew(:);
  else
    if colstoch, vnew=(v'*P)';
    else         vnew=P*v;
    end
    if expandP,  vnew=vnew(Iexpand);  end
  end
  vnew=R+d.*vnew;
  if Xindexed
    [vnew,xnew]=indexmax(vnew,Ix,ns);  % use mex version for greater speed
  else
    [vnew,xnew]=max(reshape(vnew,ns,na),[],2);
  end
end

end
