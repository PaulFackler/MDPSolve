% mdpsolve_Fin  Solves discrete-state/action finite horizon dynamic program
% USAGE
%   results = mdpsolve_Fin_stage(nstage,nrep, ...
%         R, P, d, ns, nx, Ix, Iexpand, colstoch, EV, Xindexed,expandP, T, v, keepall, print)
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

function results = mdpsolve_Fin_stage(nstage,nrep, ...
         R, P, d, ns, nx, Ix, Iexpand, colstoch, EV, Xindexed,expandP, T, v, keepall, print, checks,debug)
      
  W={}; % cell array for warnings
  
  if numel(nrep)==1, nrep=nrep+zeros(1,nstage); end
  if numel(ns)  ==1,   ns=ns  +zeros(1,nstage); end 
  vi=v;
  
  v=cell(1,nstage);
  x=cell(1,nstage);
  if keepall
    for i=1:nstage
      v{i}=zeros(ns(i),nrep(i),T);
      x{i}=zeros(ns(i),nrep(i),T);
    end
  else
    for i=1:nstage
      v{i}=zeros(ns(i),nrep(i));
      x{i}=zeros(ns(i),nrep(i));
    end
  end
 
  if print>1
    disp('Solve Bellman equation using backward iteration'); 
  end
   
  if checks, checkstage=repmat(true, 1,nstage);
  else       checkstage=repmat(false,1,nstage);
  end
  iter=T+1;  % work backwards from terminal date
% %%%%%%% MAIN ITERATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while iter>1
    iter=iter-1; 
    stage=nstage;
    while stage>0
      rep=nrep(stage);
      clear Ri Pi deltai Ixi Iexpandi   % clear stuff that might use up memory
      [errors,warnings,Ri,Pi,di,nsi,nxi,Ixi,Iexpandi,colstochi,EVi,Xindexedi,expandPi] =  ...
        mdp_getparams(stage,checkstage(stage),R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP,debug);
      if ~isempty(warnings)
        W=[W {50,stage} warnings]; %#ok<AGROW>  % warning 50 displays stage information
      end
      if ~isempty(errors)>0
        results=struct('iter',iter,'stage',stage);
        results.errors=errors; results.warnings=W;
        return
      end
      checkstage(stage)=false;
      nai=nxi/nsi;
      Ri=Ri(:);
      while rep>0
        [vi,xi] = valmax(vi); 
        % this should never happen if P is proper and R is bounded
        if ~all(abs(vi)<inf)           % NaNs or infinities in value function
          results.errors={{35,iter}};
          results.warnings=W;
          return
        end
        % store the results
        if keepall
          v{stage}(:,rep,iter)=vi;
          x{stage}(:,rep,iter)=xi;
        else
          v{stage}(:,rep)=vi;
          x{stage}(:,rep)=xi;
        end
        rep=rep-1;
      end
      stage=stage-1;
    end
  end
  % end of iteration loop
    
  % for problems with R ns x na need to transform a to x
  for i=1:nstage
    if ~Xindexed(min(i,numel(Xindexed)))
      sx=size(x{i}); sx=prod(sx(2:end));
      x{i}=ns(i)*x{i} + (1-ns(i):0)'*ones(1,sx);
    end
  end
  % collect problem information into results structure
  results.v=v;
  results.Ixopt=x;
  results.warnings=W;
    

% gets the maximized value function    
function [vnew,xnew] = valmax(v)
  if EVi
    vnew=Pi(v);
    vnew=vnew(:);
  else
    if colstochi, vnew=(v'*Pi)';
    else          vnew=Pi*v;
    end
    if expandPi,  vnew=vnew(Iexpandi);  end
  end
  vnew=Ri+di.*vnew;
  if Xindexedi
    [vnew,xnew]=indexmax(vnew,Ixi,nsi);  % use mex version for greater speed
  else
    [vnew,xnew]=max(reshape(vnew,nsi,nai),[],2);
  end
end

end
