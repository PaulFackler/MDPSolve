% mdpsolve_Inf_stage Solves discrete-state/action infinite horizon dynamic program with stages
% USAGE
%  results = mdpsolve_Inf_stage(R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP, ...
%         nstage,nrep,v,algorithm,relval,vanish,maxit,tol,nochangelim,print);
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

function results = mdpsolve_Inf_stage(R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP, ...
         nstage,nrep,v,algorithm,modpol,relval,vanish,maxit,tol,nochangelim,print,checks)
       
       % let relval be a real number on [0 1] to obtain the vanishing discount approach
       % or an integer between 1 and ns to use the relative value approach with 
       % relval indicating which state to replace
       
  nochangemin=5;
      
  if isempty(v) 
    vnew=zeros(ns,1); 
  else
    vnew=v(:);    
  end
  
  if algorithm =='p'
    policyit=true;
    spi=speye(ns(1));
  else
    policyit=false;
  end
  
  if numel(nrep)==1, nrep=nrep+zeros(1,nstage); end
  if numel(ns)==1,   ns=ns+zeros(1,nstage); end
  
  V=cell(1,nstage);
  X=cell(1,nstage);
  for i=1:nstage
    V{i}=zeros(ns(i),nrep(i))-inf;
    X{i}=zeros(ns(i),nrep(i))-inf;
  end
  
  % check for non-discounted problem
  discount1=true;
  for i=1:numel(d)
    if any(d{i}~=1)
      discount1=false;
      break
    end
  end
  % select appropriate algorithm
  if discount1
    % relative value method
    if relval>=1 
      if policyit
        if colstoch(1)
          spi=spi+sparse(relval,1:ns(1),1,ns(1),ns(1));
        else
          spi=spi+sparse(1:ns(1),relval,1,ns(1),ns(1));
        end
      end
    end
    % vanishing discount method
    if vanish>0
      d={vanish}; 
      tol=tol*(1-d{1}^sum(nrep))/2;
    end
  else
    relval=false;
    if iscell(d)
      md=max(d{1});
      for i=2:numel(d), md=max(md,max(d{i})); end
    else
      md=max(d);
    end
    tol=tol*(1-md)/2;
  end
 
  if print>1
    if algorithm == 'f'
        disp('Solve Bellman equation using function iteration'); 
    else
        disp('Solve Bellman equation using policy iteration'); 
    end
  end
  checks=false;  % set to false because checks don't work when
                 % stages have different numbers of state values
  if checks 
    checkstage=repmat(true,1,nstage);
  else
    checkstage=repmat(false,1,nstage);
  end
  done=false;
  numnochange=0;
  iter=0;
  nu=-inf;
  W={}; % cell array for warnings
% %%%%%%% MAIN ITERATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while iter<maxit 
    iter=iter+1; 
    change=0;
    xnochange=true;
    if policyit && nstage>1
      PP=speye(ns(1));
      RR=zeros(ns(1),1);
    end
    stage=nstage;
    while stage>0
      clear Ri Pi deltai Ixi Iexpandi   % clear stuff that might use up memory
      [errors,warnings,Ri,Pi,deltai,nsi,nxi,Ixi,Iexpandi,colstochi,EVi,Xindexedi,expandPi] =  ...
        mdp_getparams(stage,checkstage(stage),R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP);
      if ~isempty(warnings)
        W{end+1}={50,stage}; %#ok<AGROW>
        W=[W warnings]; %#ok<AGROW>  % warning 50 displays stage information
      end
      if ~isempty(errors)>0
        results=struct('iter',iter);
        results.v=[];results.Ixopt=[];results.pstar=[]; results.errors=errors; results.warnings=W;
        results.errors{end+1}={50,stage};
        return
      end
      checkstage(stage)=false;
      nai=nxi/nsi;
      rep=nrep(stage);
      while rep>0
      % update policy 
        [vnew,xnew] = valmax(vnew); 
        if policyit
          [pstar,rstar] = valpol(xnew);
          if colstoch
            PP=PP*mxv(pstar,deltai);
            RR=rstar+((deltai.*RR(:)')*pstar)';
          else
            PP=vxm(deltai,pstar)*PP;
            RR=rstar+pstar*(deltai.*RR);
          end
        end
        if ~all(abs(vnew)<inf) % NaNs or infinities in value function
          results.errors={{35,iter}};
          %results.warnings.W;
          return
        end
        if xnochange
          xnochange=all(xnew==X{stage}(:,rep));
        end
        changei = vnew-V{stage}(:,rep);
        change = max(change,max(abs(changei)));
        V{stage}(:,rep)=vnew; 
        X{stage}(:,rep)=xnew;
        rep=rep-1;
      end
      stage=stage-1;
    end
    if print>1
      if relval>=1
        fprintf ('%5i %10.1e %5i %1.4f\n',iter,change,numnochange,nu) % print progress
      else
        fprintf ('%5i %10.1e %5i\n',iter,change,numnochange) % print progress
      end
    end
    % update V{1}(:,1) using policy iteration algorithm
    if policyit
      if colstochi
        vnew = (RR'/(spi-PP))';  % update value
      else
        vnew = (spi-PP)\RR;     % update value
      end
      clear PP RR        % clear stuff that might use up memory
    end
    if relval>=1
      nu=vnew(relval);
      vnew=vnew-nu;      % normalize so v(relval)=0
    end
    if xnochange
      numnochange=numnochange+1;
    else
      numnochange=0;
    end
    % check for convergence
    if policyit
      if numnochange
        done=true; 
      end 
    else
      if (change<=tol && numnochange>=nochangemin) || numnochange==nochangelim
        done=true; 
      end 
    end   
    %V{1}(:,1)=vnew; 
    if done, break; end
  end
  % end of iteration loop
    
  % for non-discounted problems adjust the value function to
  % output the average reward function
  if discount1 
    if relval==0 && vanish>0
      nu=0;
      for i=1:nstage
        V{i}=(1-d{1}^sum(nrep))*V{i};
        nu=nu+mean(V{i}(:));
      end
      nu=nu/nstage;
    end
  end
  if iter>=maxit
    % Failure to converge 
    W{end+1}={53,maxit};
  end
  % get Pstar
  if any(EV)
    pstar=[];
  else
    pstar=speye(ns(1));
    for stage=nstage:-1:1
      [errors,warnings,Ri,Pi,deltai,nsi,nxi,Ixi,Iexpandi,colstochi,EVi,Xindexedi,expandPi] =  ...
          mdp_getparams(stage,checkstage(stage),R,P,d,ns,nx,Ix,Iexpand,colstoch,EV,Xindexed,expandP);
      for rep=nrep(stage):-1:1
        x=X{stage}(:,rep);
        if colstochi
          pstar=pstar*valpol(x);
        else
          pstar=pstar*valpol(x)';
        end
      end
    end 
  end
  % for problems with R ns x nx/ns need to transform a to x
  for i=1:nstage
    if ~Xindexed(min(i,numel(Xindexed)))
      sx=size(X{i}); sx=prod(sx(2:end));
      X{i}=ns(i)*X{i} + (1-ns(i):0)'*ones(1,sx);
    end
  end
  v=V; 
  x=X; 
  % collect problem information into results structure
  results=struct('pstar',pstar,'iter',iter,'change',change,'numnochange',numnochange);
  results.v=v;
  results.Ixopt=x;
  results.errors={};
  results.warnings=W;
  results.AR=nu;
  if policyit
    if relval>=1,    results.algorithm='prv';
    elseif relval>0, results.algorithm='pvd';
    else         results.algorithm='p';
    end
  else
    if relval>=1,    results.algorithm='frv';
    elseif relval>0, results.algorithm='fvd';
    else         results.algorithm='f';
    end
  end


% gets the maximized value function   
function [vnew,xnew] = valmax(v)
  if EV
    vnew=Pi(v);
    vnew=vnew(:);
  else
    if isa(Pi,'function_handle')  % may want to get Pi on the fly to limit memory use
      Pi=Pi();
    end
    if  colstochi, vnew=(v'*Pi)';
    else           vnew=Pi*v;
    end
    if expandPi,    vnew=vnew(Iexpandi);  end
  end
  vnew=Ri(:)+deltai.*vnew;
  if Xindexedi
    [vnew,xnew]=indexmax(vnew,Ixi,nsi);  % use mex version for greater speed
  else
    [vnew,xnew]=max(reshape(vnew,nsi,nai),[],2);
  end
end

function [pstar,rstar] = valpol(x)
  if Xindexedi
    ind=x;
  else
    ind=nsi*x+(1-nsi:0)';
  end
  rstar=Ri(ind);
  if expandPi
    ind=Iexpand(ind);
  end
  if colstochi
    pstar=Pi(:,ind);
  else
    pstar=Pi(ind,:);
  end
end

end
