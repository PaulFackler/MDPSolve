% mdpsolve_Inf Solves discrete-state/action infinite horizon dynamic program
% USAGE
%  results = mdpsolve_Inf(R,P,d,ns,nx,Ix,Iexpand,colstoch,EV, ...
%         Xindexed,expandP,v,algorithm,relval,maxit,tol,nochangelim,print);
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

function results = mdpsolve_Inf(R,P,d,ns,nx,Ix,Iexpand,colstoch,EV, ...
         Xindexed,expandP,v,algorithm,modpol,relval,vanish,maxit,tol,nochangelim,print,vknown)
  warnings={};
  nochangemin=5;   % safety feature to prevent early convergence with function iteration
  
  if isempty(vknown)
    constrained=false;
  else
    constrained=true;
    iknown=~isnan(vknown);
    vknown=vknown(iknown);
  end
  
  if algorithm =='p'
    policyit=true;
    spi=speye(ns);
  else
    policyit=false;
  end
  
  
  
  % check for non-discounted problem and select appropriate algorithm
  if all(d==1)
    discount1=true;
    if algorithm=='p' && relval==0 && vanish==0
      relval=1;   % need to add a warning message
    end
    if vanish>0, d=vanish; end
  else
    discount1=false;
    tol=tol*(1-max(d))/max(d);
  end
  ktol=tol;
  % relative value method
  if relval>=1 && policyit
    if colstoch
      spi=spi+sparse(relval,1:ns,1,ns,ns);
    else
      spi=spi+sparse(1:ns,relval,1,ns,ns);
    end
  end
  
  if print>1
    if algorithm == 'f'
        disp('Solve Bellman equation using function iteration'); 
    else
        disp('Solve Bellman equation using policy iteration'); 
    end
  end
  
  if algorithm=='i'
    if 0
    if EV
      if nargin(P)==1
        selfun=@(x,i) x(i);
        if relval>=1
          EVitsol=@(V,ind) V+V(ind(relval))-d.*selfun(P(V),ind);
        else
          EVitsol=@(V,ind) V-d.*selfun(P(V),ind);
        end
      else
        if relval>=1
          EVitsol=@(V,ind) V+V(ind(relval))-d.*P(V,ind);
        else
          EVitsol=@(V,ind) V-d.*P(V,ind);
        end
      end
    else
      EVitsol =@(V,P) V-d.*(P*V);
      EVitsolt=@(V,P) V-d.*(V'*P)';
    end
    end
    EVitsol=getfunc(P,d,colstoch);
    itsol=true;
    nochangemin=2;
    warning('off','MATLAB:bicgstabl:tooSmallTolerance');
  else
     itsol=false;
  end
  
  if EV, modpol=0; end  % not tested for discount1 option yet
  %if EV, modpol=0; end  % not tested for discount1 option yet
  MPI=0;                       % counts the number of modified policy iterations 
  x=zeros(ns,1);               % initialize so always do at least 2 iterations
  na=nx/ns;
  R=R(:);                      % normalize so R is a column vector
  numnochange=0;
  done=false;
  iter=0;
  nu=[];
% %%%%%%% MAIN ITERATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while iter<maxit     
    iter=iter+1; 
    % update policy 
    [vnew,xnew] = valmax(v); 
    % update value if policy iteration is used
    if policyit 
      [pstar,rstar] = valpol(xnew);
      if colstoch 
        vnew = (rstar'/(spi-mxv(pstar,d)))';  % update value
      else
        vnew = (spi-vxm(d,pstar))\rstar;      % update value
      end
      if constrained, vnew(iknown)=vknown; end
    elseif itsol
      if EV
        if Xindexed
          ind=xnew;
        else
          ind=ns*xnew+(1-ns:0)';
        end
        rstar=R(ind); rstar=rstar(:);
        if expandP
          ind = Iexpand(ind);
        end
        [vnew,ISflag] = ...
          bicgstabl(@(V) EVitsol(V,ind),rstar,ktol,[],[],[],v);  % update value
      elseif colstoch 
         [pstar,rstar] = valpol(xnew);
         [vnew,ISflag] = ...
            bicgstabl(@(V)EVitsol(V,pstar),rstar,ktol,[],[],[],v);  % update value
      else
        [pstar,rstar] = valpol(xnew);
        [vnew,ISflag] = ...
          bicgstabl(@(V)EVitsol(V,pstar),rstar,ktol,[],[],[],v);  % update value
      end
      if constrained, vnew(iknown)=vknown; end
    % modified policy iteration
    elseif modpol>0
      if EV
        if Xindexed
          ind=xnew;
        else
          ind=ns*xnew+(1-ns:0)';
        end
        rstar=R(ind); rstar=rstar(:)';
        for k=1:modpol
          vv=vnew;
          if relval>=1, nu=vnew(relval); vnew=vnew-nu; end
          if nargin(P)==1
            vnew=EV(vnew); vnew=vnew(:);
            vnew = rstar+d.*vnew(ind)';     % update value
          else
            vnew=EV(vnew,ind); vnew=vnew(:)';
            vnew = rstar+d.*vnew;     % update value
          end
          if constrained, vnew(iknown)=vknown; end
          MPI=MPI+1;
          if max(abs(vnew-vv)./max(1,abs(vv)))<1e-12, break; end
        end
      else
        [pstar,rstar] = valpol(xnew);
        if colstoch 
          for k=1:modpol
            vv=vnew;
            if relval>=1, nu=vnew(relval); vnew=vnew-nu; end
            vnew = rstar+d.*(vnew'*pstar)';     % update value
            if constrained, vnew(iknown)=vknown; end
            MPI=MPI+1;
            if max(abs(vnew-vv)./max(1,abs(vv)))<1e-12, break; end
          end
        else
          for k=1:modpol
            vv=vnew;
            if relval>=1, nu=vnew(relval); vnew=vnew-nu; end
            vnew = rstar+d.*(pstar*vnew);      % update value
            if constrained, vnew(iknown)=vknown; end
            MPI=MPI+1;
            if max(abs(vnew-vv)./max(1,abs(vv)))<1e-12,  break; end
          end
        end
      end
    end
    if ~all(abs(vnew)<inf)      % NaNs or infinities in value function
      results.errors={{35,iter}};
      return
    end
    % make adjustments if relative value algorithm used with discount1 problems
    if relval>=1
      nu=vnew(relval);    % value of the average reward
      vnew=vnew-nu;       % adjust so vnew(relval)=0
    end
    % check if policy has changed
    xnochange=all(xnew==x);
    if xnochange, numnochange=numnochange+1;
    else          numnochange=0;
    end
    change=vnew-v;
    span = max(change)-min(change);
    change=max(abs(change));
    % check for convergence
    if policyit
      if xnochange, done=true; end
    else
      if (change<=tol && numnochange>=nochangemin) || numnochange==nochangelim
        done=true; 
      end 
    end
    if print>1
      if relval>=1
        fprintf ('%5i %10.1e %10.1e %5i %10.4f\n',iter,change,span,numnochange,nu) % print progress
      else
        fprintf ('%5i %10.1e %10.1e %5i\n',iter,change,span,numnochange) % print progress
      end
    end
    v=vnew;
    x=xnew;
    if done, break; end
  end
  % end of iteration loop
  if iter>=maxit
    % Failure to converge 
    warnings{end+1}={53,maxit};
  end
  if EV
    pstar=[];
  elseif ~policyit 
    pstar = valpol(x);
    if ~colstoch, pstar=pstar'; end   % always return pstar in column stochastic form
  end
  % for problems with R ns x na need to transform a to x
  if ~Xindexed
    x=ns*x + (1-ns:0)';
  end
  % for vanishing discount method that doesn't use the relative value approach
  % adjust the value function to approximate the discounted average reward
    if discount1 && vanish>0 && relval==0, 
      v=v.*(1-d); 
      nu=mean(v);
    end
    
  % collect problem information into results structure
  results=struct('v',v,'AR',nu,'Ixopt',x,'pstar',pstar,'iter',iter,'MPI',MPI,'change',change,...
                 'numnochange',numnochange);
  results.warnings=warnings;
  % return the algorithm used
  if relval>=1, results.algorithm=[algorithm 'rv']; 
  else      results.algorithm=algorithm; 
  end
    
% gets the maximized value function
function [vnew,xnew] = valmax(v)
  if EV
    vnew=P(v);
    vnew=vnew(:);
  else
    if colstoch, vnew=P'*v;
    else         vnew=P*v;
    end
  end
  if expandP,  vnew=vnew(Iexpand);  end
  vnew=R+d.*vnew;
  if Xindexed
    [vnew,xnew]=indexmax(vnew,double(Ix),double(ns));  % use mex version for greater speed
  else
    [vnew,xnew]=max(reshape(vnew,ns,na),[],2);
  end
end

function [pstar,rstar] = valpol(x)
  if Xindexed
    ind=x;
  else
    ind=ns*x+(1-ns:0)';
  end
  rstar=R(ind);
  if expandP
    ind=Iexpand(ind);  
  end
  if colstoch
    pstar=P(:,ind);
  else
    pstar=P(ind,:);
  end
end

end

function f=getfunc(P,d,colstoch)
    if isa(P,'function_handle')
      if nargin(P)==1
        selfun=@(x,i) x(i);
        f=@(V,ind) V-d.*selfun(P(V),ind);
      else
        f=@(V,ind) V-d.*P(V,ind);
      end
    else
      if colstoch
        f=@(V,P) V-d.*(V'*P)';
      else
        f =@(V,P) V-d.*(P*V);
      end
    end
end
