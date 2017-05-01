% pomdp Converts Partially Observable MDPs to standard form
% USAGE
%   [b,Pb,Rb]=pomdp(p,P,Q,R,options);
% INPUTS
%   p        : number of belief states for each state (scalar)
%   P        : (ns x ns x na) array of state transition probabilities
%                can be (ns x (ns*na)) with na blocks each of size ns x ns stacked horizontally
%   Q        : (ny x ns x na) array of conditional observation probabilities
%                can be (ny x (ns*na)) with na blocks each of size ny x ns stacked horizontally
%   R        : (na x ny) or (nx x 1) or (nx x ny) array of rewards (see below)
%   options  : structure variable to control procedure (see below)
% OUTPUTS
%   b        : belief nodal values [N x ns]
%   Pb       : belief state transition matrix
%   Rb       : belief state reward matrix
%
% In Partially Observed Markov Decision Problems (POMDPs) 
% the state S is not directly observed. Instead Y is observed, where 
%   Prob(Y|S+,A)=Q  (set options.Qtype=1 [default])
% or
%   Prob(Y|S,A)=Q  (set options.Qtype=0)
% The first interpretation gives the signal after the state transition,
% the second gives the signal before the state transition.
% options.Qtype:
%    0: Q depends on current action and state
%    1: Q depends on current action and future state
%
% Three possible forms of the reward are handled by POMDP which are distinguished
%   by setting Rtype (this is only necessary if there is potential for ambiguity)
% options.Rtype
%    1: R (ny x na) depends on action and observation
%    2: R (nx x 1)  depends on current state/action
%    3: R (ny x nx) depends on current state/action and observation
% Note that in cases 2 and 3 the reward in unobserved
%
% If A is an na row matrix of the state action combinations under observational
%   certainty, the state/action matrix for the augmented problem can be obtained
%   using
%     Ab=rectgrid(A,b);
% The optional actions are then given by 
%     Aopt=Ab(results.Ixopt,1:size(A,2))

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

function [b,Pb,Rb,omega]=momdp(p,no,nu,P,Q,R,options)
if ~exist('options','var'), options=[]; end
getopts(options, ...
 'Qtype',   1, ...       % Q conditional on (0) current state (1) future state
 'Rtype',   []);         % R a function of (1) A,Y  (2) S,A  (3) S,A,Y

[b,Pb,Rb,omega]=pomdprun(p,no,nu,P,Q,R,Qtype,Rtype);
end

function [b,Pb,Rb,omega]=pomdprun(p,no,nu,P,Q,R,Qtype,Rtype)
[ns,nx]=size(P);
na=nx/ns;
if na~=round(na)
  error('P has improper size');
end
ny=numel(Q)/nx;
Q=reshape(full(Q),[ny,ns,na]);

mQ=abs(sum(Q,1)-1); 
if any(mQ(:)>1e-14)
  error('Q is not a proper probability array')
end

% determine Rtype (if needed) from the size of R
if isempty(Rtype)
    if all(size(R)==[ny,na])
      Rtype=1;
    elseif numel(R)==nx && any(size(R)==1)
      Rtype=2;
    elseif all(size(R)==[ny,nx])
      Rtype=3;
    else
      error('The size of R is not consistent with other input data')
    end
end

b=kron(simplexgrid(nu,p,1),simplexgrid(no,1,1));
if nargout==1, return; end

N=size(b,1);
% Get Bayesian updates for belief probabilities for each observation (tau)
% and conditional probabilities for observations
Pb=getPb;
if nargout==2, return; end
clear tau

% belief state reward matrix
switch Rtype
  case 1 % R is ny x na
    Rb=zeros(N,na);
    for j=1:na
      Rb(:,j)=reshape(omega(:,j,:),[N,ny])*R(:,j);
    end
  case 2 % R is nx x 1 (or ns x na)
    Rb=b*reshape(R,ns,na);
  case 3 % R is ny x nx
    Rb=b*reshape(sum(reshape(R,[ny,ns,na]).*Q,1),ns,na);
  otherwise
    error('Improper selection for Rtype')
end


% compute the transition probability matrix
function Pb=getPb
% This operation expands P and Q, orders both jkli, and multiplies them
if Qtype % assumes Y depends on S+
  % numerator of tau is sum_i P(ijk)Q(ikl)b(j)
  tau=repmat(permute(reshape(P,[ns ns na 1]),[2 3 4 1]),[1 1 ny 1]).* ...
      repmat(reshape(permute(Q,[3 1 2]),[1 na ny ns]),[ns 1 1 1]);
else     % assumes Y depends on S
  % numerator of tau is sum_j P(ijk)Q(jkl)b(j)
  Q=permute(Q,[2 3 1]);
  tau=repmat(permute(reshape(P,[ns ns na 1]),[2 3 4 1]),[1 1 ny 1]).* ...
    repmat(reshape(Q,[ns na ny 1]),[1 1 1 ns]);
end
% denominator of tau is sum of tau over i (future state)
taus=reshape(b*reshape(sum(tau,4),ns,na*ny),[N,na,ny]);  % sum over future states and multiply by beliefs
tau=reshape(b*reshape(tau,ns,na*ny*ns),[N na ny ns]); % updated beliefs

% observation probability given beliefs and actions
% omega(j,k,l)=probability Y=l given a=k and b=b(j) 
if Qtype
  omega=reshape(taus,[N,na,ny]);
else
  omega=reshape(b*reshape(Q,ns,na*ny),[N,na,ny]);
end
taus(taus==0)=1;  % These are action/observation pairs that can't be obtained.
                  % As the probability of these combinations is 0,
                  % they do not figure into the computations and this
                  % value is arbitrary but must be non-zero (to avoid NaNs).
% make tau sum to 1 over dimension 4 (future state)
for i=1:ns
  tau(:,:,:,i)=tau(:,:,:,i)./taus;
end
clear taus
%tau=squeeze(sum(reshape(tau,N*na,ny,no,nu),3));
tau=squeeze(sum(reshape(tau,N*na,ny,no,nu),3));
try % processing current states and actions simultaneously
  omega=reshape(omega,N*na,ny);
  Pb=sparse(N,N*na);
  ii=kron(ones(1,N*na/no),flipud(eye(no)));
  for k=1:ny
    
    Pb=Pb+kroncol(mxv(simplexbas(squeeze(tau(:,k,:)),nu,p,1),omega(:,k)),ii);
  end
  omega=reshape(omega,[N,na,ny]);
catch    %#ok<CTCH>      % otherwise process each action/y value in a loop
  disp('error encountered - switching to alternative algorithm')
  tau=reshape(tau,N,na,ny,no,nu);
  omega=reshape(omega,[N,na,ny]);
  Pb=sparse([],[],[],N,N*na);
  col=1;
  for k=1:na
    Pk=sparse([],[],[],N,N);
    for l=1:ny
      Pk=Pk+mxv(simplexbas(squeeze(tau(:,k,l,:)),nu,p,1),omega(:,k,l));
    end
    Pb=add2sparse(Pb,Pk,col,0.6,true);
    clear Pk
    col=col+N;
  end
end
end

end
