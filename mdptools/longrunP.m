% longrunP Analyzes Markov transition probability matrices
% USAGE
%   [q,F]=longrunP(p,options);
% or
%   q=longrunP(pf,ns);
% INPUT
%   P       :  ns x ns Markov transition probability matrix
%               (non-negative matrix with unit column sums)
%   options : a structure variable defining options governing the procedure
%   pf      : a function handle that returns P*q (must accept a ns-vector that 
%               sums to 1 and return an ns-vector that sums to 1)
%   ns      : the number of rows/columns in P
% OUTPUTS
%    q : ns x k matrix of invariant distributions 
%          for each recurrence class
%        q(i,k)=long-run frequency of visits to state i
%          given the system enters recurrence class k
%    F : ns x ns matrix of accessibility probabilities
%        F(i,j)= Prob[state i will be reached from state j]
% 
% Options:
%   tol : if sum(P,1)-1 falls in [-tol,tol] it is considered to equal 0
%   fast : a generally faster method if the chain is recurrent
%            [default: fast=1]
%
% Useful in analyzing long-run characteristics of the optimized
% transition probability matrix resulting from the solution of discrete
% dynamic programming problems, e.g., the matrix PSTAR returned by
% MDPSOLVE.
%
% Reference: "Introduction to Stochastic Processes" by Ethan Cinlar
%   Prentice-Hall, 1975. Chapters 5 and 6.
%
% See also: MDPSOLVE.

% Note: same as Markov from CompEcon Toolbox except P can be transposed
% CompEcon Toolbox, (c) 1997-2000, Paul L. Fackler & Mario J. Miranda 

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2011-2014, Paul L. Fackler (paul_fackler@ncsu.edu)
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

function [q,f]=longrunP(p,options)
if nargin<2, options=[]; end
if isa(p,'function_handle')
  q=longrunPf(p,options);
  f=[];
  return
end
getopts(options, ...
 'tol',   5e-15,...     % tolerance to check if probabilities sum to 1
 'fast',  1);           % use fast method
n=size(p,1);


% Error Checking to ensure P is a valid stochastic matrix
if size(p,2)~=n
  error('Transition matrix is not square');
end
if any(abs(sum(p,2)-1)>tol*n)
  try
    if nargout==2
      [q,f]=longrunP(p',options);
    else
      q=longrunP(p',options);
    end
  catch ME
    error('Matrix does not appear to be a valid probability matrix')
  end
  return
end
if any(any(p<0))
  %warning('longrunP:pnegative','Transition matrix contains negative elements');
  disp('In longrunP: Transition matrix contains negative elements - applying cleanup');
  p(p<0)=0;
  p=vxm(1./sum(p,2),p);
end

if fast && nargout<2
  [q,flag]=getlrp(p); % try Krylov method
  if flag==0
    return
  else % iterate on the probability matrix and return
    q=ones(1,n)/n;
    for i=1:1000
      q0=q;
      q=q*p;
      if max(abs(q-q0))<1e-14; 
        q=q';
        return
      end
    end
  end
end

% Determine accessibility from i to j
f=zeros(n,n);
for j=1:n
  dr=1;
  r=spones(p(:,j));            % a vector of ones where p(i,j)~=0
  while any(dr)
    dr=r;
    r=spones(p*r+r);
    dr=r-dr;
  end
  f(:,j)=r;
end

% Determine membership in recurrence classes
% Programming note:
%  f(i,:)=1 for states accessible from i
%  f(:,i)'=1 for states from which i is accessible
%  f(:,i)'.*f(i,:)=1 for states communicating with i (two-way accessibility)
%  If the set of communicating states is the same as the set of accessible 
%    states, it forms a recurrence class.
ind=zeros(n,n);
numrec=0;                   % number of recurrence classes
for i=1:n
  if all(ind(i,:)==0)
    j=f(i,:);               % states accessible from i
    if all((f(:,i)'.*j)==j) % are all accessible states communicating states?
      j=find(j);            % members in class with state i
      k=length(j);          % # of members
      if k>0 
        numrec=numrec+1; 
        ind(j,numrec)=ones(k,1); 
      end
    end
  end
end
ind=ind(:,1:numrec);        % ind(i,j)=1 if state i is in class j

% Determine recurrence class invariant probabilities
q=zeros(n,numrec);
for j=1:numrec
  k=find(ind(:,j));             % members in class j
  nk=length(k);                 % # of members
  % solve Pq=q s.t. 1'q=1
  q(k,j)=[ones(1,nk);(speye(nk)-p(k,k)')]\[1;zeros(nk,1)];
end

% Analyze transients if second output desired
if nargout>1
  if numrec>1 
    trans=find(sum(ind,2)==0);
  else
    trans=find(ind==0);
  end
  numt=length(trans);          % number of transient states
  % Determine transient absorption and reachability probabilities
  if numt>0                    
    qq=p(trans,trans);
    b=zeros(numt,n);
    for j=1:numrec
      k=find(ind(:,j));        % members of class j
      nk=length(k);            % # of members
      if nk==1                 % 1-step prob: transient states to class j
        b(:,k)=p(trans,k);
      else
        b(:,k)=sum(p(trans,k),2)'*ones(1,nk);
      end
    end;
    qq=inv(eye(numt)-qq);          
    f(trans,:)=qq*b;  %#ok<MINV>     % absorption probabilities
    d=diag(qq)';
    qq=qq./d(ones(numt,1),:);
    f(trans,trans)=qq-diag(1./d);    % transient reachability probabilities
  end
  f=f';
end
end

% getlrp Long run distribution of a Markov chain using Krylov methods
% USAGE
%   p=getlrp(P);
% INPUT
%   P    : n x n transition probability matrix (column stochastic)
% OUTPUTS
%   p    : n x 1 long run (stationary) distribution
%   flag : flag from bicgstab (0 is no problems detected)
function [p,flag]=getlrp(P)
tol=1e-8;
maxit=100;
n=size(P,1);
[p,flag]=bicgstab(@getres,[1;zeros(n-1,1)],tol,maxit,[],[],ones(n,1)/n);
p=p/sum(p);
if nargout==1 && flag ~=0
  warning('check convergence')
end

function res=getres(p)
  res=p-P'*p;
  res(1)=sum(p);
end



end

% longrunPf Analyzes Markov transition probability matrices
% USAGE
%   q=longrunPf(p,options);
% INPUT
%   p       : a function handle mappng R^n to R^n; p(x)=P*x
%               where P is a transition probability matrix
%   options : a structure variable defining options governing the procedure
% OUTPUTS
%   q : n x 1 vector containing the invariant distribution of the process
%        (assumes this exists and is unique)
%        q(i)=long-run frequency of visits to state i
function plr=longrunPf(p,options)
if isnumeric(options)
  n=options;
  maxit=1000;
  tol=1e-12;
else
end
p0=ones(n,1)/n;
plr=p(p0);
if abs(sum(plr)-1)>1e-14
  error('transition probability function incorrectly specified') 
end
for i=1:maxit
  p0=plr;
  plr=p(plr);
  r=norm(p0-plr);
  if r<tol, break; end
end
if i>=maxit
  fprintf('failure to converge in longrunPf: residual = %1.4e\n',r)
end
end
