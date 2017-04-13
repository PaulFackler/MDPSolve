% lpsolve Solves linear programming problems
%   min c'x
%   s.t.
%   Ax<=>b
%   x>=0
% USAGE
%  [x,s,message]=lpsolve(c,A,b,conind);
%
% conind has the same number of rows as A and b:
%   -1 for less than
%    0 for equal to
%    1 for greater than 
%
% OUTPUTS
%  x: the optimal value
%  s: a (m1+m2)x3 matrix.  Column 1 contains the slackness variable number
%     (if followed by .1) or the surplus variable number (if followed by .2)
%     Column 2 contains the value of the slackness or surplus variable,
%     column 3 contains the shadow price of the constraint.
%  Message: a message concerning the results
%            OK
%            Unbounded solution
%            No feasible solution
%            Maximum iterations exceeded

% Programming note: CONIND is converted to a 3-vector 
% with the # of less than, greater than and equal to constraints
% and the constraints are reordered (less, greater, equal)
% There are m1 less than constraints and m2 greater than constraints

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
% Edited 2011, Paul L. Fackler

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

function [x,s,message]=lpsolve(c,A,b,conind)

if nargin~=4
  error('4 arguments must be passed')
end

% Check input sizes
  [n,m]=size(A);
  if length(b)~=n;
    disp('RHS vector (b) inconsistent with constraint matrix (A)');
    disp(['LENGTH(b): ' num2str(length(b)) '  not equal to ROWS(A): ' num2str(n)])
    error(' ')
  end
  if length(c)~=m;
    disp('Objective function vector (c) inconsistent with constraint matrix (A)')
    disp(['LENGTH(c): ' num2str(length(c)) '  not equal to COLS(A): ' num2str(m)])
    error(' ')
  end

  if length(conind)~=n;
     disp('Constraint indicator vector (CONIND) inconsistent with constraint matrix (A)')
     disp(['LENGTH(CONIND): ' num2str(length(conind)) '  not equal to ROWS(A): ' num2str(n)]);
     error(' ')
  end

  indl=find(conind==-1); inde=find(conind==0); indg=find(conind==1);
  ind=[indl;indg;inde];
  conind=[length(indl);length(indg);length(inde)];
  A=A(ind,:);
  b=b(ind);

  n1=n+1;
  numless=conind(1,1);
  nummore=conind(2,1);
  numeq=conind(3,1);
  c=c(:)';
  basis=0;
  nonbasis=(1:m)';
  if numless; basis=[basis;m+(1:numless)']; end
  if nummore;
      basis=[basis;m+numless+nummore+(1:nummore)'];
      temp=sparse(n,nummore);
      temp(numless+1:numless+nummore,:)=-speye(nummore);
      A=[A temp];
      c=[c zeros(1,nummore)];
      nonbasis=[nonbasis;m+numless+(1:nummore)'];
  end
  if numeq; basis=[basis;m+numless+2*nummore+(1:numeq)']; end
  basis=basis(2:end,:);

% Adjust for artificial variables) 
if (nummore+numeq)>0
  bignum=10000;
  c=c-bignum*sum(A(numless+1:n,:));
  z=-bignum*sum(b(numless+1:n,1));
else
  z=0;
end

A=[A b(:);c z];
[basis,nonbasis,A,message]=lpx(basis,nonbasis,A,nummore+numeq);
x=zeros(m,1);
i=find(basis<=m);
% i';basis(i)';m;n;ord(a);mm1;
x(basis(i))=A(i,end);
if (numless+nummore)>0
  s=[];
  if numless>0; s=[s;(1:numless)'+0.1]; end
  if nummore;s=[s;(1:nummore)'+0.2]; end
  s=[s zeros(nummore+numless,2)];
  i=find((basis>m) & (basis<=m+numless+nummore));
  if ~isempty(i), s(basis(i)-m,2)=A(i,size(A,2));end
  i=find((nonbasis>m) & (nonbasis<=m+numless+nummore));
  if ~isempty(i), s(nonbasis(i)-m,3)=-A(n1,i)'; end
else
  s=[];
end
if nargout<3 && ~strcmp(message,'OK')
  disp(['Warning: ' message])
end


%
% Performs simplex steps with a starting basis A for which the
% basis and nonbasic variables are listed in the the vectors
% BASIS and NONBASIC.  The number of artificial variables is passed
% as the scalar ARTIF.
%
function [basis,nonbasis,A,message]=lpx(basis,nonbasis,A,artif)
maxcount=5000;
% Iterate until convergence @
count=0;
message='OK';
tol=4*eps;
[n1,mm1]=size(A);
n=n1-1;
mm=mm1-1;
m=mm+n-artif;
if artif>0, phase=0; else phase=1; end
[c,j]=min(A(n1,1:mm));
while (c<-tol) && (count<=maxcount);
  count=count+1;
  i=find(A(1:n,j)>tol);
  if isempty(i);
    message='Unbounded solution'; break;
  else
    r=A(i,mm1)./A(i,j);
    i=i(r==min(r));      
    if size(i,1)>1                       % in degenerate case
      i=i(ceil(rand(1,1))*size(i,1),:);  % pick randomly
    end
    pivot=1/A(i,j);
    tempj=A(:,j)*(-pivot);
    tempi=A(i,:)*pivot;
    A=A-A(:,j)*tempi;
    A(i,:)=tempi;
    A(:,j)=tempj;
    A(i,j)=pivot;
    tempi=basis(i);
    basis(i)=nonbasis(j);
    nonbasis(j)=tempi;
    if phase==0;
      itemp=basis>m;
      if all(itemp==0)
        phase=1;
        itemp=find(nonbasis<=m);
        nonbasis=nonbasis(itemp);
        A=A(:,[itemp;mm1]);
        mm1=size(A,2);
        mm=mm1-1;
      end % if
    end % if
    [c,j]=min(A(n1,1:mm));
  end % if
end % while

if count>maxcount, message='Maximum iterations exceeded'; end
if phase==0; message='No feasible solution'; end
% Rearrange the final tableau to aid reading 
[basis,rowind]=sort(basis);
A=A([rowind;n1],:);
[nonbasis,colind]=sort(nonbasis);
A=A(:,[colind;mm1]);

