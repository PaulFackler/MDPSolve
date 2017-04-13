% betatrans Transform the sufficient statistics of the Beta distribution to usual form
% f(x;a,b) is proportional to x^(a-1)*(1-x)^(b-1)
% eta=E[[ln(X);ln(1-X)]] (note: exp(eta1)+exp(eta2)<1)
% USAGE
%   [a,b]=betatrans(eta1,eta2,method);
% or
%   [ab]=betatrans(eta,[],method);

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

function [theta1,theta2]=betatrans(eta1,eta2,method,start1,start2)
  global failcount
  failcount=0;
if nargin<3 || isempty(method), method=1; end
% organize the input to have two vectors
% for eta1 and eta2
  if nargin==1 || isempty(eta2)
    if size(eta1,1)==2
      eta2=eta1(2,:); eta1=eta1(1,:);
    else
      eta2=eta1(:,2); eta1=eta1(:,1);
    end
  end
  % get output vectors
  n=size(eta1);
  theta1=zeros(n);
  theta2=zeros(n);
  % loop over all input values
  n=max(n);
  switch method
    case 1
      optset('newton','tol',1e-12)
      getf=@getthetaN;
    case 2
      getf=@getthetaF;
  end 
  warning off MATLAB:singularMatrix
  warning off MATLAB:illConditionedMatrix
  warning off MATLAB:nearlySingularMatrix
  if nargin<4
    for i=1:n
      [theta1(i),theta2(i)]=getf(eta1(i),eta2(i));
    end
  else
    for i=1:n
      [theta1(i),theta2(i)]=getf(eta1(i),eta2(i),start1(i),start2(i));
    end 
  end
  % organize the output
  if nargout==1
    if size(theta1,1)==1
     theta1=[theta1;theta2];
    else
      theta1=[theta1 theta2];
    end
  end
  warning on MATLAB:singularMatrix
  warning on MATLAB:illConditionedMatrix
  warning on MATLAB:nearlySingularMatrix 
  %fprintf('failcount: %1i\n',failcount)
  
% uses Newton's method with log transform to avoid negative values
function [outvar1,outvar2]=getthetaN(eta1,eta2,start1,start2) 
  global failcount
  eta=[eta1 eta2];
  expeta=exp(eta);
  if sum(expeta)>=1, outvar1=NaN; outvar2=NaN; return; end
  % starting values
  %theta= eta*[0.3541 -0.1446;-0.1446 0.3541];
  if nargin<3
    theta=[1 expeta sum(expeta)^4]* ...
     [-1.3532 -1.3532;1.3974 -0.7118;-0.7118 1.3974;3.2078 3.2078];
     start1=theta(1); start2=theta(2);
  else
    theta=log([start1 start2]);
  end
  for i=1:1000
    etheta=exp(theta);
    pp=[etheta sum(etheta)];
    pp0=psi(pp);
    res=eta-pp0(1:2)+pp0(3);
    if max(abs(res))<1e-8 || any(isnan(theta)), break; end
    pp1=psi(1,pp);
    J=diag(etheta)*(pp1(3) - diag(pp1(1:2)));
    theta=theta-res/J;
  end
  i
  if i>=1000
    %disp('maximum iterations exceeded')
    outvar1=NaN; outvar2=NaN;
    failcount=failcount+1;
  end
  if any(isnan(theta))
    theta=newton(@(x)betares(x,eta'),[start1;start2]);
    outvar1=theta(1);
    outvar2=theta(2);
    if any(isnan(theta))
      [outvar1,outvar2]=getthetaF(eta1,eta2,start1,start2);
    end
    if any(isnan(theta))
      failcount=failcount+1;
    end
  else
    theta=exp(theta);
    outvar1=theta(1);
    outvar2=theta(2);
  end
  
 
  function [res,J]=betares(theta,eta)
    if any(theta<0), res=[1e100;1e100]; J=eye(2); return; end
    pp=[theta' sum(theta)];
    pp0=psi(pp);
    pp1=psi(1,pp);
    res=eta-pp0(1:2)'+pp0(3);
    J=pp1(3) - diag(pp1(1:2));
  
  

% use fminsearch
function [theta1,theta2,fval,flag]=getthetaF(eta1,eta2,start1,start2)
  options=struct('MaxIter',1000);
  if nargin<3
    expeta=exp([eta1 eta2]);
    theta=[1 expeta sum(expeta)^4]* ...
     [-1.3532 -1.3532;1.3974 -0.7118;-0.7118 1.3974;3.2078 3.2078];
  else
    theta=[start1 start2];
  end
  [theta,fval,flag,junk]=fminsearch(@(theta) f_klmin_fmin(theta,eta1,eta2),theta,options);
  theta1=theta(1);
  theta2=theta(2);
  
function F = f_klmin_fmin(x,intlnsbar,intlnomsbar)
% klmin: Kullback-Leibler divergence minimization
% Specify FOCs for the Kullback-Leibler minimization:
%F = [psi(x(1)) - psi(x(1)+x(2)) - intlnsbar;
%     psi(x(2)) - psi(x(1)+x(2)) - intlnomsbar];
 
if any(x<0)  % penalty for non-feasible levels of x
    F=[inf;inf];
else
    cc = psi([x sum(x)]);
    F = [cc(1) - cc(3) - intlnsbar;
         cc(2) - cc(3) - intlnomsbar];
end
F=F'*F;  % for fminsearch must square the output and combine.

  
  
