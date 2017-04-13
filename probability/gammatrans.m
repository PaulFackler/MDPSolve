% gammatrans Gets the parameters of the Gamma distribution from the sufficient statistics
% USAGE
%  [alpha,lambda]=gammatrans(M1,M2);
%  theta=gammatrans(M);
% INPUTS
%    M1 : n-vector of values of E[X]    
%    M2 : n-vector of values of E[ln(X)]
%    M  : n x 2 matrix of [M1 M2]
% OUTPUTS
%    alpha  : n-vector of values of the shape parameter
%    lambda : n-vector of values of the scale parameter
%    theta  : n x 2 matrix of [alpha,lambda]
%
% Note: log(M1)>M2 (a NaN is returned if this is not true)
%
%  X ~ cx^(alpha-1)exp(-x/lambda) (c is the constant of integration)
%
% Uses Newton's method on ln(alpha):
%   ln(alpha)-psi(alpha)-(ln(M1)-M2)=0
%
% Residual accuracy is about 1e-13

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

function [alpha,lambda]=gammatrans(M1,M2)
tol=1e-11;
maxit=200;
if nargin==1
  s=log(M1(:,1))-M1(:,2);
else
  s=log(M1)-M2;
end
ind=s>0;
s=s(ind);
if ~isempty(s)  
  alpha=log(0.5./s);
  for i=1:maxit
    ealpha=exp(alpha);
    res=alpha-psi(ealpha)-s;
    da=res./(1-psi(1,ealpha).*ealpha);
    alpha=alpha-da;
    if max(abs(res))<tol; break; end
  end
  alpha=ealpha;
else
  alpha=[];
end
if ~all(ind)
  x=alpha;
  alpha=NaN(size(ind));
  alpha(ind)=x;
end
if nargout<=1
  alpha=[alpha M1./alpha];
else
  lambda=M1./alpha;
end
disp(i)