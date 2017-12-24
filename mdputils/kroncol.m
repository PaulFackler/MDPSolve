% kroncol Columnwise Kronecker product
% USAGE
%   C=kroncol(A,B);
% INPUTS
%   A : mxn matrix
%   B : pxn matrix
% OUTPUT
%   C : mpxn matrix
%
% C(:,j)=kron(A(:,j),B(:,j))
% 
% This is sometimes known as the Khatri-Rao product.

% Coded as a MEX file

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

function C=kroncol(A,B)
[ma,na]=size(A);
[mb,nb]=size(B);
if nb~=na
  error('Matrices must have the same number of columns')
end
if ~(issparse(A) || issparse(B))
   C=bsxfun(@times,reshape(A,[1 ma na]),reshape(B,[mb 1 nb]));
   return
end
try
  C=A(ones(mb,1)*(1:ma),:).*B((1:mb)'*ones(1,ma),:);
catch %#ok<CTCH>
  if 1 % 2 approaches - not clear which is better in terms of speed
    N=sum(sum(A~=0,1).*sum(B~=0,1));
    C=sparse([],[],[],ma*mb,na,N);
    for j=1:na
      cj=B(:,j)*A(:,j).';
      C(:,j)=cj(:);
    end
  else
    C=cell(1,na);
    for j=1:na
      C{j}=vec(B(:,j)*A(:,j).');
    end
    C=[C{:}];
  end
end