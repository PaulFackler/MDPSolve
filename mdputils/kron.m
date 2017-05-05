%KRON   Kronecker tensor product.
%   KRON(X,Y) is the Kronecker tensor product of X and Y.
%   The result is a large matrix formed by taking all possible
%   products between the elements of X and those of Y.   For
%   example, if X is 2 by 3, then KRON(X,Y) is
%
%      [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%        X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
%
%   If either X or Y is sparse, only nonzero elements are multiplied
%   in the computation, and the result is sparse.
%
%   Class support for inputs X,Y:
%      all numeric and logical types
%   Result is double unless both X and Y are logical, in which case the 
%     result is logical.

% This is a replacement for the builtin MATLAB kron function
% It executes faster and uses less memory.
% Help comments are based on builtin kron function.
%
% m-file version is based on kron function 
% created by Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%   Version: 06/02/2011
% Copyright (c) 2010, Laurent Sorber
% All rights reserved.
% Sorber's function is redefined to reduce memory usage

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

function C=kron(A,B)
  [I J] = size(A);
  [K L] = size(B);
  if issparse(A) || issparse(B)
    % get indices of non-zero elements of A and B
    [ia,ja,C] = find(A);
    [ib,jb,c] = find(B);
    % combine to get indices for C
    % ia & ja are overwritten as they are no longer needed
    ia = bsxfun(@plus,K*(ia(:)-1).',ib(:));  
    ja = bsxfun(@plus,L*(ja(:)-1).',jb(:));
    clear ib jb
    % get the non-zeros elements of C
    if islogical(A) && islogical(B)
      % The @and operator is slightly faster for logicals. 
      C=nonzeros(bsxfun(@and,c(:),C(:).'));
    else
      % otherwise converts to double to ensure compatability
      C=nonzeros(double(c(:))*double(C(:)).');  
    end
    clear c 
    C = sparse(ia,ja,C,I*K,J*L);
  else
    % Both matrices are dense.
    A = reshape(A,[1 I 1 J]);
    B = reshape(B,[K 1 L 1]);
    C = reshape(bsxfun(@times,A,B),[I*K J*L]);
  end
    
