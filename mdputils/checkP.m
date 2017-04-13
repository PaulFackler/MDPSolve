% checkP Checks to determine if a matrix is column or row stochastic
% USAGE
%   code=checkP(P,ns,nx);
% INPUT
%   P    : an m x n matrix of real numbers
%   ns   : number of next period states
%   nx   : number of this period state/action combinations
% OUTPUT
%   code : -2) can't determine - P is square and rows and columns sum to 1
%          -1) matrix is invalid 
%           0) matrix is row stochastic
%           1) matrix is column stochastic
%
% Tests are:
%     code = -2 P is square and rows and columns sum to (approximately) 1
%     code =  1 if m=ns and n=nx and columns sum to (approximately) 1 (and not code=-2)
%     code =  0 if m=nx and n=ns and rows    sum to (approximately) 1 (and not code=-2)
%     code = -1 otherwise
%
% Note: tests for positively are not conducted
% The sum check takes anything number on (1-tol,1+tol) to equal 1
% where tol=1e-15*min(m,n)
%
% If called without assignment to a variable a message is displayed.

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

function code=checkP(P,ns,nx)
  [m,n]=size(P);
  tol=1e-15*min(m,n);
  code = -1;
  if m==n && ns==nx && m==ns     % P is square and the right size
    spc=max(abs(sum(P,1)-1));
    spr=max(abs(sum(P,2)-1));
    if     spc<tol && spr>tol,  code =  1;
    elseif spc>tol && spr<tol,  code =  0;
    elseif spc<tol && spr<tol,  
      if all(all(P==P')), code =  1;  % symmetric is okay - colstoch is arbitrary
      else                code = -2;
      end
    end
  elseif m==nx && n==ns           % should be row stochastic
    spr=max(abs(sum(P,2)-1));
    if spr<tol, code=0;  end
  elseif m==ns && n==nx           % should be column stochastic
    spc=max(abs(sum(P,1)-1));
    if spc<tol, code=1;  end
  else 
    code=-1;                      % not the right size
  end
  
  if nargout==0
    switch code
      case -2
        disp('Can''t determine whether matrix is column or row stochastic')
      case -1
        disp('Matrix is not stochastic')
      case 0
        disp('Matrix appears to be row stochastic')
      case 1
        disp('Matrix appears to be column stochastic')
    end
  end