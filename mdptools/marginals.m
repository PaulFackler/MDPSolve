% marginals Computes marginal distributions from a probability vector
% USAGE
%   M=marginals(P,n);
% or
%   M=marginals(P,S);
% INPUTS
%   P : Nx1 vector of probabilities (summing to 1)
%   n : d-vector of variable sizes (prod(n)=N)
%   S : Nxd matrix of variable values
% OUTPUTS
%   M : d element cell array. Element i contains an n(i)-vector of
%         probabilities associated with variable i

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

function M=marginals(P,input2)
if all(size(P)>1)
  error('P must be a vector')
end
% S matrix passed 
if size(input2,1)==numel(P)
  S=input2;
  d=size(S,2);
  M=cell(1,d);
  for i=1:d
    [Si,temp,ind]=unique(S(:,i));
    %M{i}=indexsum(P,ind,length(Si));
    M{i}=accumarray(ind,P,[length(Si) 1]);
  end
% vector of variable sizes passed
else
  n=input2;
  if numel(P)~=prod(n)
    error('inputs are incompatible')
  end
  d=length(n);
  n=flipud(n(:))';
  P=reshape(P,n);
  M=cell(1,d);
  for i=1:d
    p=P;
    for j=d:-1:1
      if j~=i
        p=sum(p,j);
      end
    end
    M{d-i+1}=full(p(:));
  end
end