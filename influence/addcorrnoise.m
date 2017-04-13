% addcorrnoise Adds a set of correlated noise terms to an influence diagram
%   D=addcorrnoise(D,C,n);
% INPUTS
%   D : an influence diagram structure (empty to start a new diagram)
%   C : a d x d correlation matrix
%   n : # number of nodel values for discretization (optional)
% OUTPUTS
%   D : an updated diagram with d new hidden variables and d correlated
%         noise variables

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

function D=addcorrnoise(D,C,n)
if nargin<3, n=0; end
d=size(C,1);
Q=sqrtm(C);
z=cell(1,d);
s='z1';
for i=1:d
  z{i}=['z' num2str(i)];
  D=add2diagram(D,z{i},'h',0,{},rvdef('n',[0;1],n));
  if i>1, s=[s ',' z{i}]; end
end
for i=1:d
  f=eval(['@(' s ')[' s ']*Q(:,' num2str(i) ')']);
  D=add2diagram(D,['e' num2str(i)],'c',0,z,rvdef('f',f));
end
  