% ifthenelse One line if/then/else function
% USAGE
%   x=ifthenelse(condition,trueres,falseres);
% Inputs
%   condition : logical array
%   trueres   : vector with x(i)=trueres(i) when condition(i)=true
%   falseres  : vector with x(i)=falseres(i) when condition(i)=false
% OUTPUT
%   x         : array of the same size as condition
%
% Note: trueres and falseres can be scalrs or arrays of equal size as
% condition
%
% Performs like the C expression c?t:f
%
% Example:
%   a=randn(3,2,2); b=2; x=ifthenelse(a<b,a,b);
% Returns the same value as x=min(a,b);

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

function x=ifthenelse(condition,trueres,falseres)
  if ~islogical(condition)
    condition=logical(condition);
  end
  if isscalar(falseres)
    x=falseres(ones(size(condition)));
  else
    x=falseres;
  end
  if isscalar(trueres)
    x(condition)=trueres;
  else
    x(condition)=trueres(condition);
  end