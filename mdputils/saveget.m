% saveget Save a variable to disk and provides a function handle to retreive it.
% USAGE
%   fhandle=saveget(X);
% This function saves the variable X to disk in a MAT file with a randomly
%   generated name, clears the variable from the caller's workspace and 
%   creates a function that can be used to retreive it.
% Example:
%   P=randn(3,4); f=saveget(P);
% At this point P is not available. It can be reteived using 
%   x=f();
% Notice that it can be assigned to whatever variable name is convenient 
%   (in this case x).
% 
% This function leaves randomly named MAT files on you hard disk. These can
% (should) be erased after you finish you session (a way to do automatic cleanup
% would be nice).

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

function fhandle=saveget(X) %#ok<INUSD>
if nargin~=1
  error('A single variable must be passed')
end

% generate a name for a MAT file
matfilename='';
while isempty(matfilename)
  matfilename=genvarname(char(ceil(rand(1,8)*128)));
  try
    eval(['load ' matfilename]);
    matfilename='';
  catch %#ok<CTCH>
     % continue on
  end
end
% get the variable name in the caller's workspace (if it exists)
vname=inputname(1);
if isempty(vname)  % isn't a variable in the caller's workspace
  % get made up name
  vname=genvarname(char(ceil(rand(1,8)*128)));
  % assign X to the made up name
  eval([vname '=X;'])
  save(matfilename,vname);
else
  % save the variable using its original name
  evalin('caller',['save(''' matfilename ''',''' vname ''')'])
  % clears the original variable from the caller's workspace
  evalin('caller',['clear ' vname]) 
end
% create a function handle to retrieve the vaiable from disk
fhandle = @() getfield(load(matfilename,vname),vname);  %#ok<GFLD>