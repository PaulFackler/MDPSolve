% simplexbasc Creates a basis matrix for a grid on a simplex
% A simplex is a space of q non-negative numbers that sum to C
% A regular grid on a simplex is defined by p evenly spaced intervals
%   (implying p+1 grid points for each dimension)
% Only q-1 values is needed to specify a point on the simplex because
%   the qth value is equal to C less the sum of the first q-1 value
% USAGE
%   B=simplexbasc(x,q,p,C);
% or
%   [b,ir]=simplexbasc(x,q,p,C);
% INPUTS
%   x    : mx(q-1) matrix of evaluation points (each column is a point)
%            Can be mxq but the last column is ignored
%   q    : dimension of the simplex (positive integer)
%   p    : number of subintervals in each dimension (positive integer)
%   C    : grid values are non-negative and sum to C or less  (positive number)
% OUTPUTS
%   B    : nxm matrix of basis values 
%          where n is the number of grid points:
%             n = (p+q-1)!/p!/(q-1)!
% or
%   b    : qxm matrix of values
%   ir   : qxm matrix of row indices
% The two outputs are related by B=sparse(ir,ones(q,1)*(1:m),b,m,n)

% MEX file version of simplexbas  
% This function is called by simplexbas and is best if not called directly
% because sparsebas included error checking

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

function varargout=simplexbasc(varargin)
error('This function should not be called')
