% Rouwenhorst grid and transition probability for discrete AR(1) model
% USAGE
%   [z,P]=Rouwenhorst(n,mu,rho,sigma);
% INPUTS
%   n     : number of grid points
%   mu    : long-run mean
%   rho   : AR(1) coefficient
%   sigma : noise variance
% OUTPUTS
%   z : grid vector
%   P : nxn transition probability matrix
% 
% This provides a discretization of the model
%   z(t+1) = mu + rho(z(t)-mu) + e
% where e~N(0,sigma^2)
%
% The method ensures that the discrete process has the same mean, variance and 
% AR(1) coefficient as the continuous process.

% Rouwenhorst, K. Geert, 1995. Asset Pricing Implications of Equilibrium Business Cycle Models. 
% In: Cooley, T.F. (Ed.), Frontiers of Business Cycle Research. Princeton University Press, Princeton, NJ, 294-330.

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

function [z,P]=Rouwenhorst(n,mu,rho,sigmae)
psi=sigmae*sqrt((n-1)/(1-rho^2));
z=mu+linspace(-psi,psi,n)';
p=(1+rho)/2;
q=p;
P=[p 1-p;1-q q];
for i=3:n
  P =    p*[P zeros(i-1,1); zeros(1,i)]     + ...
     (1-p)*[zeros(i-1,1) P; zeros(1,i)]     + ...
     (1-q)*[zeros(1,i)    ; P zeros(i-1,1)] + ...
         q*[zeros(1,i)    ; zeros(i-1,1) P];
  P(:,2:i-1)=P(:,2:i-1)/2;
end