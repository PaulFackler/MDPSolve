% Called by tprodm
% This algorithm works by permuting x and y so the matched dimensions are
% in the columns and the unmatched are in the rows. It then uses a
% columnwise Kronecker product. Any summed dimensions are then summed
% using a postmultiplication operator. The result is then permuted to have the
% appropriate dimensions.

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

function [z,nzz]=tprods(x,x2z,nx,y,y2z,ny,nzout,zsparse)
  if nargin<8 || isempty(zsparse), zsparse=true; end
  if nargin<7, nzout=[]; end
  if numel(x)~=prod(nx), error('x and nx are incompatible'); end
  if numel(y)~=prod(ny), error('y and ny are incompatible'); end
  d=max(abs([x2z y2z]));
  nz=ones(1,d);
  for i=1:d 
    ii=find(abs(x2z)==i); if ~isempty(ii), nz(i)=nx(ii); end
    ii=find(abs(y2z)==i); if ~isempty(ii), nz(i)=ny(ii); end
  end
  % q has the matched dimensions in sorted order
  [q,ix,iy]=intersect(x2z,y2z);
  nomatchx=find(~ismember(x2z,q));
  nomatchy=find(~ismember(y2z,q));
  % put summed dimensions (q<0) at the end
  nsum=sum(q<0);
  nmatch=length(q)-nsum;
  if nsum>0 && nmatch>0, 
     q=[ q(nsum+1:end)  q(nsum:-1:1)];
    ix=[ix(nsum+1:end) ix(nsum:-1:1)];
    iy=[iy(nsum+1:end) iy(nsum:-1:1)];
  end
  
  nxx=[prod(nx(nomatchx)) prod(nx(ix))];
  nyy=[prod(ny(nomatchy)) prod(ny(iy))];
  nm=prod(nx(ix(1:nmatch)));        % size of matched dimensions
  ns=prod(nx(ix(nmatch+1:end)));    % size of summed dimensions
  ix=[nomatchx ix(:)'];
  iy=[nomatchy iy(:)'];
  
  x=sppermute(x,ix,nx,nxx);
  y=sppermute(y,iy,ny,nyy);
  z=kroncol(y,x);   
  
  if nsum>0  % sum dimensions with negative indices  
    if nmatch>0 
      Q=kron(ones(ns,1),speye(nm));
      z=z*Q;
    else
      z=sum(z,2);
    end
  end
  
  % permute z to have the desired order
  zvar=[x2z(nomatchx) y2z(nomatchy) q];
  nzz=nz(abs(zvar));
  nzz(zvar<0)=1;
  zvar=abs(zvar);
  iz(zvar)=1:length(zvar);
  iz(iz==0)=[];
  if zsparse 
    if isempty(nzout), nzout=[nzz(iz(1)) prod(nzz(iz(2:end)))]; end
    if any(diff(iz)<0)
      z=sppermute(z,iz,nzz,nzout);
    else
      z=reshape(z,nzout);
    end
  else 
    if isempty(nzout), nzout=nzz(iz); end
    if any(diff(iz)<0)
      z=permute(reshape(full(z),nzz),iz);
    end
    z=reshape(z,nzout);
  end
 return
 
 function [ind,vals]=ismember(a,b)
   m=length(a);
   ind=false(1,m);
   for i=1:m
     if any(a(i)==b)
       ind(i)=true;
     end
   end
   if nargout>1
     vals=a(ind);
   end
   
 function  [c,ia,ib]=intersect(a,b)
   a=a(:);
   b=b(:)';
   x=bsxfun(@eq,a,b);
   ia=find(any(x,2));
   [c,ind]=sort(a(ia));
   ia=ia(ind)';
   ib=find(any(x,1));
   [c,ind]=sort(b(ib));
   ib=ib(ind);   
   return