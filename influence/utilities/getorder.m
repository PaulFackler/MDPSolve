% getorder Obtains an ordering for the sum-product algorithm
%   order=getorder(V,n,sumvar,penalty,orderalg);
% INPUTS
%   V        : mxn logical matrix with V(i,j)=true if factor j is a function of
%                variable i
%   n        : n-vector of variable sizes
%   sumvar   : logical n-vector with ith element true if variable can be
%                summed out
%   penalty  : logical m-vector with jth element true if factor should be
%                down in the order
%   orderalg :  1 - hybrid algorithm that eliminates subset factors
%               2 - optimal order
%               3 - greedy algorithm 
% OUTPUT
%   order    : (m-1)x2 matrix with the order that pairs of factors are
%                 combined. The combined factor is placed back in the
%                 position of the first of the pair of factors

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

function [order,cost]=getorder(V,n,sumvar,penalty,orderalg)
  m=size(V,1);
  switch orderalg
    case 0
      if  m>7 
        if m>20
          [order,cost]=minsize(V,n,sumvar,penalty);  
        else
          [o1,v1,remaining,cost1]=initialorder(V,n,sumvar,penalty);
          [order,cost]=optfactororder(v1,n,sumvar,penalty(remaining));
          order=remaining(order);
          order=[o1;order];
          cost=cost+cost1;
        end
      else
        [order,cost]=optfactororder(V,n,sumvar,penalty);
      end
    case 1
      [order,cost]=minsize(V,n,sumvar,penalty); 
    case 2
      [order,cost]=optfactororder(V,n,sumvar,penalty); 
    case 3
      [order,cost]=greedy(V,n,sumvar,penalty); 
  end
  % move factor 1 down in the order to the last possible position
  % all subsequent joins will have factor 1 first
  if penalty(1)
    zz=any(order==1,2);
    [zz,ind]=sort(zz);
    order=order(ind,:);
  end
  
  

% gets an initial ordering by eliminating any factors that do not increase the
% size of the factor being combined with
% If penalty(i)=1 the ith factor is not altered
function [order,V,remaining,cost]=initialorder(V,n,sumvar,penalty)
m=size(V,1);
remaining=true(1,m)&~penalty;
cost=0;
order=zeros(0,2);
i=1;
while i<m
  if remaining(i)
    vi=V(i,:);
    if all(vi==0)
      error('shouldn''t happen')
    end
    mincost=inf; 
    jj=[];
    rv=find(remaining);
    for j=rv(rv>i)
      vj=V(j,:);
      if all(vj(vi)) || all(vi(vj)) % all elements in i are in j or vice versa
        costij=prod(n(vi | vj));
        if costij<mincost
          mincost=costij;
          jj=j;
        end
      end
    end
    if ~isempty(jj)
      order=[order;i jj]; %#ok<*AGROW>
      vi=vi | V(jj,:);
      cost=cost+mincost;
      V(jj,:)=false;
      remaining(jj)=false;
      V(i,:)=vi;
      V(i,vi & sumvar & sum(V,1)==1)=false; % sum singleton dimensions
      i=1;
    else
      i=i+1;
    end
  else
    i=i+1;
  end
end
remaining=remaining | penalty;
V=V(remaining,:);
remaining=find(remaining);


%% 

function [order,cost]=optfactororder(V,n,sumvar,penalty)
m=size(V,1);
V=[V;false(m-2,size(V,2))];
remaining=[true(1,m) false(1,m-2)];
order=[];
vi=1;
[order,cost]=recursivecall(vi,remaining,order,[0 0],[inf inf],[],V,n,sumvar,penalty);
% renumber the factors so new factors replace the first input factor
for i=1:m-1
  if order(i,2)>m, order(i,2)=order(order(i,2)-m,1); end
  if order(i,1)>m, order(i,1)=order(order(i,1)-m,1); end
end
for i=1:m-1
  if order(i,2)<order(i,1)
    temp=order(i,1);
    order(i,1)=order(i,2);
    order(i,2)=temp;
    oi=order(i+1:end,:); oi(oi==temp)=order(i,1);
    order(i+1:end,:)=oi;
  end
end
cost=sum(cost);

function [minorder,mincost]=recursivecall(vi,remaining,order,cost,mincost,minorder,V,n,sumvar,penalty)
rv=find(remaining);
if length(rv)==2
  minorder=[order;rv];
  if any(penalty(rv)>0)
    mincost=[cost(1) cost(2)+prod(n(any(V(rv,:),1)))];
  else
    mincost=[cost(1)+prod(n(any(V(rv,:),1))) cost(2)];
  end
  return
end
nextin=rv(end)+1;
for i=rv(1:end-1)
  nextvi=max(vi,i);
  if i<vi
    jj=rv(end);
  else
    jj=rv(rv>i);
  end
  for j=jj
    vij=V(i,:)|V(j,:);
    if cost(2)>0 || penalty(i)>0 || penalty(j)>0
      costij=[cost(1) cost(2)+prod(n(vij))];
    else
      costij=[cost(1)+prod(n(vij)) 0];
    end
    if any(costij<mincost)
      nextr=remaining;
      nextr([i j])=false; nextr(nextin)=true;
      orderij=[order;i j];
      Vij=V;
      Vij([i j],:)=false;
      vij(vij & sumvar & sum(Vij,1)==0)=false;
      Vij(nextin,:)=vij;
      if costij(2)>0, penalty(nextin)=true; else penalty(nextin)=false; end
      [orderij,costij]=recursivecall(nextvi,nextr,orderij,costij,mincost,minorder,Vij,n,sumvar,penalty);
      if costij(2)<=mincost(2)
        if costij(2)<mincost(2) || costij(1)<mincost(1)
          mincost=costij; minorder=orderij;
        end
      end
    end
  end
end

%%


function [order,cost]=minsize(V,n,sumvar,penalty)
m=size(V,1);
pencost=zeros(1,m); for i=1:m, pencost(i)=prod(n(V(i,:))); end
pencost=max(pencost);
cost=0;
remaining=1:m;
order=zeros(m-1,2);
for k=1:m-1
  ck=inf;
  increasek=inf;
  rv=find(remaining);
  for i=rv(1:end-1)
    vi=V(i,:);
    for j=rv(rv>i)
      vij=vi|V(j,:);
      increaseij=sum(vij)-max(sum(vi),sum(V(j,:)));
      if increaseij<=increasek
        cij=prod(n(vij))+(penalty(i) | penalty(j))*pencost;
        if (increaseij<increasek) || (increaseij==increasek && cij<ck)
          increasek=increaseij;
          ck=cij;
          ii=i; jj=j;
          vk=vij;
        end
      end
    end
  end
  V(ii,:)=vk;
  V(jj,:)=false;
  vk(vk & sumvar & sum(V,1)==1)=false;
  V(ii,:)=vk;    
  order(k,:)=[ii jj];
  remaining(jj)=false;
  cost=cost+ck;
end


function [order,cost]=greedy(V,n,sumvar,penalty)
[m,d]=size(V);
objord=[1 3 2];
do=length(objord);
remaining=1:m;
order=zeros(m-1,2);
cost=0;
for k=1:m-1
  costk=inf+zeros(1,do);
  for i=1:length(remaining)-1
    ifact=remaining(i);
    for j=i+1:length(remaining);
      jfact=remaining(j);
      costs=joincosts(ifact,jfact,V,n,sumvar);
      for ii=1:do
        if costs(objord(ii))<costk(ii)
          costk=costs(objord);
          order(k,:)=[ifact jfact];
          costkeep=costs;
          break;
        elseif costk(ii)==costs(objord(ii))
        else
          break
        end
      end
    end
  end
  cost=cost+prod(abs(costkeep(1:3)));
  V(order(k,1),:)=(V(order(k,2),:) | V(order(k,2),:));
  V(order(k,2),:)=false;
  V(i,V(i,:) & sumvar & sum(V,1)==1)=false; % eliminate summed variables
  remaining(remaining==order(k,2))=[];
end
  
 
function [costs]=joincosts(i,j,V,n,sumvar)
vi=V(i,:);
vj=V(j,:);
V(i,:) =vi | vj;
V(j,:)=false;
ind=V(i,:) & sumvar & sum(V,1)==1;
vm=vi & vj;  % matched or summed
vu=(vi | vj) & ~vm;
vs=vm & ind;
vm=vm & ~ind;
costs=[prod(n(vu)) prod(n(vm)) -prod(n(vs)) ...
        sum(vu) sum(vm) -sum(vs)];
      
      
      
%% legacy code
      

% finds the variable with the least number of 
function [order,cost]=sumproductorder(V,n,sumvar,penalty)
[m,d]=size(V);
cost=[0 0];
order=zeros(m-1,2);
k=1;
sumvar=sumvar & any(V,1);
while true
  vi=zeros(1,d);
  for i=1:d
    if sumvar(i)
      vi(i)=sum(any(V(V(:,i),:),1));
    end
  end
  [dd,ii]=min(vi);
  ii=find(V(:,ii));
  for i=2:length(ii)
    order(k,1)=ii(1);
    order(k,2)=ii(i);
    ind=V(ii(1),:) | V(ii(2),:);
    c=prod(n(ind));
    if ii(1)==1
      cost(1)=cost(1)+c;
    else
      cost(2)=cost(2)+c;
    end
    k=k+1;
  end
  V(ii(1),:)=any(V(ii,:),1);
  V(ii(2:end),:)=false;
  ind=V(ii(1),:) & sumvar & sum(V,1)==1;
  V(ii(1),ind)=false;
  sumvar(ind)=false;
  if k>=m, break; end
end

      