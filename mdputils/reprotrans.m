%reprotrans transition probability matrix for reproduction stage
% USAGE
%   [P,p]=reprotrans(N,f);
% INPUTS
%   N  : population size (scalar positive integer)
%   f  : M+1 probability vector, f(i+1)=Pr(i offspring)
%          or an m-element cell array composed of probability vectors
%
% OUTPUTS
%   P  : transition matrix when each member of a population of N or less
%           is assigned to one of m age/stage classes
%   p  : transition matrix for a single individual
%
% Offspring are added to the lowest category
% Implicit category m+1 is defined as absent/dead (makes total # equal N)
% If the number of offspring causes the total population to exceed N then
%   category 1 is truncated
% To get the number of individuals in each row/column of P use
%   S=simplexgrid(m+1,N,N,0);
% This returns an m-column matrix with rows summing to N or less
%
% Examples: 
% single class: reprotrans(3,[0.4;0.6]);
%   1     0     0     0  
%   0  0.40     0     0  
%   0  0.60  0.16     0  
%   0     0  0.84     1
% juvenile/adult classes: reprotrans(3,{1,[0.4;0.6]})
%   1     0     0     0     0     0     0     0     0     0  
%   0  0.40     0     0     0     0     0     0     0     0  
%   0     0  0.16     0     0     0     0     0     0     0  
%   0     0     0     1     0     0     0     0     0     0  
%   0     0     0     0     1     0     0     0     0     0  
%   0  0.60     0     0     0  0.40     0     0     0     0  
%   0     0  0.84     0     0     0     1     0     0     0  
%   0     0     0     0     0     0     0     1     0     0  
%   0     0     0     0     0  0.60     0     0     1     0  
%   0     0     0     0     0     0     0     0     0     1  
function [P,p,theta,Q]=reprotrans(N,f)
if isnumeric(f)
    f={f};
end
try
  m=length(f);  % # of stages
  n=cellfun(@length,f);
  p=[blkdiag(f{:});zeros(1,m)];
  p=[p [zeros(sum(n),1);1]];
  P=catcountP(N,size(p,1),size(p,2),p);
  theta=zeros(0,m);
  for i=1:m
    temp=zeros(n(i),m);
    if i==1
      temp(:,1)=(1:n(i))';
    else
      temp(:,1)=(0:n(i)-1)';
      temp(:,i)=ones(n(i),1);
    end
    theta=[theta;temp]; %#ok<AGROW>
  end
  S0=simplexgrid(size(theta,1)+1,N,N,0);
  St=S0*theta;
  % truncate the first stage if population has maxed out
  St(:,1)=min(N-sum(St(:,2:end),2),St(:,1));

  S=simplexgrid(m+1,N,N,0);
  [~,Q]=matchrows(St,S);
  P=Q*P;
catch
  [P,p]=reprotransSmallMemory(N,f);
end
return


function [P,p]=reprotransSmallMemory(N,f)
[M,m]=size(f);  % max offspring +1, # of stages
p=[f zeros(M,1);zeros(1,m) 1];  % add in the residual category
P=catcountP(N,size(p,1),size(p,2),p);
S0=simplexgrid(m+1,N,N,0);
S1=simplexgrid(M+1,N,N,0);
ns0=size(S0,1);
[ii,jj,pp]=find(P);
x=S1(ii,:)*(0:M-1)' + S0(jj,1);
x=min(N-sum(S0(jj,2:end),2),x);
ii=matchrows([x S0(jj,2:end)],S0);
P=sparse(ii,jj,pp,ns0,ns0);
