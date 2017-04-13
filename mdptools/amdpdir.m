% amdpdir Adaptive Management using Dirichlet priors
% USAGE
%   [alpha,P,S,X,Ix]=amdpdir(p,inc,M,S,X,Ix);
% INPUTS
%   p       : ns x nx column stochastic transition matrix
%               with uncertain values set to NaN
%   inc     : positive integer, define the grid increment to be 1/inc
%   M       : maximal value for sum(alpha_i)
%   S       : state matrix for no model uncertainty case
%   X       : state/action matrix for no model uncertainty case
%   Ix      : state index for no model uncertainty case
% OUTPUTS
%   alpha   : matrix of alpha grid values
%   P       : transition matrix for the augmented problem
%   S       : state matrix for the augmented problem
%   X       : state/action matrix for the augmented problem
%   Ix      : state index for the augmented problem
%
% Set up an adaptive management problem with one or more columns of p containing 
% uncertain values that are represented by a Dirichlet distribution. Each of the
% uncertain parameter values are represented by a value of alpha. Use
%   Ep=amdpdirp(p,alpha);
% to obtain the expected value of the p matrix for a specified alpha vector
function [alpha,P,S,X,Ix]=amdpdir(p,inc,M,S,X,Ix)
[ns,nx]=size(p);
inc=inc+zeros(1,nx);
M=M+zeros(1,nx);
iM=inc.*M;
na=ones(1,nx);     % # of alpha vectors for column j
uind=zeros(1,nx);  % # of uncertain probability values in column j
alpha=cell(1,nx);  % alpha values of each column   
pb=cell(1,nx);     % associated probability values
for j=1:nx
  uind(j)=sum(isnan(p(:,j)));
  if uind(j)>1
    iM(j)=iM(j)-uind(j);
    pp=1-sum(p(~isnan(p(:,j)),j));  % uncertain values must sum to pp
    alpha{j}=simplexgrid(uind(j)+1,iM(j),iM(j),0,1);
    pb{j}=vxm(pp./(double(sum(alpha{j},2))+uind(j)),double(alpha{j}+1));
    na(j)=size(alpha{j},1);
  elseif uind(j)==1
    error('Only 1 unknown parameter in a column of p is not permitted')
  end
end

Na=prod(na);
N=ns*Na;
% define initial allocation of space for P
P=sparse([],[],[],N,0,(ns*nx-sum(p(:)==0))*Na);
inds=cell(1,nx); for j=1:nx, inds{j}=(1:na(j))'; end 
% loop over the state/action combinations
for j=1:nx
  l=0;
  Pi=sparse(N,Na);
  % loop over the future states
  for i=1:ns
    if isnan(p(i,j))
      l=l+1;
      pij=pb{j}(:,l);
      bj=alpha{j}; 
      ii=sum(alpha{j},2)<=(iM(j)-inc(j));
      bj(ii,l)=bj(ii,l)+inc(j); 
      inds{j}=simplexindex(bj,uind(j)+1,iM(j),iM(j));  
      ind=index(na,rectgrid(inds));
      Pi=Pi+sparse(ind+(i-1)*Na,1:Na,pij,N,Na); 
    else
      pij=p(i,j);
      if pij~=0
        Pi=Pi+sparse((1:Na)'+(i-1)*Na,1:Na,pij,N,Na); 
      end
    end
  end
  P=[P Pi]; %#ok<AGROW>
  inds{j}=(1:na(j))';
  clear ind ii pij Pi
end
clear inds Pi bj pb

for j=1:nx
  alpha{j}=(alpha{j}+1)/inc(j);
end  
alpha=rectgrid(alpha);
% expand the state and state/action matrices to include the new belief states
if nargin>3 && nargout>2
  S=rectgrid(S,alpha);
  if nargin>4 && nargout>3
    X=rectgrid(X,alpha);
    if nargin>5 && nargout>4
      nb=size(alpha,1);
      Ix=ones(nb,1)*(nb*(Ix(:)'-1)) + (1:nb)'*ones(1,numel(Ix));
      Ix=Ix(:);
    end
  end
end

function ind=index(ns,x)
d=size(x,2);
ind=x(:,d);
nn=ns(d);
for i=d-1:-1:1
  ind=ind+(x(:,i)-1)*nn;
  nn=nn*ns(i);
end