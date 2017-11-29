% EVindices Creates index vectors for EV functions (called by EVcreate)
% USAGE
%   [Ip,Iy,Jp,Jy,Xp,yn]=EVindices(X,parents,pn);
% INPUTS
%   X       : matrix of state/action combinations
%   parents : m-element cell array of conditioning variables (parents)
%   pn      : m-vector of column sizes of the CPTs (optional to check compatibility)
% OUTPUTS
%   Ip,Iy,Jp,Jy : index vectors to guide EV function operations
%   Xp          : matrix of conditioning variables for each CPT
%   yn          : m-vector of # of columns in the output at each evaluation step 
%
% parents{i} contains the columns of X that condition the ith variable
% Thus Xpi=unique(X(:,parents{i}),'rows') contains the unique elements of the
%   conditioning variables that condition output variable i.
% It must be the case that the number of rows in Xi matches the number 
%   of columns in p{i}: size(Xpi,1)=np(i)
  
function [Ip,Iy,Jp,Jy,Xp,yn]=EVindices(X,parents,pn)
  ind64=1;
  m=length(parents);
  parentscombined=cell(1,m);
  parentscombined{1}=parents{1};
  for i=2:m
    parentscombined{i}=union(parentscombined{i-1},parents{i});
  end
  Ip=cell(1,m); Iy=cell(1,m); Jp=cell(1,m); Jy=cell(1,m);
  Xp=cell(1,m);
  yn=zeros(1,m);
  Xi=X;
  J=(1:size(X,1))';
  for i=m:-1:2
    yn(i)=size(Xi,1);
    if isequal(parents{i},parentscombined{i})
      Ip{i}=[]; 
      Jp{i}=[];
      Xp{i}=Xi;
    else
      ind=matchindices(parents{i},parentscombined{i});
      [Ip{i},Xp{i}]=getI(Xi,ind);
      Jp{i}=Ip{i}(J);
    end
    
    if nargin>=3 && size(Xp{i},1)~=pn(i)
      error(['parents{' num2str(i) '} is incompatible with p{' num2str(i) '}'])
    end
    if isequal(parentscombined{i-1},parentscombined{i})
      Iy{i}=[]; 
      Jy{i}=[];
    else
      ind=matchindices(parentscombined{i-1},parentscombined{i});
      [Iy{i},Xi]=getI(Xi,ind);
      Jy{i}=Iy{i}(J);
      J=Jy{i};
    end
  end
  ind=matchindices(parents{1},parentscombined{1});
  [Ip{1},Xp{1}]=getI(Xi,ind);
  Jp{1}=Ip{1}(J);
  yn(1)=size(Xi,1);
  if ind64
    for i=1:m
      Ip{i}=uint64(Ip{i});
      Iy{i}=uint64(Iy{i});
      Jp{i}=uint64(Jp{i});
      Jy{i}=uint64(Jy{i});
    end
  end
    
% matchindices Finds the values of ind1 that match those of ind0
function ind=matchindices(ind0,ind1)
  [~,ind]=ismember(ind0,ind1);
  ind=ind(ind>0);