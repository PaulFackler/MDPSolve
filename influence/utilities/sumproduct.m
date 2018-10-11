% sumproduct Applies the sumproduct operator to a set of factors
% USAGE
%   [f,ordinfo,cost,Iexpand]=sumproduct(F,Fvars,n,rowlist,collist,options);
% INPUTS
%   F       : a cell array of m factors
%   Fvars   : a 1xm cell array with indices of the variables in each factor
%   n       : a d-vector with the sizes of each variable
%   rowlist : list of variables associated with rows of the output
%   collist : list of variables associated with columns of the output
%   options : structure variable with control options (described below)
% OUTPUTS
%   f       : the final factor, a prod(n(rowlist)) x prod(n(collist)) matrix
%   ordinfo : (m-1)x7 cell array containing information on the order
%                factors are processed and the variables involved
%                Use with orddisp function.
%   cost    : total processing cost (number of multiply operations)
%   Iexpand : if f has less columns than prod(n(collist)) then
%               f should be calculated as f=f(:,Iexpand) 
%             otherwise Iexpand is empty
%             If Iexpand is not requested it is used prior to returning f
%
% Options:
%   order       :  the order of processing of the variables
%   orderalg    : algorithm to determine elimination order
%                     0 for default
%                     1 forces greedy
%                     2 forces optimal
%   orderonly   : 1 to return order info and skip processing (f and Iexpand set to [])
%   forcefull   : 0 use sparse factors
%                 1 convert sparse factors to full
%   spthreshold : number on [0,1] converts final factor to sparse if sparsity
%                   ratio is less than threshold
%   feasible    : logical vector set to 1 for elements of collist that
%                   should be retained
%   print       : print level, 0: none, 1: moderate, 2: heavy
%
% F{i} is an array with prod(n(Fvars{i}) elements. Implicitly it is a d_i
%   dimensional array (where d_i is the length of Fvars{i}) with n(Fvars{i}) 
%   values in dimension i. In practice F{i} can have any number of dimensions 
%   between 1 and d_i. Note that elements of F{i} are arranged in reverse
%   lexicographic order. The result of this operation is always returned as
%   a matrix with output indices ordered according to rowlist and collist.
%   Any indices included in the Fvars vectors and not included rowlist or
%   collist are summed out of the final result.

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

function [f,ordinfo,cost,Iexpand]=sumproduct(F,Fvars,n,rowlist,collist,options)
  order           = [];  % order for processing
  orderalg        = 0;   % algorithm to determine elimination order
                         % 0 for default
                         % 1 forces greedy
                         % 2 forces optimal
  orderonly       = 0;   % 1 to return order info (f and Iexpand set to [])
  forcefull       = 0;   % 1 to force sumproduct to use full factors
spthreshold       = 0.1; % number on [0,1] converts to sparse if sparsity
                         %   ratio is less than threshold
feasible          =[];   % vector of indices to contract the output
  print           = 0;   % print level, 0: none, 1: moderate, 2: heavy
  if nargin>=6 && ~isempty(options)
    if isfield(options,'order'),       order=options.order;             end
    if isfield(options,'orderalg'),    orderalg=options.orderalg;       end
    if isfield(options,'orderonly'),   orderonly=options.orderonly;     end
    if isfield(options,'forcefull'),   forcefull=options.forcefull;     end
    if isfield(options,'spthreshold'), spthreshold=options.spthreshold; end
    if isfield(options,'feasible'),    feasible=options.feasible;       end
    if isfield(options,'print'),       print=options.print;             end
  end
d=length(n);
m=length(F);
V=false(m,d);
for i=1:m, V(i,Fvars{i})=true; end
outorder=[rowlist collist];
sumvar=true(1,d); sumvar([rowlist collist])=false;

makefunc=false;
Iexpand=[];
reorder=[];
cost=0;
if m==1
  i1=1;
  cost=0;
  ordinfo=[];
  if orderonly; 
    f=[];
    return; 
  end
else
  penalty=false(1,m); for i=1:m, if isempty(F{i}), penalty(i)=true; end; end
  % may one day allow user supplied processing order
  if print>0
    disp('determining processing order')
    start=cputime;
  end
  % get processing order
  if ~isempty(order)
    warning('user supplied orders are not yet implemented; getorder is used')
    [order,cost]=getorder(V,n,sumvar,penalty,orderalg);
  else
    [order,cost]=getorder(V,n,sumvar,penalty,orderalg);
  end
  if print>0
    fprintf('time taken: %8.3f\n',cputime-start)
    fprintf('order cost: %16i\n',cost)
  end
  ordinfo=getvarorder(V,order,sumvar);
  if orderonly
    f=[];
    return
  end
  
  if print>0
    if print>1
      disp('factor process order:')
      disp(order)
      disp('processing factors (with factor size): ')
    else
      disp('processing factors')
    end
    start=cputime;
  end
  
  % factor processing loop
  for i=1:size(order,1)
    i1=order(i,1); 
    i2=order(i,2);
    x2z=Fvars{i1}; 
    nx=n(x2z);
    y2z=Fvars{i2}; 
    ny=n(y2z);
    if i<m-1
      outvars=[ordinfo{i,5:7}];
      nzout=[prod(n([ordinfo{i,5}])) prod(n([ordinfo{i,6}]))];
      Fvars{i1}=[ordinfo{i,5:6}];
    else
      temp=collist(ismember(collist,[ordinfo{i,5:6}]));
      outvars=[rowlist temp ordinfo{i,7}];
      nzout=[prod(n(rowlist)) prod(n(temp))];
      Fvar1in=Fvars{i1};
      Fvars{i1}=[rowlist temp];
    end
    for j=1:length(x2z)
      if any(abs(x2z(j))==ordinfo{i,7})
        x2z(j) = -find(abs(x2z(j))==outvars);
      else
        x2z(j) = find(x2z(j)==outvars);
      end
    end
    for j=1:length(y2z)
      if any(abs(y2z(j))==ordinfo{i,7})
        y2z(j) = -find(abs(y2z(j))==outvars);
      else
        y2z(j) =  find(y2z(j)==outvars);
      end
    end  
    if print>1
      fprintf('%3i  %3i  %15i\n',[i1 i2 prod(n(outvars))])
    end
    cost=cost+prod(n(outvars));
    if isempty(F{i1})
      if ~makefunc
        makefunc=true;
        startfunc=i;
        control=cell(m-i,5);
      end
    end
    if makefunc
      % permutes the second factor so it need not be latter
      % not needed if there are only two factors as it will be done latter
      if isempty(F{i1})  && size(control,1)>1    %%%%%%%%%%%%%%%%  CHECK
        [q,ix,iy]=intersect(x2z,y2z); 
        nsum=sum(q<0);
        nmatch=length(q)-nsum;
        if nsum>0 && nmatch>0, 
           q=[q(nsum+1:end) q(1:nsum)];
           iy=[iy(nsum+1:end)' iy(1:nsum)'];
        end
        iy=[find(~ismember(y2z,q)) iy(:)'];
        nyy=[prod(ny(~ismember(y2z,q))) prod(ny(ismember(y2z,q)))];
        if forcefull==2, F{i2}=full(F{i2}); end
        F{i2}=sppermute(F{i2},iy,ny,nyy); 
        y2z=y2z(iy);
        ny=ny(iy);
        Fvars{i2}=Fvars{i2}(iy);
      end
      control(i-startfunc+1,:)={i1,[x2z;nx],i2,[y2z;ny],nzout};
    else
      F{i1}=tprodm(F{i1},[x2z;nx],F{i2},[y2z;ny],nzout);
      F{i2}=[];
    end
  end
  if print>0
    fprintf('time taken to process factors: %8.3f\n',cputime-start)
  end
end

% If all non-summed variables are not present in the final factor
% create in Iexpend index to expand the factor to the correct size.
% This operates on the columns of the factor, i.e., F=F(:,Iexpand)
% (this assumes that the rows are all represented - if not an error occurs).
Fvars1=Fvars{i1};
if length(Fvars1)~=length(outorder) || any(Fvars1~=outorder)
   if m>1 && length(Fvars1)==length(outorder)
    warning('shouldn''t happen')
   end
   [F{i1},Fvars1,nout1,Iexpand]=adjustfactor(F{i1},n,Fvars1,rowlist,collist);
end

% if feasibility constraints are present incorporate them into Iexpand
if ~isempty(feasible)
  if ~isempty(Iexpand)
    Iexpand=Iexpand(feasible);
  else
    Iexpand=feasible;
    if islogical(Iexpand), Iexpand=find(Iexpand); end
  end
end

if isempty(F{i1})  % an EV function requested
  if size(control,1)>1 
    f=getfunc(F,control,reorder,Iexpand);
  else   % only one remaining factor so form transition matrix
    F=F{i2};
    [F,Fvars1,nout1,Iexpand]=adjustfactor(F,n,Fvars{i2},Fvar1in,collist);
    if nnz(F)/numel(F)<spthreshold
       F=sparse(F);
    end
    if ~isempty(feasible)
      if ~isempty(Iexpand)
        Iexpand=Iexpand(feasible);
      else
        Iexpand=feasible;
        if islogical(Iexpand), Iexpand=find(Iexpand); end
      end
    end
    f=getfunctranmat(F,Iexpand);
  end
else     % a transition matrix requested
  F=F{i1};
  if nnz(F)/numel(F)<spthreshold
    F=sparse(F);
  end
  f=F;
  if ~isempty(Iexpand) && nargout<4
    f=f(:,Iexpand); 
  end
end

function [F,Fvars,nout,Iexpand]=adjustfactor(F,n,varorder,rlist,clist)
  if ~isempty(F) && numel(F)~=prod(n(varorder))
    error('size of F does not match sizes of varorder')
  end
  dF=length(varorder);
  outorder=[rlist clist];
  [~,iin,iout]=intersect(varorder,outorder);
  if length(iin)~=dF
    error('F contains dimensions that are not in output lists')
  end
  reorder=zeros(1,length(outorder));
  reorder(iout)=iin;
  reorder(reorder==0)=[];
  if length(outorder)==dF
    nout=[prod(n(rlist)) prod(n(clist))];
    if ~isempty(F)
      F=sppermute(F,reorder,n(varorder),nout);
    end
    Fvars=varorder(reorder);
    Iexpand=[];
  else
    if all(ismember(rlist,varorder))
      nout=prod(n(rlist));
      nout=[nout prod(n(varorder))/nout];
      if ~isempty(F)
        F=sppermute(F,reorder,n(varorder),nout);
      end
      Fvars=varorder(reorder);
      ic=ismember(clist,Fvars);
      newsize=ones(1,length(clist));
      newsize(ic)=n(Fvars(length(rlist)+1:end));
      Iexpand=reshape((1:nout(2))',newsize);
      expand=n(clist);
      expand(ic)=1;
      Iexpand=repmat(Iexpand,expand);
      Iexpand=Iexpand(:);
    else
      error('not implemented to handle missing variables in rlist')
    end
  end
    
% Creates a 7 column cell array with columns 1 & 2 containing
% the two factors (x and y) in the join operation and, 
% for each join operation, columns 3-7 indicate which variables are in:
%   3)x 4)y 5)unmatched 6)matched 7)matched & summed
function ordinfo=getvarorder(V,order,sumvar)
  m=size(V,1);
  ordinfo=cell(m-1,7);
  for k=1:m-1
    o1=order(k,1); 
    o2=order(k,2);
    v1=V(o1,:); 
    v2=V(o2,:); 
    ordinfo{k,1}=o1;
    ordinfo{k,2}=o2;
    ordinfo{k,3}=find(v1);
    ordinfo{k,4}=find(v2);
    ordinfo{k,5}=[find(v1 & ~(v1&v2)) find(v2 & ~(v1&v2))];
    ordinfo{k,6}=v1 & v2;
    v1=v1 | v2;
    V(o1,:)=v1;
    V(o2,:)=false;
    ss=v1 & sumvar & sum(V,1)==1;
    ordinfo{k,6}=find(ordinfo{k,6} & ~ss);
    ordinfo{k,7}=find(ss);
    v1(ss)=false;
    V(o1,:)=v1;
  end
  
  
% this is included to use a transition matrix in the EV function 
% if only one factor remains at the end of processing.
% not yet implemented
function func=getfunctranmat(F,Iexpand)
  func=@(x)  evalfunctranmat(x,F,Iexpand);

function f=evalfunctranmat(x,F,Iexpand)
  f=F'*x(:);
  if ~isempty(Iexpand);
    f=f(Iexpand);
  end


% gets a function of F1 when F1 is passed as empty
function func=getfunc(F,control,reorder,Iexpand)
  func=@(x)  evalfunc(x,F,control,reorder,Iexpand)';

% the evaluation routine called to evaluate the function
function f=evalfunc(x,F,control,reorder,Iexpand)
  F{control{1,1}}=x;   % set the input factor as the first factor
  % run through the processing order 
  for i=1:size(control,1)
    i1=control{i,1};
    F{i1}=tprodm(F{i1},control{i,2},F{control{i,3}},control{i,4},control{i,5},0);
  end
  f=F{i1};
  % reorder the output if needed
  if ~isempty(reorder)
    error('shouldn''t happen')
    if issparse(f)
      f=sppermute(f,reorder,n(Fvars),size(f));
    else
      f=reshape(permute(reshape(f,n(Fvars)),reorder),size(f));
    end
  end
  % expand the output if needed
  if ~isempty(Iexpand), f=f(Iexpand); end
  


  
