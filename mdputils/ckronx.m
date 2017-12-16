% CKRONX The product of repeated Kronecker products and a matrix. 
% USAGE      
%   [C,order]=ckronx(A,B,options);
% INPUTS
%   A       : a d-element cell array with element i an m(i) x n(i) matrix
%   B       : a compatible matrix prod(n) x p (prod(m) x p if transpose=1)
%   options : an options structure (see below)
%
% Fields in options structure:
%   ind       : a selection vector of numbers on 1..d that reorders 
%                 the elements of A [default: 1:d]
%   transpose : d-vectors of logicals: 1 to use the transpose of A(i) 
%                 a scalar entry will be expanded to a d-vector [default: 0]
%   reorder   : use optimal ordering of operations
%   forward   : 1 forces use of the forward algorithm, 0 forces the backward
%   print     : print information about the operation
% OUTPUTS  
%   C         :  prod(m) x p matrix (prod(n) x p if transpose=1)
%   order     :  1 x d vector with the optimal order of variables
% Solves (A1 x A2 x...x Ad)*B
% where x denotes Kronecker (tensor) product.
% The Ai are passed as a cell array A. 
% A must be a vector cell array containing 2-D numerical arrays (matrices).
% If A is a matrix the function returns A*B (or A'*B if transpose=1)
%
% Alternative input syntax (for backward compatability):
%   C=ckronx(A,B,ind,transpose);

% Adapted from the ckronx function in the CompEcon Toolbox
%(c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% Modified in 2017 (c) by Paul L. Fackler

function [C,order]=ckronx(A,B,options,transpose)

if nargin<2, error('At least two parameters must be passed'), end
if ~exist('transpose','var'), transpose=false; end
reorder=0;
useforward=[];
print=0;
order=[];

% handle case if a single numeric matix is passed as A
if ~iscell(A)                         % A is a matrix: return A*B
  if exist('options','var') && ~isempty(options) && isfield(options,'transpose')
    transpose = options.transpose;
  end
  if transpose
    if size(A,1)~=size(B,1)
      error('A and B are not conformable')
    end
    C=A'*B;
  else
    if size(A,2)~=size(B,1)
      error('A and B are not conformable')
    end
    C=A*B;
  end
  return
end

% handle cell array case
d=numel(A);
% get default order
if ~exist('options','var')  || isempty(options)
  ind=1:d;
else
  if isstruct(options)
    ind=1:d;       
    if isfield(options,'ind'),       ind=options.ind;              end  
    if isfield(options,'transpose'), transpose=options.transpose;  end   
    if isfield(options,'reorder'),   reorder=options.reorder;      end   
    if isfield(options,'forward'),   useforward=options.forward;   end   
    if isfield(options,'print'),     print=options.print;          end           
  else
    ind=options;
  end
end
if length(transpose)==1, transpose=repmat(transpose,1,d); end
A=A(ind);
transpose=transpose(ind);
m=zeros(1,d);  % # of rows (cols if transpose(i)=1)
n=zeros(1,d);  % # of cols (rows if transpose(i)=1)
q=zeros(1,d);  % # of non-zeros
for i=1:d, 
  if transpose(i), [n(i),m(i)]=size(A{i}); 
  else             [m(i),n(i)]=size(A{i}); 
  end
  if issparse(A{i}), q(i)=nnz(A{i}); 
  else               q(i)=numel(A{i}); 
  end
end
if prod(n)~=size(B,1)
  error('A and B are not conformable')
end
if d==1, 
  if transpose(1), C=A{1}'*B; 
  else             C=A{1}*B; 
  end
  return; 
end

cost=@(m,n,q)sum(q.*cumprod([1 m(1:end-1)]).*fliplr(cumprod([1 fliplr(n(2:end))])));
fcost=cost(m,n,q);
bcost=cost(n,m,q);

if reorder
  [~,order]=sort((m-n)./q);
  gcost=cost(m(order),n(order),q(order));
  if gcost+prod(m)+prod(n) < min(fcost,bcost) 
    p=size(B,2);
    C=permute(reshape(B,[fliplr(n) p]),[d+1 d+1-fliplr(order)]);
    n=n(order);
    transpose=transpose(order);
    for i=1:d
      C=reshape(C,numel(C)/n(i),n(i));
      if transpose(i), C=A{order(i)}'*C';
      else             C=A{order(i)}*C';
      end
    end
    C=ipermute(reshape(C,[fliplr(m(order)) p]),[d+1-fliplr(order) d+1]);
    C=reshape(C,numel(C)/size(B,2),size(B,2));
    if print
      fprintf('order: '); fprintf('%1.0f ',order); fprintf('\n');
      disp('reordered cost')
      fprintf('%25.0f\n',gcost);
      disp('reorder costs')
      fprintf('%25.0f\n',[prod(n)*p; prod(m)*p]);
    end
    return
  end
end

% if reordering is not done
if print
  disp('forward & backward costs')
  fprintf('%25.0f\n',[fcost;bcost])
  disp('full Kronecker cost')
  fprintf('%25.0f\n',prod(q));
  if ~reorder
    [~,order]=sort((m-n)./q);
    gcost=cost(m(order),n(order),q(order));
    fprintf('order: '); fprintf('%1.0f ',order); fprintf('\n');
    disp('reordered cost')
    fprintf('%25.0f\n',gcost);
    disp('reorder costs')
    p=size(B,2);
    fprintf('%25.0f\n',[prod(n)*p; prod(m)*p]);
  end
end
if isempty(useforward)
  if fcost<bcost || (fcost==bcost && all(~transpose)) % use forward approach
    C=forward(A,B,d,n,transpose);
  else  % use backward approach
    C=backward(A,B,d,n,transpose);
  end
else
  if useforward~=0 % use forward approach
    C=forward(A,B,d,n,transpose);
  else  % use backward approach
    C=backward(A,B,d,n,transpose);
  end
end

function C=forward(A,B,d,n,transpose)
  C=B';
  for i=1:d
    C=reshape(C,numel(C)/n(i),n(i));
    if transpose(i), C=A{i}'*C';
    else             C=A{i}*C';
    end
  end
  C=reshape(C,numel(C)/size(B,2),size(B,2));
    
    
function C=backward(A,B,d,n,transpose)
  C=B;
  for i=d:-1:1
    C=reshape(C,n(i),numel(C)/n(i));
    if transpose(i), C=C'*A{i};
    else             C=C'*A{i}';
    end
  end
  C=reshape(C,size(B,2),numel(C)/size(B,2))';