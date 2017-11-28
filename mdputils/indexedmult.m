% indexedmult An indexed multiplcation of two arrays
%   z=indexedmult(x,xind,y,yind,checks);
% INPUTS
%   x      : m x n matrix
%   xind   : length w index vector of values on 1,...,n
%   y      : p x q matrix (p<=m and m/p integer)
%   yind   : length w index vector of values on 1,...,q
%   checks : 0/1 check the correctness of the inputs
%              set to 0 to avoid time consuming checks if inputs
%              are known to be correct
% OUTPUT
%   z      : m/p x w matrix
%
% First reshape x to be m/p x p x n
% Then z(i,j) = x(i,:,xind(j)) * y(:,yind(j))
%
% Note that size checks are made on inputs (these checks are fast)
% but user must ensure that 1<=xind(i)<=size(x,2) and 1<=yind(i)<=size(y,2)

function z=indexedmult(x,xind,y,yind,checks,usebsxfun)
if nargin<6 || isempty(usebsxfun), usebsxfun=false; end
if nargin<5 || isempty(checks), checks=true; end
% check compatibility of arrays and index vectors
if checks, checkinputs(x,xind,y,yind); end

if exist('EVmergefunc','file') || usebsxfun
  z=EVmergefunc(x,uint64(xind),y,uint64(yind));
else
  ny=size(y,1);
  x=reshape(x,[size(x,1)/ny,ny,size(x,2)]);
  if ~isempty(xind), x=x(:,:,xind); end
  if ~isempty(yind), y=y(:,yind); end
  z=bsxfun(@times,x,reshape(y,[1 size(y)]));
end
z=squeeze(z);

function checkinputs(x,xind,y,yind)
if isempty(xind)
  if isempty(yind)
    if size(x,2)~=size(y,2)
      error('arrays are not column-compatible')
    end
  else
    if size(x,2)~=length(yind)
      error('# of columns in x must equal length(yind)')
    end
    if ~isindexvector(yind,size(y,2))
      error('yind contains invalid values')
    end
  end
else
  if isempty(yind)
    if length(xind)~=size(y,2)
      error('# of columns in y must equal length(xind)')
    end
    if ~isindexvector(xind,size(x,2))
      error('xind contains invalid values')
    end
  else
    if length(xind)~=length(yind)
      error('length(xind) must equal length(yind)')
    end
    if ~isindexvector(xind,size(x,2))
      error('xind contains invalid values')
    end
    if ~isindexvector(yind,size(y,2))
      error('yind contains invalid values')
    end
  end
end

% isindexvector Checks that a vector has integer values between 1 & n
function okay=isindexvector(ind,n)
okay=true;
if any(ind<1) || any(ind>n) || any(ind~=round(ind))
  okay=false;
end