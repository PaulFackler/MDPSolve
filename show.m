% SHOW Displays a matrix x with a specified number of decimal places 
% USAGE
%   show(x,d)
% INPUTS
%   x - a matrix of numbers
%   d - a vector indicating the # of decimal places for each column of x
% d is optional, if omitted 0 is used for columns of integers and 4 is used 
% otherwise.  2 blank spaces are used to separate columns.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function show(x,d)

space='  ';
x=squeeze(x);
if length(size(x))>2
  error('Not defined for multidimensional arrays')
end

x=full(x);

if nargin<2 | isempty(d)
  d=4*(any(x~=fix(x)));
end
n=size(x,2);
if size(d,2)<n
  if length(d)==1, d=d+zeros(1,n);
  else,  d=[d zeros(1,n-size(d,2))];
  end
end

fmt=[];
for i=1:n
  if size(x,1)>1
    m=max(x(x(:,i)<inf & ~isnan(x(:,i)),i));
    s=min(x(x(:,i)>-inf & ~isnan(x(:,i)),i))<0;
  else
    m=x(1,i);
    s=x(1,i)<0;
  end
  m=fix(log(m+eps)./log(10))+1;
  temp=real(m+s+d(i)+(d(i)>0));
  fmt=[fmt  space '%' num2str(temp) '.' num2str(d(i)) 'f'];
end
fmt=[fmt '\n'];
fprintf(fmt,x')