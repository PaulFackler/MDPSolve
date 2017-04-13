function n=size(a,dim)
if nargin>1
  if dim==1
    n=a.m;
  elseif dim==2
    n=a.s(end);
  else
    n=1;
  end
else
  n=[a.m a.s(end)];
end