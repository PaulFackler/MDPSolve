function c=mtimes(a,b)
if isa(b,'cellsparse')
  if isa(a,'cellsparse')
    error('multiplication of two cellsparse objects is not supported.')
  end
  if size(a,2)~=b.m
    error('Inner matrix dimensions must agree.')
  end
  c=zeros(size(a,1),b.s(end));
  for j=1:b.n
    c(:,b.s(j)+1:b.s(j+1))=a*b.data{j};
  end
else
  if size(b,1)~=a.s(end)
    error('Inner matrix dimensions must agree.')
  end
  c=zeros(a.m,size(b,2));
  for j=1:a.n
    c=c+a.data{j}*b(a.s(j)+1:a.s(j+1),:);
  end    
end