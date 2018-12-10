function [order,mergevec,c]=EVoptorderMC(p,parents,X,e)
m=length(parents);
itmax=50;
best = inf;
order = 1:m;
bestorder = order;  
itfail = 0;
while 1
  itfail = itfail +1;
  [mergevec,c] = EVoptgroups(p(order),parents(order),X,e);
  if c < best
    best = c;
    bestorder = order;
    itfail = 0;
  end
  mm=nan(1,m); mm(1:length(mergevec))=mergevec;
  fprintf('%1.0f  %1.0f  %1.0f       %1.0f  %1.0f  %1.0f       %16.0f    %16.0f \n',order, mm, c, best)
  if itfail > itmax, break; end
  if length(mergevec)==1
    neworder = randperm(m);
    while isequal(order,neworder)
      neworder = randperm(m);
    end
    order = neworder;
  else
    k = length(mergevec);
    neworder = cell(1,k);
    g = randperm(k);
    while isequal(g,1:k)
      g = randperm(k);
    end
    j = 0;
    for i = 1:k
      neworder{g(i)} = order(j+1:j+mergevec(i));
      j = j + mergevec(i);
    end
    order = [neworder{:}];
  end
end
order = bestorder;
[mergevec,c] = EVoptgroups(p(order),parents(order),X,e);      