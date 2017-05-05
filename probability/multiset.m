% multiset Multiset coefficient ((q,p))
% USAGE
%   m=multiset(q,p);
% m=(q+p-1)!/(q!(p-1)!)
function m=multiset(q,p)
  m=1:p+1;
  for i=1:q-2
    m=cumsum(m);
  end
  m=m(end);
end
