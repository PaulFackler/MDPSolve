% EVgreedyorder Greedy algorithm to determine order for an EV function
% USAGE
%   [I,cost] = EVgreedyorder(q,n,m);
% INPUTS
%   n  : d-vector of state variable sizes
%   q  : d-element cell array with parent vectors for each state variable
%          note: noise variables appear as negative integers (see EVcreate)
%   mx : vector of conditioning variable sizes
%   me : vector of noise variable sizes
% OUTPUT
%   I    : order vector
%   cost : number of flops used to compute the EV function

function [I,cost] = EVgreedyorder(n,q,mx,me)
if nargin < 4; me = []; end
if isempty(me)
  % check forward and backward versions 
  [I,cost]   = EVgreedyorderF(q,n,mx);
  [Ib,costb] = EVgreedyorderB(q,n,mx);
else % model includes noise terms
  d = length(n);
  dimx = length(mx);
  eused = zeros(1,length(me));
  for i = 1:d
    ind = q{i} < 0; 
    qe = -q{i}(ind);
    eused(qe) = eused(qe) + 1;
    q{i}(ind) = dimx + qe;
  end
  m = [mx(:)' me(:)'];
  % check forward and backward versions 
  [I,cost]   = EVgreedyorderFe(q,n,m,eused,dimx);
  [Ib,costb] = EVgreedyorderBe(q,n,m,eused,dimx);
end
% pick the best of the forward and backwards orders
if costb < cost, cost = costb; I = Ib; end

% forward version with noise - adds the variable that increases
% the cost the least to the next in the order
function [I,cost]=EVgreedyorderFe(q,n,m,eused,dimx)
d = length(q);
I = zeros(1,d);
combined = [];
remaining = 1:d;
N = prod(n);
cost = 0;
for i = 1:d
  wi = inf;
  for j = 1:d-i+1
    var = remaining(j);
    combinedj = union(combined,q{var});
    wj = N*prod(m(combinedj));
    if wi > wj
      jopt = j;
      wi = wj;
    end
  end
  % I(i) is the variable entering at ith step
  I(i) = remaining(jopt);
  qi = q{I(i)};
  combined = union(combined,qi);
  N = N/n(I(i));
  cost = cost + wi;
  remaining(jopt) = [];
  % determine any noise terms that get summed out 
  % and eliminate these from combined
  ii = qi(qi>dimx) - dimx;
  eused(ii) = eused(ii)-1;
  ii = find(eused==0) + dimx;
  combined(ismember(combined,ii)) = []; 
end

% backward version with noise - adds the variable that increases
% the cost the most to the end of the order
function [I,cost] = EVgreedyorderBe(q,n,m,eused,dimx)
I = []; cost = inf; return % not currently implemented
d = length(q);
I = zeros(1,d);
remaining = 1:d;
N = 1;
cost = 0;
for i = d:-1:1
  wi = 0;
  combined = unique([q{remaining}]);
  M = prod(m(combined));
  for j = 1:i
    rj = remaining; rj(j) = [];
    combined = unique([q{rj}]);
    wj = M/prod(m(combined));
    if wi < wj
      jopt = j;
      wi = wj;
    end
  end
  I(i) = remaining(jopt);
  N = N*n(I(i));
  cost = cost + N*M;
  remaining(jopt) = [];
end


% forward version - adds the variable that increases
% the cost the least to the next in the order
function [I,cost] = EVgreedyorderF(q,n,mx)
d = length(q);
I = zeros(1,d);
combined = [];
remaining = 1:d;
N = prod(n);
cost = 0;
for i = 1:d
  wi = inf;
  for j = 1:d-i+1
    var = remaining(j);
    combinedj = union(combined,q{var});
    wj = N*prod(mx(combinedj));
    if wi > wj
      jopt = j;
      wi = wj;
    end
  end
  I(i) = remaining(jopt);
  combined = union(combined,q{I(i)});
  N = N/n(I(i));
  cost = cost + wi;
  remaining(jopt) = [];
end

% backward version - adds the variable that increases
% the cost the most to the end of the order
function [I,cost] = EVgreedyorderB(q,n,mx)
d = length(q);
I = zeros(1,d);
remaining = 1:d;
N = 1;
cost = 0;
for i = d:-1:1
  wi = 0;
  combined = unique([q{remaining}]);
  M = prod(mx(combined));
  for j = 1:i
    rj = remaining; rj(j) = [];
    combined = unique([q{rj}]);
    wj = M/prod(mx(combined));
    if wi < wj
      jopt = j;
      wi = wj;
    end
  end
  I(i) = remaining(jopt);
  N = N*n(I(i));
  cost = cost + N*M;
  remaining(jopt) = [];
end