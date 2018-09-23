% EVcombine Creates an EV function with both a full and indexed evaluation form
% USAGE
%   EV=EVcombine(EV0,EV1);
% INPUTS
%   EV0 : a function handle for a function of the form EV0=@(V)
%   EV1 : a function handle for a function of the form EV1=@(V,I)
% OUTPUT
%   EV  : a function handle for a function that calls EV0 if passed
%           a single input and calls EV1 if passed two inputs
%
% Example:
%   P=normalizeP(rand(3,6));
%   EV0=@(V) P'*V; 
%   EV1=@(V,I) P(:,I)'*V; 
%   EV=EVcombine(EV0,EV1);
%   EV(randn(3,1))
%   v=randn(3,1);
%   ev0=EV(v);
%   I=[1;3;6];
%   ev1=EV(v,I);
function EV=EVcombine(EV0,EV1)
  if nargin==2
    EV=@EVeval2;
  else
    EV=@EVeval1;
  end
  function v=EVeval1(V,I)
    v=EV0(V);
    if nargin==2, v=v(I); end
  end
  function v=EVeval2(V,I)
    if nargin==1
      v=EV0(V);
    else
      v=EV1(V,I);
    end
  end
end
