% performs checks on the MDPSolve model variable
% USAGE
%   checkmdpmodel(model)
function checkmdpmodel(model)
P=model.P;
if isnumeric(P)
  cerr=max(abs(sum(P,1)-1));
  if  cerr > 5e-14   % columns don't sum to 1
    rerr = max(abs(sum(P,2)-1));
    if rerr > 5e-14 %    rows don't sum to 1
      fprintf('column sum divergence: %6.4e\n',cerr)
      fprintf('   row sum divergence: %6.4e\n',rerr)
      error('either columns or rows of P must sum to 1')
    else
      colstoch=false;
    end
  else
    rerr = max(abs(sum(P,2)-1));
    if rerr > 4e-14
      colstoch=true;
    end
  end
end
try
  R=model.R;
catch
  R=model.reward;
end
[ns,na]=size(R);
if na==1
  nx=ns;
  if ~isfield(model,'Ix')
    if ~(isfield(model,'X') && isfield(model,'svar'))
      error('Either X and svars or Ix must be defined if R is a vector')
    else
      X=model.X;
      svars=model.svars;
      [Ix,S]=getI(X,svars);
    end
  else
    Ix=model.Ix;
  end
  if length(Ix)~=nx
    error('Ix and R are not compatible')
  end
else
  nx=ns*na;
end

