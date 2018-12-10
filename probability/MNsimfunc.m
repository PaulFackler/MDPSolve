% MNsimfunc Creates a simulation function for a Multinomial(p) distribution
% USAGE
%   simfunc = MNsimfunc(p);
% INPUT
%   p       : a d-vector of probability wieghts (summing to 1)
% OUPUT
%   simfunc : a function handle of the form x=simfunc(N) that accepts
%               a scalar positive integer and returns a 1xd vector of
%               non-negative integers summing to N with E[x]/N = p
% Example:
% N=3; d=5; p=normalizeP(rand(d,1))'; simfunc = MNsimfunc(p); reps=10000;  
% x=zeros(reps,d); for i=1:reps, x(i,:)=simfunc(N); end; disp([p;mean(x,1)/N])
% 
% Note that this simulator uses a random number of draws on a uniform 
% random number generator.
function simfunc = MNsimfunc(p)
p=p(:)';
pvec = p./([1,1-cumsum(p(1:end-1))]);
simfunc = @(N) MNsim(N,pvec);

function x = MNsim(N,pvec)
d = length(pvec);
x=zeros(1,d);
for i=1:d-1
  x(i) = sum( rand(N,1)<pvec(i) );
  N = N - x(i);
  if N==0, break; end
end
x(end) = N;