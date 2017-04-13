% Lubow's coin toss game as a POMDP
% This problem was described in
% "Introduction to Adaptive Stochastic Dynamic Programing Theory
%  for Adaptive Resource Management" by Bruce C. Lubow, pp.22-31.
clc
disp('Coin toss game')
% states are (1) coin is biased to heads & (2) coin is biased to tails 
p=1000;       % # number of belief state intervals
delta=0.95;   % discount factor
P=[1 0        % state does not change
   0 1]; 
Q=[.6 .4;     % Information variable is the outcome of the current flip
   .4 .6];
R=[1 -1;      % reward is based on the outcome of the current flip
   1 -1];

% this creates a belief state model
options=struct('Qtype',1,'Rtype',3);
[b,Pb,Rb]=pomdp(p,P,Q,R,options);

% supplement the model with a stopped state
% define action to be (1) continue (2) stop
n=size(b,1);
Pb=[Pb sparse(n,n+2);zeros(1,n) ones(1,n+2)];
Rb=[Rb zeros(n,1);0 0];
Xb=[zeros(n+1,1);ones(n+1,1)];

clear model
model.discount=delta;
model.P=Pb;
model.R=Rb;
results=mdpsolve(model,options);
v=results.v; x=results.Ixopt; pstar=results.pstar;


figure(1);  clf
plot(b(:,1),Xb(x(1:end-1)),'*')
title('Optimal Action (0=Play, 1=Stop)')
xlabel('Belief Weight that Coin is Biased Towards Heads')
set(gca,'ytick',[0 1],'ylim',[-.25 1.25])


if 1-delta<0.00001
  figure(2);  clf
  plot(b(:,1),(1-delta)*v(1:end-1))
  title('Average Value Function')
else
  figure(2);  clf
  plot(b(:,1),v(1:end-1))
  title('Value Function')
end
xlabel('Belief Weight that Coin is Biased Towards Heads')


if p<=10
disp('State and Longrun probabilities')
disp([[b(:,1);0] longrunP(pstar)])
end