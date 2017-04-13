% Lubow's coin toss game as an adaptive management problem
% This problem was described in
% "Introduction to Adaptive Stochastic Dynamic Programing Theory
%  for Adaptive Resource Management" by Bruce C. Lubow, pp.22-31.
clc
disp('Coin toss game')
p=1000;
delta=0.95;
T = inf;
P1=[.6 .6;   % states are outcomes of the toss (1) heads (2) tails
    .4 .4];  % P1 is for model 1
P2=[.4 .4;   % P2 is for model 2
    .6 .6];
R=[ 1;
   -1];

 % Define a belief state for models 1 and 2
[b,Pb,Rb,Sb]=amdp(p,{P1,P2},R,[1;2]);

% add a stopped state to the problem
n=size(Pb,1);
Pb=[Pb sparse(n,n+2);sparse(1,n) ones(1,n+2)];
Rb=[[Rb;-inf] [Rb;0]];
Sb=[Sb;3 0 0];
Xb=[zeros(n+1,1);ones(n+1,1)];

clear model
model.discount=delta;
model.P=Pb;
model.R=Rb;
model.T=T;
clear options
options.algorithm='p';
results=mdpsolve(model,options);
v=results.v; x=results.Ixopt; pstar=results.pstar;


figure(1); clf
plot(Sb(1:end-1,2),Xb(x(1:end-1)),'*')
title('Optimal Action (0=Play, 1=Stop)')
xlabel('Belief Weight that Coin is Biased Towards Heads')
set(gca,'ytick',[0 1],'ylim',[-.25 1.25])

if 1-delta<0.00001
  figure(2);  clf
  plot(Sb(1:p+1,2),(1-delta)*v(1:p+1))
  title('Average Value Function')
else
  figure(2);  clf
  plot(Sb(1:p+1,2),reshape(v(1:end-1),p+1,2))
  title('Value Function')
  legend({'Current Toss is Heads','Current Toss is Tails'},'location','northwest')
end
xlabel('Belief Weight that Coin is Biased Towards Heads')

if p<=10
[plr,F]=longrunP(pstar);
disp('Coin State, Model Weights and Longrun Probabilities')
disp([Sb(:,[1 2]) plr])
disp(' ')
disp('Coin State, Model Weights and Reachability Probabilities')
disp([Sb(:,[1 2]) abs(F(:,1:p+1))])
end