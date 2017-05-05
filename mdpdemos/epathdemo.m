 % illustrates the use of epath
% as well as catcountP
clc
disp('Demonstrating the use of epath to compute expected time paths')
p=[0.8 0.1;
   0.2 0.9];
N=10;
T=25;
S=(N:-1:0)';
P=catcountP(N,2,2,p);
Ef=epath(S,P,[],T,1);

figure(1);
plot(0:T,Ef,'k')
title('Expected Time Paths of the # of Successes')
xlabel('time')
ylabel('E[# of successes]')