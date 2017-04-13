% noisy states
clc
disp('Discretize S and display transition matrix')
s=linspace(0,1,11)';
y=s;
e=[0.9;1;1.1];
w=[1;1;1]/3;
options=struct('cleanup',0);
Q=g2P(@(x,e) x.*e,y,s,e,w,options);
disp('S+ = S*e')
full(Q)
%textable(full(Q)',3)


e=[-0.1;0;0.1];
Q=g2P(@(x,e) x+e,y,s,e,w,options);
disp('S+ = S+e (low variance)')
full(Q)
%textable(full(Q)',3)

e=[-0.5;0;0.5];
Q=g2P(@(x,e) x+e,y,s,e,w,options);
disp('S+ = S+e (moderate variance)')
full(Q)
%textable(full(Q)',3)


e=[-1.5;0;1.5];
Q=g2P(@(x,e) x+e,y,s,e,w,options);
disp('S+ = S+e (high variance)')
full(Q)
%textable(full(Q)',3)