%% demonstrates discrete ('d') type
p=[0.1 0.5;
   0.3 0.4;
   0.6 0.1]; 
e=rvdef('d',p); 
m=size(p,2); q=10; 
u=rand(q,1); 
ind=randi(m,q,1); 
x=e.simfunc(u,ind); 
disp('discrete type')
disp('  u   index   x')
fprintf('%5.3f  %2.0f    %2.0f\n',[u ind x]')


%% demonstrates binary ('b01') type
p=[0.1 0.9;
   0.9 0.1];
e=rvdef('b01',p);
m=size(p,2); q=10;
u=rand(q,1);
ind=randi(m,q,1);
x=e.simfunc(u,ind);
disp('binary type')
disp('  u   index   x')
fprintf('%5.3f  %2.0f    %2.0f\n',[u ind x]')



%% demonstrates multinomial logit ('logit') type
n=4;                       % # of categories
m=2;                       % # of conditioning variables
beta=randn(n,m+1);         % random coefficient matrix
e=rvdef('logit',beta);     % define the rv structure
q=10;                      % # of variates to generate
u=rand(q,1);               % generate a vector of uniform [0,1] variates
X={randn(q,1) randn(q,1)}; % generate random conditioning variable values
x=rvgen(q,e,X,u);          % generate the rv values
disp('multinomial logit type type')
disp('  u        X1       X2     x')
fprintf('%5.3f  %7.3f  %7.3f   %2.0f\n',[u X{:} x]')

%% demonstrates simple binomial ('bin') type
p=rand(1,1);
N=8;
e=rvdef('bin',p,N);  % alternately use e=rvdef('bin',{p,N}); 
q=10;
u=rand(q,1);
x=e.simfunc(u);
disp('simple binomial type')
disp('  u     x')
fprintf('%5.3f  %2.0f\n',[u x]')

%% demonstrates conditional binomial ('bin') type
p=rand(1,3);
N=[4 5];
e=rvdef('bin',{p,N});
u=rand(q,1); 
m=length(p)*length(N);
ind=randi(m,q,1); 
x=e.simfunc(u,ind);
disp('conditional binomial type')
disp('  u   index   x')
fprintf('%5.3f  %2.0f    %2.0f\n',[u ind x]')

%% demonstrates hypergeometric type type
N=5;
e=rvdef('hypgeo',N,N);
q=10;
K=randi(N,q,1);
n=randi(N,q,1);
u=rand(q,1); 
k=e.simfunc(u,K,n);
disp('hypergeometric type')
disp('  u     N   K   n   k')
fprintf('%5.3f  %2.0f  %2.0f  %2.0f  %2.0f\n',[u N+zeros(q,1) K n k]')

