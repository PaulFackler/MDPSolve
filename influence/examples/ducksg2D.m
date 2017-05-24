clearclose
disp('Anderson (1975) duck harvest model using g2EV')

% model parameters
delta  = 0.98;               % discount factor 
mu     = 16.46;              % rainfall mean
sigma2 = 4.41;               % rainfall variance
% pond transition function
Pplus = @(Pn,Rn) -2.76 + 0.391*Pn + 0.233*Rn;  
% population transition function
Nplus = @(N,P,D,Hn) Anderson75tran(N,P,D,Hn,0);

% harvest shock
hshock=[0.9;1;1.1];
hprob=[0.2;0.6;0.2];

% solution parameters
nmin=2;   nmax=18;    ninc=0.1;
pmin=0.5; pmax=3.5;   pinc=0.1;
dmin=0;   dmax=12;    dinc=0.05;
nr=11;          % number of random rainfall noise values
cleanup=2;      % project to boundary

% get variable values and indices
rrv=rvdef('n',[mu;sqrt(sigma2)],nr);
hrv=rvdef('d',hprob,hshock);
e={rrv,hrv};

n=(nmin:ninc:nmax)';     % population values
p=(pmin:pinc:pmax)';     % pond number values
d=(dmin:dinc:dmax)';     % harvest values
nn=length(n); np=length(p); nd=length(d);

xelist={[1 2 3 -2],[2 -1]};
D=g2D({Nplus Pplus},{n,p},{n,p,d},e,xelist,1);