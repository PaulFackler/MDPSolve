% creates a simple Markov process and obtains the long run distribution
% analytically and using Monte Carlo simulation

close all

n=10;
S=(1:n)';
ps=normalizeP(rand(n,n));

D=[];
D=add2diagram(D,'s','s',1,{},S);
D=add2diagram(D,'s+','f',1,{'s'},rvdef('d',ps,S));

%%
tt1=tic;
lrpa=longrunP(ps);               % analytic long run probability
tt1=toc(tt1);
% simulate process for T periods and get the empirical distribution 
% of the period T values
tt2=tic;
S0=1;
rep=10000;
T=100;
X=dsim(D,S0,rep,T);
lrps=histc(X{2}(:,end),S)/rep;   % simulated long run probability
tt2=toc(tt2);

disp('Markov transition matrix')
disp(ps)

disp('analytic and simulated long run distribution')
disp([S lrpa lrps])

disp('time taken to compute long run probabilities')
disp([tt1 tt2])