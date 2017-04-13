close all
clear variables

r       = 0.1;                                 % intrinsic growth rate
K       = 1;                                   % carrying capacity
sigma   = 0.02;                                % noise standard deviation 
delta   = 0.98;                                % discount factor
Snext   = @(E,z) (1+r)*E./(1+(r/K)*E).*exp(z); % next period state
Utility = @(S,E) ifthenelse(E<=S,S-E,-inf);    % utility function
Smax    = 1.5*K;                               % maximum state value
ns      = 100;                                 % number of state values
ne      = 10;                                  % number of noise values
S       = linspace(0,Smax,ns)';                % population values
pe      = rvdef('ne',[-sigma^2/2;sigma],ne);   % noise distribution 
D=add2diagram([],'population', 's',1,{}                         ,S);
D=add2diagram(D ,'escapement', 'a',1,{}                         ,S);
D=add2diagram(D ,'noise',      'c',1,{}                         ,pe);
D=add2diagram(D ,'population+','f',1,{'escapement','noise'}     ,Snext);
D=add2diagram(D ,'Utility',    'r',1,{'population','escapement'},Utility);



%%
t=cputime;
options=struct('ptype',1,'forcefull',0,'orderalg',0,'cleanup',0,'orderdisplay',-1,'reps',1000);
model=d2model(D,options);
fprintf('time taken to set up model using d2model: %6.3f\n',cputime-t)

if isnumeric(model.P)
  if  nnz(model.P)/numel(model.P)>0.35, model.P=full(model.P);
  else                                  model.P=sparse(model.P);
  end
  PP=model.P;
else
  EV=model.P;
end

model.d=0.98;

options.algorithm='i'; options.print=1;
results=mdpsolve(model,options);
Xopt=model.X(results.Ixopt,:);

%%
D.locs=[ ...
0.3  0.3 0.5  0.7  0.7;
0.65 0.5 0.25 0.65 0.4]';
D.attachments=[ ...
 2  1  3  2;
 5  5  4  4;
 5  5  7  5;
 1  1  1  1]';
figure(1); clf
drawdiagram(D);

%%
figure(2)
plot(Xopt(:,2),Xopt(:,1),'k')
xlabel('population','fontname','Rockwell','Fontsize',12)
ylabel('optimal escapement','fontname','Rockwell','Fontsize',12)


