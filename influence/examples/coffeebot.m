% coffee robot example with perfect observability
% Boutilier, C., R. Dearden  and M. Goldszmidt
% “Stochastic Dynamic Programming with Factored Representations”
% Artificial Intelligence, 121: 49–107
2000.
[BRY 86] BRYANT
x=[0;1];
sI=speye(2);

pDrops=0.2; % probability coffee drops at delivery

% actions
A=(0:4)';
DoNothing = 0;
Go        = 1;
Buy       = 2;
Del       = 3;
GetU      = 4;

Wet =@(Wet,Rain,HasU,A)       Wet | (Rain & ~HasU & A==Go);
HasU=@(HasU,InO,A)            HasU | (InO & A==GetU);
InO =@(InO,A)                 (InO & A~=Go) | (~InO & A==Go);
HCU =@(HCU,HCR,InO,A,Drops)   HCU | (HCR & InO & A==Del & ~Drops);
HCR =@(HCR,InO,A)             (HCR | A==Buy) & A~=Del;

Ufunc    = @(HCU,Wet) 0.9*(HCU) + 0.1*(~Wet);
feasible = @(A,HasU,InO,HCU,HCR) (A==DoNothing) | (A==Go) | (A==Buy & InO==0 & HCR==0) | ...
           (A==Del & InO==1 & HCR==1 & HCU==0) | (A==GetU & InO==1 & HasU==0);
Utility  = @(A,Wet,HasU,InO,HCU,HCR) ...
           ifthenelse(feasible(A,HasU,InO,HCU,HCR),Ufunc(HCU,Wet),-inf);

%%
D=[];
D=add2diagram(D, 'A',      'a',1,{},                             A);
D=add2diagram(D, 'HCU',    's',1,{},                             x);
D=add2diagram(D, 'HCR',    's',1,{},                             x);
D=add2diagram(D, 'InO',    's',1,{},                             x);
D=add2diagram(D, 'Wet',    's',1,{},                             x);
D=add2diagram(D, 'Rain',   's',1,{},                             x);
D=add2diagram(D, 'HasU',   's',1,{},                             x);
D=add2diagram(D, 'Drops',  'c',1,{},                             rvdef('d',[pDrops;1-pDrops],x));
D=add2diagram(D, 'HCU+',   'f',1,{'HCU','HCR','InO','A','Drops'},HCU);
D=add2diagram(D, 'HCR+',   'f',1,{'HCR','InO','A'},              HCR);
D=add2diagram(D, 'InO+',   'f',1,{'InO','A'},                    InO);
D=add2diagram(D, 'Wet+',   'f',1,{'Wet','Rain','HasU','A'},      Wet);
D=add2diagram(D, 'Rain+',  'f',1,{'Rain'},                       rvdef('d',sI,x));
D=add2diagram(D, 'HasU+',  'f',1,{'HasU','InO','A'},             HasU);
D=add2diagram(D, 'Utility','u',1,{'A','Wet','HasU','InO','HCU','HCR'}, Utility);
D.locs=[ ...
0.261 0.261 0.261 0.261 0.261 0.261 0.261 0.329 0.712 0.712 0.712 0.712 0.712 0.712 0.712;
0.916 0.789 0.689 0.589 0.490 0.391 0.291 0.144 0.789 0.689 0.589 0.490 0.391 0.291 0.090]';
D.attachments=[ ...
 3  2  4  7  5  1  1  4  7  6  1  7  6  5  1  4  1  4  3  8  1  4  3  2;
15 15 15 15 15 15 14 14 14 13 12 12 12 12 11 11 10 10 10  9  9  9  9  9;
 4  4  4  4  4  4  5  5  5  5  5  5  5  5  5  5  5  5  5  6  5  5  5  5;
 1  1  1  1  1  1  8  1  1  1  8  1  1  1  8  1  8  1  1  2  8  1  1  1]';

%%
figure(1); clf
set(gcf,'units','normalized','position',[.45  .05 .45  .85])
drawdiagram(D)

%%
t=cputime;
options=struct('ptype',1,'forcefull',0,'orderalg',1,'orderdisplay',-1,'reps',0);
model=d2model(D,options);
fprintf('time taken to set up model using d2model: %6.3f\n',cputime-t)
model.d=1;
options=struct('algorithm','i','vanish',0.999999);
results=mdpsolve(model,options);
Xopt=model.X(results.Ixopt,:);

%%
%disp('columns of X:')
%disp(' 1) action 2) Wet 3) HasU 4) Rain 5) InO 6) HCU 7) HCR')

disp(' 1) HCU 2) HCR 3) InO 4) Wet 5) Rain 6) HasU 7) action ')

Xopt(:,[2:end 1])
return
disp('optimal policy in office and not raining')
Xopt(Xopt(:,6)==0 & Xopt(:,5)==1 & Xopt(:,4)==0,:)
disp('optimal policy in office and raining')
Xopt(Xopt(:,6)==0 & Xopt(:,5)==1 & Xopt(:,4)==1,:)
disp('optimal policy in cafe and not raining')
Xopt(Xopt(:,6)==0 & Xopt(:,5)==0 & Xopt(:,4)==0,:)
disp('optimal policy in cafe and raining')
Xopt(Xopt(:,6)==0 & Xopt(:,5)==0 & Xopt(:,4)==1,:)

