% Christine M. Hunter and Michael C. Runge. 
% The Importance of Environmental Variability and Management Control Error to Optimal Harvest Policies. 
% The Journal of Wildlife Management Vol. 68, No. 3, July, 2004, 585-594 

close all
clear variables

disp('Christine M. Hunter and Michael C. Runge') 
disp('The Importance of Environmental Variability and Management Control Error to Optimal Harvest Policies') 
disp('The Journal of Wildlife Management Vol. 68, No. 3, July, 2004, 585-594')


plots=true;

sas = 0.960; % Adult summer survival rate 
saw = 0.853; % Adult winter survival rate
sys = 0.930; % Yearling summer survival rate 
syw = 0.846; % Yearling winter survival rate 
sfs = 0.744; % Fawn summer survival rate 
sfw = 0.823; % Fawn winter survival rate 
k   = 0.903; % Proportion of yearlings to adults in the harvest
ma  = 0.75;  % maximal adult fecundity rate
la  = 19.3;  % lower bound in adult fecundity
ua  = 65.6;  % upper bound in adult fecundity
my  = 0.10;  % maximal yearling fecundity rate
ly  =  7.7;  % lower bound in yearling fecundity
uy  = 30.9;  % upper bound in yearling fecundity

sf=sfs*sfw;

% maximum population levels
maxA=60;
maxY=40;
inc=0.5;     % population increment
hinc=0.01; % harvest increment

a=(0:inc:maxA)';
y=(0:inc:maxY)';
h=(0:hinc:1)';

%%
stochastic=false;
nn=7;
% Adult female fecundity (female fawns per female) 
bta = @(Nt) min(ma,max((ma/(ua-la))*(ua - Nt),0));
% Yearling female fecundity 
bty = @(Nt) min(my,max((my/(uy-ly))*(uy - Nt),0));

fa=@(A,sas) A.*sas;  % fall adults
fy=@(Y,sys) Y.*sys;  % fall yearlings

% total harvest
H  = @(h,fa,fy) h.*(fa + fy);
% adults harvested
ha = @(h,fa,fy) (fa>(((1+k)/k./h-1).*fy)).*(H(h,fa,fy)-fy) + ...
                (fy>(((1+k)./h-1).*fa)).*fa  + ...
                (fa<=(((1+k)/k./h-1).*fy) & fy<=(((1+k)./h-1).*fa)).*(H(h,fa,fy)/(1+k));
% yearlings harvested
hy = @(h,fa,fy) (fa>(((1+k)/k./h-1).*fy)).*fy + ...
                (fy>(((1+k)./h-1).*fa)).*(H(h,fa,fy)-fa)  + ...
                (fa<=(((1+k)/k./h-1).*fy) & fy<=(((1+k)./h-1).*fa)).*(H(h,fa,fy)*(k/(1+k)));
Atran=@(fa,fy,h,saw,syw) (fa-ha(h,fa,fy)).*saw + (fy-hy(h,fa,fy)).*syw;
Ytran=@(F,sf) F.*sf;
fawns=@(A,Y) A.*bta(A+Y)+Y.*bty(A+Y);
TH=@(AH,YH) AH+YH;  % total harvest

if stochastic
  psys=rvdef('ne',[sys;(0.15*sys)],nn);
  psas=rvdef('ne',[sas;(0.15*sas)],nn);
  psf =rvdef('ne',[sf;(0.15*sf)],nn);
  psyw=rvdef('ne',[syw;(0.15*syw)],nn);
  psaw=rvdef('ne',[saw;(0.15*saw)],nn);
else
  psys=rvdef('ne',[sys;0],1);
  psas=rvdef('ne',[sas;0],1);
  psf =rvdef('ne',[sf; 0],1);
  psyw=rvdef('ne',[syw;0],1);
  psaw=rvdef('ne',[saw;0],1);  
end

D=add2diagram([],'adults',             's',1,{},                                               a);
D=add2diagram(D, 'yearlings',          's',1,{},                                               y);
D=add2diagram(D, 'harvest rate',       'a',1,{},                                               h);
D=add2diagram(D, 'fawns',              'c',0,{'adults','yearlings'},                           fawns);
D=add2diagram(D, 'sas',                'c',0,{},                                               psas);
D=add2diagram(D, 'sys',                'c',0,{},                                               psys);
D=add2diagram(D, 'fall adults',        'c',0,{'adults','sas'},                                 fa);
D=add2diagram(D, 'fall yearlings',     'c',0,{'yearlings','sys'},                              fy);
D=add2diagram(D, 'adults harvested',   'c',0,{'harvest rate','fall adults','fall yearlings'},  ha);
D=add2diagram(D, 'yearlings harvested','c',0,{'harvest rate','fall adults','fall yearlings'},  hy);
D=add2diagram(D, 'sf',                 'c',0,{},                                               psf);
D=add2diagram(D, 'saw',                'c',0,{},                                               psaw);
D=add2diagram(D, 'syw',                'c',0,{},                                               psyw);
D=add2diagram(D, 'adults+',            'f',1,{'fall adults','fall yearlings','harvest rate','saw','syw'},Atran);
D=add2diagram(D, 'yearlings+',         'f',1,{'fawns','sf'},                                   Ytran);
D=add2diagram(D, 'harvest',            'u',0,{'adults harvested','yearlings harvested'},       TH);

D.locs=[ ...
0.114 0.114 0.114 0.360 0.210 0.239 0.360 0.360 0.664 0.664 0.496 0.588 0.678 0.863 0.863 0.863;
0.635 0.449 0.898 0.449 0.220 0.100 0.319 0.180 0.319 0.180 0.920 0.920 0.920 0.635 0.449 0.250]';
D.attachments=[ ...
 1  2  1  5  2  6  3  7  8  3  7  8  7  8  3 12 13  4 11  9 10;
 4  4  7  7  8  8  9  9  9 10 10 10 14 14 14 14 14 15 15 16 16;
 5  5  4  5  4  5  4  5  5  4  5  5  6  6  5  4  4  5  4  5  5;
 8  1  8  2  8  2  8  1  1  8  1  1  2  2  1  8  8  1  8  1  1]';

figure(1); clf
set(gcf,'units','pixels','position',[600  500  1200 500])
drawdiagram(D)
text(0.05,0.1,'Hunter and Runge (2004)')

%%
t=cputime;
options=struct('ptype',1,'forcefull',0,'orderalg',1,'orderdisplay',-1,'orderonly',0,'print',1,'reps',0,'chunk',10000);
model=d2model(D,options);
fprintf('time taken setup using d2model: %6.3f\n',cputime-t)

if isnumeric(model.P)
  if  nnz(model.P)/numel(model.P)>0.35, model.P=full(model.P);
  else                                  model.P=sparse(model.P);
  end
  PP=model.P;
else
  EV=model.P;
end

model.d=1;
options=struct('vanish',0.999999,'print',2,'algorithm','i');
results=mdpsolve(model,options);
x=results.Ixopt;
X=model.X;

%%
cc=[0 0.2 0.4 0.6];
figure(2), clf
contour(y,a,reshape(X(x,1),length(y),length(a))',cc);
ylabel('adults')
xlabel('yearlings')
title('Optimal Harvest Rate')
%legend({'a=0','a=0.2','a=0.4','a=0.6'},'location','northeast')
