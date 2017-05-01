% Christine M. Hunter and Michael C. Runge. 
% The Importance of Environmental Variability and Management Control Error to Optimal Harvest Policies. 
% The Journal of Wildlife Management Vol. 68, No. 3, July, 2004, 585-594 

clc
disp('Christine M. Hunter and Michael C. Runge') 
disp('The Importance of Environmental Variability and Management Control Error to Optimal Harvest Policies') 
disp('The Journal of Wildlife Management Vol. 68, No. 3, July, 2004, 585-594')
disp('Determinstic Model')
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

% maximum population levels
maxA=60;
maxY=40;
inc=1;     % population increment
hinc=0.02; % harvest increment

H    =@(A,Y,h) h.*(A*sas + Y*sys);
% Adult female fecundity (female fawns per female) 
bta = @(Nt) min(ma,max((ma/(ua-la))*(ua - Nt),0));
% Yearling female fecundity 
bty = @(Nt) min(my,max((my/(uy-ly))*(uy - Nt),0));

ha =@(A,Y,h) ((A*sas)>(((1+k)/k./h-1).*Y*sys)).*(H(A,Y,h)-Y*sys) + ...
             ((Y*sys)>(((1+k)./h-1).*A*sas)).*A*sas  + ...
             ((A*sas)<=(((1+k)/k./h-1).*Y*sys) & (Y*sys)<=(((1+k)./h-1).*A*sas)).*(H(A,Y,h)/(1+k));
hy =@(A,Y,h) ((A*sas)>(((1+k)/k./h-1).*Y*sys)).*Y*sys + ...
             ((Y*sys)>(((1+k)./h-1).*A*sas)).*(H(A,Y,h)-A*sas)  + ...
             ((A*sas)<=(((1+k)/k./h-1).*Y*sys) & (Y*sys)<=(((1+k)./h-1).*A*sas)).*(H(A,Y,h)*(k/(1+k)));
Atran=@(A,Y,h) (A*sas-ha(A,Y,h))*saw + (Y*sys-hy(A,Y,h))*syw;
Ytran=@(A,Y,h) (A.*bta(A+Y)+Y.*bty(A+Y))*(sfs*sfw);


fawns=@(A,Y) A.*bta(A+Y)+Y.*bty(A+Y);


a=(0:inc:maxA)';
y=(0:inc:maxY)';
h=(0:hinc:1)';

D=add2diagram([],'adults', 's',a,{},[],'l',[.2 .8]);
D=add2diagram(D,'yearlings', 's',y,{},[],'l',[.2 .5]);
D=add2diagram(D,'Harvest rate', 'a',h,{},[],'l',[.2 .2]);
D=add2diagram([],'fawns', 's',[],{'adults','yearlings'},fawns,'l',[.4 .65]);

figure(1); clf
drawdiagram(D)