% Transition functions for Anderson's duck harvesting model
% USAGE
%   Splus=Anderson75tran(Nn,Pn,Dn,Hn,additive);
% INPUTS
%   Nn : population size (state 1)
%   Pn : pond numbers    (state 2)
%   Dn : harvest size in birds (action)
%   Hn : harvest random variable
%   additive : 0/1  (0) compensatory model (1) additive model
% OUTPUTS
%   Splus: 2-column matrix of values of next period's state: [N+ P+]
function Nplus=Anderson75tran(Nn,Pn,Dn,hshock,additive)
Yn    = 1./( (1./(12.48*Pn.^0.851))+(0.519./Nn) ) ;   % Eqn. 3
Fn    = 0.92*Nn + Yn;                                 % Eqn. 5 
%Dn = min(Nn,hshock.*Dn);
Dn = hshock.*Dn;
Hn    = Dn./Fn ;                                      % p 1286 
Hn(Fn==0)=0;
Hn(Hn>1)=1;
if additive                                           % Additive model
  Nplus = Nn.*(1-0.37*exp(2.78*Hn)) + ...
          Yn.*(1-0.49*exp(0.90*Hn));                  % Eqns.7,8 
  %Nplus = Nn.*(1-0.27*exp(2.08*Hn)) + ...
  %        Yn.*(1-0.40*exp(0.67*Hn));                  % Eqns.9,10 
else                                                  % Compensatory model 
  Nplus = (Nn*0.57 + Yn*0.50).*(Hn < 0.25) + ...      % p 1292 
          (max(0, Nn.*(0.57-1.2*(Hn-0.25))) + ...
          (max(0, Yn.*(0.50-1.0*(Hn-0.25))))).* ...
          (Hn >= 0.25);                               % Eqn.14 
end
Nplus(Fn<Dn)=0;
Nplus=max(0,Nplus);
