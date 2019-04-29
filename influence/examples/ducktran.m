% duck transition function (see ducks)
function Nplus=ducktran(Nn,Pn,Dn)
Yn    = 1./( (1./(12.48*Pn.^0.851))+(0.519./Nn) ) ;   % Eqn. 3
Fn    = 0.92*Nn + Yn;                                 % Eqn. 5 
Hn    = Dn./Fn ;                                      % p 1286 
Hn(Fn==0)=0;
Hn(Hn>1)=1;
% Additive model
Nplus = Nn.*(1-0.37*exp(2.78*Hn)) + ...
          Yn.*(1-0.49*exp(0.90*Hn));                  % Eqns.7,8 
