function P=pests_catcountfunc(Xj,b0u,b0t,b1u,b1t,P1,P2,N)
 cu=1./(1+exp(b0u+b1u*(Xj(:)/N))); 
 ct=1./(1+exp(b0t+b1t*(Xj(:)/N)));
 P1([1 2],1)=[1-cu cu];
 P2([1 2],1)=[1-ct ct];
 P=[P1 P2];