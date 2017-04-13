s=[0;1];
a=[1;2];
P=normalizeP(rand(2,4));

D=add2diagram([],'State',  's',1,{},                s);
D=add2diagram(D, 'Action', 'a',1,{},                a);
D=add2diagram(D, 'State+', 'f',1,{'State','Action'},rvdef('d',P,s));
D=add2diagram(D, 'Utility','u',1,'State',           @(S)S);

figure(1); clf
drawdiagram(D,struct('nodecolor',[0.85 0.85 0.85],'backgroundcolor',[1 1 1]))

Pfunc1=@(S,A) P(:,gridmatch({S,A},{s,a}));
Pfunc2=@(S,A) P(:,A+2*S);
X=dvalues(D,{'State','Action'});
disp('Three ways to obtain the transition matrix')
disp(P)
disp(' ')
disp(Pfunc1(X{:}))
disp(' ')
disp(Pfunc2(X{:}))
