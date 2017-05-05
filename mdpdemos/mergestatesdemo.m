 clc 
 disp('Demonstrates use of mergestates')
 S1=[0 2;1 1;2 0]; P1=rand(5,3); P1=vxm(1./sum(P1,2),P1);
 A1=[0;1;0;0;1]; Ix1=[1;1;2;3;3];
 S2=[0;1]; P2=rand(2,2); P2=vxm(1./sum(P2,2),P2);
 [S,P,A,Ix]=mergestates({S1,P1,A1,Ix1},{S2,P2});
 disp('states')
 disp(S)
 disp('states (S), actions (A) and index (Ix)')
 disp([S(Ix,:) A Ix])
 P1,P2,P