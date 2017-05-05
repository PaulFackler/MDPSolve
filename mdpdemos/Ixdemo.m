% demonstrates the effect of alternative variable orderings
clc
disp('Demonstrates the effect of alternative variable orderings')
S=(0:3)'; A=[0;1]; 
ns=length(S); na=length(A); 
X=rectgrid(A,S); 
Ix=(1:ns)'*ones(1,na); Ix=Ix(:);
disp('Action is first in order')
disp('     j    Ix     A     S')
disp([(1:ns*na)' Ix X]);
%textable([(1:ns*na)' Ix X],0,{'$j$','$I^x$','$A$','$S$'})


disp('Action is second in order')
disp('     j    Ix     A     S')
X=rectgrid(S,A); 
Ix=ones(na,1)*(1:ns); Ix=Ix(:);
disp([(1:ns*na)' Ix X]);
%textable([(1:ns*na)' Ix X],0,{'$j$','$I^x$','$S$','$A$'})