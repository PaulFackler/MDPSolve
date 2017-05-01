% checkQ determines whether 1 message service dominates another
% USAGE
%   res=checkQ(p,Q1,Q2);
% INPUTS
%   p : n-element stochastic vector
%   Q1 : m1xn column stochastic matrix
%   Q2 : m2xn column stochastic matrix
% OUTPUT
%   res : 0 if neither is dominant
%         1 if 1 is dominant
%         2 if 2 is dominant
function res=checkQ(p,Q1,Q2)
q1=Q1*p; q2=Q2*p;
P1=diag(q1)\(Q1*diag(p1));
P2=diag(q2)\(Q2diag(p2));

m1=size(P1,1);
m2=size(P2,1);

A=[kron(eye(m1),P2);kron(q1',eye(m2);kron(eye(m1),ones(1,m2))];
b=[P1(:);q2;ones(m1,1)];
A=A\b;

if all(A(:)>=0), res=2; return; end

A=[kron(eye(m2),P1);kron(q2',eye(m1);kron(eye(m2),ones(1,m1))];
b=[P2(:);q1;ones(m2,1)];
A=A\b;

if all(A(:)>=0), res=1; return; end

res=0;