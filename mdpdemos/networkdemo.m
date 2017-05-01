% network demo
%% setup the network
n=[10 5 8 5 10];
children=[1 3 5];
nc=length(children);
parents={[4 1], [5 3 4], 5};
% generate random conditional transition matrices
p=cell(1,nc);
for i=1:nc
  p{i}=normalizeP(rand(n(children(i)),prod(n(parents{i}))));
end
v=randn(prod(n(children)),1);

%% demonstration that the order or operations matters
order=[2 1 3];
EV1=networkgetEV(n,children,parents,p,order);
tic
ev1=EV1(v);
fprintf('time to process order 1:    %10.4f\n', toc)
order=[3 1 2];
EV2=networkgetEV(n,children,parents,p,order);
tic
ev2=EV2(v);
fprintf('time to process order 2:    %10.4f\n', toc)
fprintf('maximum difference in the result for the alternative orders: %1.4e\n',max(abs(ev2-ev1)))

%% compare the factored approach to computing and using the transition matrix
ns=prod(n(children));
nx=false(size(n)); 
nx([parents{:}])=true;
nx=prod(n(nx));
fprintf('size of the transition matrix: %1.0f by %1.0f (%1.0f elements)\n',[ns nx ns*nx])
k=memory; 
if k.MaxPossibleArrayBytes < nx*ns*8
  disp('transition matrix does not fit into memory')
else
  tic
  [P,X]=networkgetP(n,children,parents,p);
  fprintf('time to construct P matrix: %10.4f\n', toc)
  tic
  ev3=P'*v;
  fprintf('time to process P''v:        %10.4f\n', toc)
  fprintf('maximum difference in the result for order1 and P''v: %1.4e\n',max(abs(ev3-ev1)))
end