Jmax=12;
reps=20;

ns=3^Jmax;
v=randn(ns,1);
timeres=zeros(7,5);
coptions=struct('forward',1,'transpose',1);
pp=cell(1,Jmax);
for c=3:9
  for i=1:reps
    for j=1:Jmax
      zz=randperm(9,c);
      pj=sparse([],[],[],3,3,c);
      pj(zz)=rand(c,1);
      pp{j}=normalizeP(pj);
      if c>=1, pp{j}=full(pp{j}); end
    end
    pp2={ckron(pp(1:2)), ckron(pp(3:4)), ckron(pp(5:6)), ckron(pp(7:8)), ckron(pp(9:10)), ckron(pp(11:12))};
    EV2=@(v)ckronx(pp2,v,coptions);
    pp3={ckron(pp(1:3)), ckron(pp(4:6)), ckron(pp(7:9)), ckron(pp(10:12))};
    EV3=@(v)ckronx(pp3,v,coptions);
    pp4={ckron(pp(1:4)), ckron(pp(5:8)), ckron(pp(9:12))};
    EV4=@(v)ckronx(pp4,v,coptions);
    if c<=7
      pp6={ckron(pp(1:6)), ckron(pp(7:12))};
      EV6=@(v)ckronx(pp6,v,coptions);
    end
    for j=1:Jmax
      pp{j}=full(pp{j});
    end
    EV1=@(v)ckronx(pp,v,coptions);

    tt=tic; x=EV1(v); timeres(c-2,1)=timeres(c-2,1)+toc(tt);
    tt=tic; x=EV2(v); timeres(c-2,2)=timeres(c-2,2)+toc(tt);
    tt=tic; x=EV3(v); timeres(c-2,3)=timeres(c-2,3)+toc(tt);
    tt=tic; x=EV4(v); timeres(c-2,4)=timeres(c-2,4)+toc(tt);
    if c<=7
      tt=tic; x=EV6(v); timeres(c-2,5)=timeres(c-2,5)+toc(tt);
    else
      timeres(c-2,5)=inf;
    end
  end
  fprintf('.')
end
fprintf('\n')
timeres