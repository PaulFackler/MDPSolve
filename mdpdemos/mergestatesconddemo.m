function mergestatesconddemo
 clc
 disp('Demonstrates use of mergestates - 8 forms of the combined state are displayed')

 S1=[0 2;1 1;2 0]; 
 P1m{1}=rand(5,3); P1m{1}=vxm(1./sum(P1m{1},2),P1m{1});
 P1m{2}=rand(5,3); P1m{2}=vxm(1./sum(P1m{2},2),P1m{2});
 A1=[0;1;0;0;1]; Ix1=[1;1;2;3;3];
 S2=[0;1]; 
 P2 =rand(2,2); P2(1,2)=0; P2=vxm(1./sum(P2,2),P2);
 P2s=sparse(P2);
 P1f =@(S2) P1m{S2+1};
 P1ft=@(S2) P1m{S2+1}'; 
 options=struct('transposed',1);
 
 disp('Prob(S1|S2=0)')
 disp(P1m{1})
 disp('Prob(S1|S2=1)')
 disp(P1m{2})
 disp('Prob(S2)')
 disp(P2)
 
 % regular form - S1 first
 disp('S1 first in order')
 disp('regular form')
 
 [S,P,A,Ix]=mergestates({S1,P1f,A1,Ix1},{S2,P2});
 disp('P full')
 showres
 %textable([NaN+zeros(3,4) S'],0)
 %textable([S(Ix,:) A P],[0 0 0 0 2 2 2 2 2 2])
 
 [S,P,A,Ix]=mergestates({S1,P1f,A1,Ix1},{S2,P2s});
 disp('P sparse')
 showres
 
 % transposed form - S1 first
 disp('transposed form')
 
 [S,P,A,Ix]=mergestates({S1,P1ft,A1,Ix1},{S2,P2'},options);
 disp('P full')
 P=P'; showres
 
 [S,P,A,Ix]=mergestates({S1,P1ft,A1,Ix1},{S2,P2s'},options);
 disp('P sparse')
 P=P'; showres
  
 % regular form - S2 first
 disp('S2 first in order')
 disp('regular form')
 
 [S,P,A,Ix]=mergestates({S2,P2},{S1,P1f,A1,Ix1});
 disp('P full')
 showres
 %textable([NaN+zeros(3,4) S'],0)
 %textable([S(Ix,:) A P],[0 0 0 0 2 2 2 2 2 2])
 
 [S,P,A,Ix]=mergestates({S2,P2s},{S1,P1f,A1,Ix1});
 disp('P sparse')
 showres
 
 % transposed form - S2 first
 disp('transposed form')
 
 [S,P,A,Ix]=mergestates({S2,P2'},{S1,P1ft,A1,Ix1},options);
  disp('P full')
 P=P'; showres
 
 [S,P,A,Ix]=mergestates({S2,P2s'},{S1,P1ft,A1,Ix1},options);
 disp('P sparse')
 P=P'; showres
 
 return
 
 % prints in LaTex format
 textable(P1m{1},2)
 textable(P1m{2},2)
 textable(P2,2)
 
 
 
function showres
  for i=1:3
    fprintf('        ')
    for j=1:6, fprintf('   %6i',S(j,i)); end
    fprintf('\n')
  end
  for i=1:10
    fprintf('%2i %2i %2i %2i',S(Ix(i),:),A(i));
    for j=1:6, 
      if P(i,j)==0
        fprintf('   0     '); 
      else
        fprintf('   %6.4f',P(i,j));
      end
    end
    fprintf('\n')
  end
end

end