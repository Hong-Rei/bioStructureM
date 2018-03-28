function [N,model,ref,MASS] = RMSD1(no_of_trj_mat)
NT=str2num(no_of_trj_mat);%these matrices are stored as 3N by m matrices
load TRJ1.mat
[N3,model]=size(TRJ1);
MASS=textread('MASS.txt');
N=N3/3;

load 50trj/avg_coord.mat

checkrmsd=10000;
s1=zeros(N3,1);
for nt = 1 : NT
NTstr=['load TRJ' num2str(nt) '.mat'];
eval(NTstr);
ss=['TRJnt=TRJ' num2str(nt) ';'];
eval(ss);
[N3,model]=size(TRJnt)
['In trajecotory'  num2str(nt) ]

for k=1:model;
     vector=reshape(TRJnt(:,k),3,N)';
     [R,T,eRMSD,oRMSD,fromXYZ1] = rmsdfitm(ref,vector,MASS);
     if eRMSD<checkrmsd
       s1=vector;
       checkrmsd=eRMSD;
       eRMSD
       tt=['save struct1.mat s1'];
       eval(tt);
     end
end
clear TRJnt;
uu=['clear TRJ' num2str(nt)];
eval(uu);
end



return


