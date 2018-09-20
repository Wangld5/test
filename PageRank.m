function [Rank,index]=PageRank(A)
N=size(A,1);
q=0.85;
d=sum(L,2);
D=diag(d);
%构造转移矩阵
M=L'*inv(D);
e=ones(N,1);
a=(d==0);
%构造概率矩阵的随机修正,处理悬挂网页的随机修正
S=M+e*a'/N;
%构造最终的概率转移矩阵
G=q*S+(1-q)*e*e'/N;

%初始化H0,H1
H0=zeros(N,1);
H1=ones(N,1);
while norm(H1-H0)>=1e-6
    H0=H1;
    H1=G*H0;
end
[Rank, index]=sort(H1,'descend');
end