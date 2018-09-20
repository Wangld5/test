function [Rank,index]=PageRank(A)
N=size(A,1);
q=0.85;
d=sum(L,2);
D=diag(d);
%����ת�ƾ���
M=L'*inv(D);
e=ones(N,1);
a=(d==0);
%������ʾ�����������,����������ҳ���������
S=M+e*a'/N;
%�������յĸ���ת�ƾ���
G=q*S+(1-q)*e*e'/N;

%��ʼ��H0,H1
H0=zeros(N,1);
H1=ones(N,1);
while norm(H1-H0)>=1e-6
    H0=H1;
    H1=G*H0;
end
[Rank, index]=sort(H1,'descend');
end