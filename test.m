function [x1,x2,x3,x4,n1,n2,n3,n4]=test(A,b,x0,wb)

x=A\b;
[x1,n1,t1]=jaccbi(A,b,x0);
[x2,n2,t2]=gauss(A,b,x0);
[x3,n3,t3]=SOR(A,b,x0,wb);
[x4,n4,t4]=CG(A,b,x0,100);
iterNum_UP=500;
iterNum_Down=1;
iterNums=iterNum_Down:iterNum_UP;
Xs_j = [];
Xs_G = [];
Xs_S = [];
Xc_C = [];
Ts_j = [];
Ts_G = [];
Ts_S = [];
Tc_C = [];
for iterNum = iterNums
[x_j,out_j,t1] = jaccbi(A,b,x0,'',iterNum);
Ts_j(end+1) = out_j;
Xs_j(:,end+1) = x_j;
Time_j(iterNum) = t1;
[x_G,out_G,t2] = gauss(A,b,x0,'',iterNum);
Ts_G(end+1) = out_G;
Xs_G(:,end+1) = x_G;
Time_G(iterNum) = t2;
[x_S,out_S,t3] = SOR(A,b,x0,wb,'',iterNum);
Ts_S(end+1) = out_S;
Xs_S(:,end+1) = x_S;
Time_S(iterNum) = t3;
[x_C,out_C,t4] = CG(A,b,x0,iterNum);
Tc_C(end+1) = out_C;
Xc_C(:,end+1) = x_C;
Time_C(iterNum) = t4;
end
for k = 1:iterNum_UP - iterNum_Down + 1
error_j(k)= norm((Xs_j(:,k) - x), 2);
error_G(k)= norm((Xs_G(:,k) - x), 2);
error_S(k)= norm((Xs_S(:,k) - x), 2);
error_C(k)= norm((Xc_C(:,k) - x), 2);
end
error_m = norm(A\b - x);
plot(iterNums,error_j,'b--o',iterNums,error_G,'--c*',iterNums,error_S,'r',iterNums,error_C,'--g+',iterNums,error_m,'-mo');
legend('jaccbi迭代精度变化','Gauss迭代精度变化','SOR迭代精度变化','CG迭代精度变化','matlab内置函数精度');
title('不同迭代方式的收敛曲线');
end
