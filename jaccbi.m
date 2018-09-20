function [x,out,t_j] = jaccbi(A,b,x0,eps, iterNum)
%输入参数eps表示精度，iterNum表示迭代步数
%当输入参数为五个时，eps失效，以迭代步数作为终止条件
n = length(b);
tic
D = diag(diag(A));
invD = diag(1./diag(A));
J = invD * (D - A);
if max(abs(eig(J))) >= 1
    error('Jacobi Algorithm Cannot Convergence！')
end
invDb = invD * b;
f = @(x)(J*x + invDb);
x = x0;
if nargin == 5
    for k = 1:iterNum-1
        x = f(x);
    end
    x_0 = x;
    x = f(x);
    out = norm(x-x_0,2);
else
    if nargin == 3
        eps = 1.0e-6;%缺省精度
    end
    out = 0;
    while 1
        x_0 = x;
        x = f(x);
        out = out + 1 ;
        if norm( x - x_0 ,2 ) < eps
            break;
        end
    end
end
toc
t_j = toc;
end
