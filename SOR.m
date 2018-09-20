function [x,out,t_S] = SOR(A,b,x0,w,eps,iterNum)
n = length(b);
tic
D = diag(diag(A));
invD = diag(1./diag(A));
J = invD * (D - A);
w_ = 1 - w;
invDb = invD * b;
L = invD * (D - tril(A));
B_w = D*(eye(n) - w*L)/w;
H_w = eye(n) - B_w^-1*A;
if max(abs(eig(H_w))) >= 1
    error('SOR Algorithm Cannot Convergence!!')
end

    function x = f(x)
        for l = 1:n
            temp_x = J(l,:)*x + invDb(l);
            x(l) = w_*x(l) + w*temp_x;
        end
    end
x = x0;
if nargin == 6
    for k = 1:iterNum-1
        x = f(x);
    end
    x_0 = x;
    x = f(x);
    out = norm(x-x_0,2);
else
    if nargin == 4
        eps = 1.0e-6;
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
t_S = toc;
end