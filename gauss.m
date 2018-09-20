function [x,out,t_G] = gauss(A,b,x0,eps,iterNum)

n = length(b);
tic
D = diag(diag(A));
invD = diag(1./diag(A));
J = invD * (D - A);
invDb = invD * b;
H = eye(n) - tril(A)^-1*A;
if max(abs(eig(H))) >= 1
    error('Gauss_Seidel Algorithm Cannot Convergence£¡')
end

    function x = f(x)
        for l = 1:n
            temp_x = J(l,:)*x + invDb(l);
            x(l) = temp_x;
        end
    end
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
t_G = toc;
end
