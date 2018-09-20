function [x,steps,t_C] = CG(A,b,x0,iterNum, eps)
r0 = b - A*x0;
p0 = r0;
if nargin == 4
    eps = 1.0e-6;
end
steps = 0;
tic
for k=1:iterNum-1
    if abs(p0) < eps
        break;
    end
    steps = steps + 1;
    a0 = r0'*r0/(p0'*A*p0);
    x1 = x0 + a0*p0;

    r1 = r0 -a0*A*p0;

    b0 = r1'*r1/(r0'*r0);
    p1 = r1 + b0*p0;

    x0 = x1;
    r0 = r1;
    p0 = p1;
end
x = x0;
toc
t_C = toc;
end