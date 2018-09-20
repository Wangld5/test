function [x5, T5] = gaus(A2, b)
    tic
    B=[A2 b];
    n=length(b);
    RA=rank(A2);
    RB=rank(B);
    zhica=RB-RA;
    if zhica>0
        return
    end
    if RA==RB
        if RA==n
            x5=zeros(n,1);
            for p=1:n-1
                for k=p+1:n
                    m = B(k,p)/B(p,p);
                    B(k,p:n+1)=B(k,p:n+1)-m*B(p,p:n+1);
                end
            end
            b = B(1:n,n+1);
            A2=B(1:n,1:n);
            x5(n)=b(n)/A2(n,n);
            for q=n-1:-1:1
                x5(q)=(b(q)-sum(A2(q,q+1:n)*x5(q+1:n)))/A2(q,q);
            end
        else
            disp('no answer');
        end
    end
    toc
    T5=toc;
end