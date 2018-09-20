function comTime()
    testNum=[10,50,100,200];
    i=1;
    for n = 1:4
        N=testNum(n);
        wb=1.6;
        V=diag(rand(N,1));
        M=orth(rand(N));
        A=M*V*M';
        b=normrnd(0,1,N,1);
        A2=normrnd(0,1,N,N);
        x0=zeros(N,1);
        [x3,n3,T3]=SOR(A,b,x0,wb);
        [x1,n1,T1]=jaccbi(A,b,x0);
        [x2,n2,T2]=gauss(A,b,x0);
        [x4,n4,T4]=CG(A,b,x0,100);
        [x5, T5] = gaus(A, b);
        [x6, T6]=liezhu(A,b);
        test1(i)=T1;
        test2(i)=T2;
        test3(i)=T3;
        test4(i)=T4;
        test5(i)=T5;
        test6(i)=T6;
        i = i+1;
    end
    plot(testNum,test1,'r',testNum,test2,'g',testNum,test3,'b',testNum,test4,'c',testNum,test5,'k', testNum,test6,'y');
    legend('jaccbiʱ��仯','Gaussʱ��仯','��������SORʱ��仯','CGʱ��仯','��˹��ȥ��ʱ��仯','����Ԫ��ȥ��ʱ��仯');
    title('��ͬ��ʽ���ʱ����ά���仯');
end