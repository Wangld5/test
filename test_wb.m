function test_wb(A,b,x0)
    testNum=1:0.1:1.9;
    step = [];
    for test=testNum
        [x3,n3,T3]=SOR(A,b,x0,test);
        step(:,end+1)=n3;
    end
    plot(testNum,step,'r');
    title('不同松弛因子迭代步数');
end