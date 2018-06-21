clear

n = 2048;
m = 512;
alpha = 0.05;

Data= Model;

Data.n = n;
Data.m = m;
Data.getWeights;



for p = 10:10
    start = tic;
    p
    M = SparseStatVar1(p);
    Data.generateVAR1(M);
    Data.StopStatEarly = true;
    
    [T,E] = MultHypTest(Data, m , alpha);
    %E = MatsudaTest(Data,alpha);
    elapsed(p) = toc(start)
    
    %We must now delete the periodogram in Data, as the model has changed
    Data.clearDynamicVariables();
end

save('./MandT.mat')

