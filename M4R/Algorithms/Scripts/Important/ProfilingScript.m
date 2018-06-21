% Run this script using the profiler

clear

n = 1024;
m = 128;
alpha = 0.05;

Data= Model;

Data.n = n;
Data.m = m;
Data.getWeights;



for p = 30:30
    start = tic;
    p
    M = SparseStatVar1(p);
    Data.generateVAR1(M);
    
    [T,E,Divergence] = MultHypTest(Data, m , alpha, 'sequential');
    %E = MatsudaTest(Data,alpha);
    elapsed(p) = toc(start)
    
    %We must now delete the periodogram in Data, as the model has changed
    Data.F = [];
end
