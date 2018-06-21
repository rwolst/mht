clear

p = 5;
M = SparseStatVar1(p);

Data= Model;
alpha = 0.05;


count = 0;

for n = 200:100:10000
    m = 2*round(n/32); %Need m to be even and want it to be closest to n/16   
    count = count + 1  
    S(count,1) = n;
    S(count,2) = m;

    Data.n = n;
    Data.m = m;
    Data.getWeights;

    Data.generateVAR1(M);

    start = tic;
    %[T,E] = MultHypTest(Data);
    E = MatsudaTest(Data,alpha);
    S(count,3) = toc(start)
end