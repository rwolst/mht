clear

%This script allows us to run either algorithm quickly
%Define Data as the 'Model' class, that our functions take as an input
Data = Model;
%Enter sample size and bandwidth and significance level
Data.n = 2048;
Data.m = 256;
alpha = 0.05;

testType = 'FDR';  % FDR, FDR_dep or FWER

%Enter Model
x = 0;
A(:,:,1) = [0.2  0   -0.1 0   -0.5; 
            0.4 -0.2 x    0.2 0; 
            -0.2 x   0.3  0   0.1;
            0.3  0.1 0    0.3 0;
            0    0   0    0.5 0.2];

%Define weight sequence making sure we use the correct C and D
Data.getWeights;

%Create data
Data.generateVAR1(A);

%Choose the algortihm to run
[T_MultHyp, E_MultHyp, ~] = MultHypTest(Data, Data.m, alpha, testType);
% [T_Matsuda, E_Matsuda] = MatsudaTest(Data, alpha);
