clear

%This script allows us to run either algorithm quickly
%Define Data as the 'Model' class, that our functions take as an input
Data = Model;
%Enter sample size and bandwidth
Data.n = 2048;
Data.m = 128;

%Enter Model
x = 0.1;
A(:,:,1) = [0.2 0 0.3 0 0.3; 0.3 -0.2 x 0 0; 0.2 x 0.3 0 0;0.2 0.3 0 0.3 0;0.2 0 0.2 0.2 0.2];
%A = [0 1 0 0; 0 0 1 1; 0 0 0 0; 0 0 0 0];
%A(:,:,1) = 0.1*[2 0 -1 0 -5; 4 -2 x 2 0; -2 x 3 0 1; 3 1 0 3 0; 0 0 0 5 2];


%Define weight sequence making sure we use the correct C and D
Data.getWeights;

%Create data
Data.generateVAR1(A)

%Choose the algortihm to run
%E=MatsudaTest(Data,0.05);
%[T,E,Divergence] = MultHypTest(Data);
[T,E,Divergence] = AdaptiveMHT(Data);
