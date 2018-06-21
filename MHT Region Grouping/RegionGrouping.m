p=15;
Edges = [1 2; 2 3; 3 4; 3 6; 4 6; 7 10; 9 10; 8 9; 8 12; 12 13; 14 15];


A = zeros(p,p);
for i = 1:p
    A(i,i) = 0.1;
end

for i = 1:p-1
    for j=i+1:p
        if ismember([i,j],Edges,'rows');
            A(i,j) = 0.5;
        end
    end
end

A = A';

addpath('M4R');

%This script allows us to run either algorithm quickly
%Define Data as the 'Model' class, that our functions take as an input
Data = Model;
%Enter sample size and bandwidth
Data.n = 2048;
Data.m = 128;


%Define weight sequence making sure we use the correct C and D
Data.getWeights;

%Create data
Data.generateVAR1(A)

%Choose the algortihm to run
%E=MatsudaTest(Data,0.05);
[T,E,Divergence] = MultHypTest(Data);

S = [1 1 1 1 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 1 1 1 0 1 1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1];
Y = S*Data.X;

Data2 = Model(Y);
Data2.m = 128;

Data2.getWeights;

[T2,E2,Divergence] = MultHypTest(Data2);

%This shows that with the independent group summing we still get the
%correct answers.

%If we now swap over one of the nodes into a wrong group, 7 to 6
S2 = [1 1 1 1 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 1 0 1 1 1 0 1 1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1];
Y2 = S2*Data.X;

Data3 = Model(Y2);
Data3.m = 128;

Data3.getWeights;

[T3,E3,Divergence] = MultHypTest(Data3);

%It works!!!!!!!!

