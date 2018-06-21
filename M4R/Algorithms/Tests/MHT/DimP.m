clear

%This script allows us to run either algorithm quickly
%Define Data as the 'Model' class, that our functions take as an input
Data = Model;
%Enter sample size and bandwidth and significance level
Data.n = 2048;
Data.m = 256;
alpha = 0.05;
p = 50;

testType = 'FDR_dep';  % FDR, FDR_dep or FWER
statType = 'Eichler';

%Enter Model
x = 0;
A = SparseStatVar1(p);

%Define weight sequence making sure we use the correct C and D
if strcmp(statType, 'Eichler')
    Data.w = quadratic_weights(Data.m);
else
    Data.getWeights();
end

%Create data
Data.generateVAR1(A);

%Choose the algortihm to run
[T_MultHyp, E_MultHyp, ~] = MultHypTest(Data, Data.m, alpha, testType, statType);
% [T_Matsuda, E_Matsuda] = MatsudaTest(Data, alpha);

% Now print the type 1 and type 2 errors
% First get true model edges
E_True = ones(p,p);
for i = 1:p-1
    for j = (i+1):p
        if (abs(A(:,i)'*A(:,j)) > 0.0001) || (abs(A(i,j) > 0.0001)) || (abs(A(j,i) > 0.0001))
            E_True(i,j) = 1;
            E_True(j,i) = 1;
        else
            E_True(i,j) = 0;
            E_True(j,i) = 0;
        end
    end
end

% Now we can get type I and type II errors and empirical FDR

% Type I is when we reject the null hypothesis incorrectly i.e.
% E_MultHyp(i,j) == 1 but E_True(i,j) == 0
TypeI_errors  = sum(sum(E_MultHyp - E_True > 0))/2;  % Divide by 2 as graph is symmetric

% FDR is false discoveries (TypeI errors) over total discoveries (edges
% found)
FDR = TypeI_errors/(sum(sum(E_MultHyp))/2);

% Type II is when we accept the null hypothesis incorrectly i.e.
% E_MultHyp(i,j) == 0 but E_True(i,j) == 1
TypeII_errors = sum(sum(E_MultHyp - E_True < 0))/2;

TypeI_error_rate  = TypeI_errors/(p*(p-1)/2);
TypeII_error_rate = TypeII_errors/(p*(p-1)/2);

fprintf('\nType I  error rate %g', TypeI_error_rate)
fprintf('\nType II error rate %g\n', TypeII_error_rate)
fprintf('\nFDR %g\n', FDR)


