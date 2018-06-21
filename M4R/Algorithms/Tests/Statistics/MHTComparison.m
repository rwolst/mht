% This tests the normality of the Matsuda and Eichler statistics where we
% choose the m for the model by sure
fprintf('-------------------------------------------')
fprintf('\nTesting the normality of our statistics\n')
fprintf('-------------------------------------------')

% Define weight generating function

Z         = @(m) cos(pi * (-m/2:m/2)/m);  % Matsuda
W_Matsuda = @(m) Z(m)/sum(Z(m));

W_Eichler = @(m) quadratic_weights(m);  % Eichler


fprintf('\nHow many observations n do you want for each realisation (1024)?\n')
n = input('');

fprintf('\nWhat dimension do you want the model to be (10)?\n')
p = input('');

fprintf('\nHow many realisations of the model (number of m we find) do you want\n to use for SURE/PURE (25)?\n')
k = input('');

fprintf('\nHow many Multiple Hypothesis tests do you want to use?\n')
total_tests = input('');

% We generate a random VAR model such that the 1 and 2 channel are
% conditionally independent
A = rand(p);

while true
    % Make it so dot product of columns 1 and 2 equals 0
    A(round(p/2)+1:p, 1)   = 0;
    A(1:round(p/2), 2)     = 0;
    A(1,2)                 = 0;
    A(2,1)                 = 0;
    assert(dot(A(:,1), A(:,2)) == 0)
    
    % Make it so that all eigenvalues are in the unit circle
    [V,D] = eig(A);
    ind = abs(D) > 1;
    if sum(ind) == 0
        break
    else
        % Invert the eigenvalues outside the unit circle
        D(ind) = 1./D(ind);
        A = V*D*inv(V);
    end
end

% Now generate a number of realisations from this process to find our
% SURE/PURE m's histogram


% The matrices we store our m values for each channel and realisation
M_Matsuda = zeros(p, k);  
M_Eichler = zeros(p, k);
for i = 1:k
    % Generate our realisation
    X = Model.VAR1(A, n);
    I = Model.calculatePdgm(X, n);
    
    % Perform SURE on each channel of our series
    for j = 1:p
        % Matsuda
        [m, fval]         = SURE_opt(I(j,j,:), W_Matsuda);
        M_Matsuda(j,i)    = m;
        
        % Eichler
        [m, fval]         = SURE_opt(I(j,j,:), W_Eichler);
        M_Eichler(j,i)    = m;
    end
end

% We now take the median m (50%), the m+10% (60%) and the m-10% (40%) to test our normality
% Matsuda
M_dist_Matsuda  = sort(mean(M_Matsuda, 1));
hist(M_dist_Matsuda)  % We plot the histogram of our M values
m0_Matsuda      = 2*round(M_dist_Matsuda(round(k/2))/2);
m_plus_Matsuda  = 2*round(M_dist_Matsuda(round(6*k/10))/2);
m_minus_Matsuda = 2*round(M_dist_Matsuda(round(4*k/10))/2);

% Eichler
M_dist_Eichler  = sort(mean(M_Eichler, 1));
hist(M_dist_Eichler)  % We plot the histogram of our M values
m0_Eichler      = 2*round(M_dist_Eichler(round(k/2))/2);
m_plus_Eichler  = 2*round(M_dist_Eichler(round(6*k/10))/2);
m_minus_Eichler = 2*round(M_dist_Eichler(round(4*k/10))/2);

% We now test the edge (1,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matsuda Stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edge = [1,2];
E1 = ones(p, p);
E2 = E1;
E2(edge(1), edge(2)) = 0;
E2(edge(2), edge(1)) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eichler Stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define percent (of cosine I think)
percent=0.2;
deltat=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Matsuda_I  = zeros(total_tests, 1);
Matsuda_II = zeros(total_tests, 1);
Eichler_I  = zeros(total_tests, 1);
Eichler_II = zeros(total_tests, 1);

% We now find total_statistics for the edge edge for each algorithm
for t = 1:total_tests
    disp(t)
    % 1: Matsuda
    % 2: Eichler
    
    % Create new data
    Data = Model;
    Data.n = n;
    Data.generateVAR1(A);
    alpha = 0.05;
    testType = 'FDR_dep';
    
    
    for stat = 1:2 
        % Clear any dynamic variables as they need to be recalculated
        Data.clearDynamicVariables()
        
        % Define weight sequence making sure we use the correct C and D
        if stat == 1
            statType = 'Matsuda';
            Data.m = m0_Matsuda;
        elseif stat == 2
            statType = 'Eichler';
            Data.m = m0_Eichler;
        end
        Data.getWeights;

        %Choose the algortihm to run
        [T_MultHyp, E_MultHyp, ~] = MultHypTest(Data, Data.m, alpha, testType, statType);

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

        % Now we can get type I and type II errors
        % Type I is when we reject the null hypothesis incorrectly i.e.
        % E_MultHyp(i,j) == 1 but E_True(i,j) == 0
        TypeI_errors  = sum(sum(E_MultHyp - E_True > 0))/2;  % Divide by 2 as graph is symmetric

        % Type II is when we accept the null hypothesis incorrectly i.e.
        % E_MultHyp(i,j) == 0 but E_True(i,j) == 1
        TypeII_errors = sum(sum(E_MultHyp - E_True < 0))/2;
        
        if stat == 1
            Matsuda_I(t)  = TypeI_errors;
            Matsuda_II(t) = TypeII_errors;
        elseif stat == 2
            Eichler_I(t)  = TypeI_errors;
            Eichler_II(t) = TypeII_errors;
        end
    end
end

total_edges = (p*(p-1)/2);

fprintf('\nAverage Type I  error rate (Matsuda vs Eichler):\t%g vs %g',   mean(Matsuda_I)/total_edges, mean(Eichler_I)/total_edges)
fprintf('\nAverage Type II error rate (Matsuda vs Eichler):\t%g vs %g\n', mean(Matsuda_II)/total_edges, mean(Eichler_II)/total_edges)


        
    