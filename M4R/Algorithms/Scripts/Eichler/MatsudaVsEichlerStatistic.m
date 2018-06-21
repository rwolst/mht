% Test the normality of the statistics produced by the Matsuda and the
% Eichler method
edge = [2,3];
total_statistics = 100;
T_matsuda = zeros(total_statistics, 1);
T_eichler = zeros(total_statistics, 1);

%Define Data as the 'Model' class, that our functions take as an input
Data = Model;
%Enter sample size and bandwidth and significance level
Data.n = 2048;
alpha = 0.05;

%Enter Model
x = 0;
A(:,:,1) = [0.2  0   -0.1 0   -0.5; 
            0.4 -0.2 x    0.2 0; 
            -0.2 x   0.3  0   0.1;
            0.3  0.1 0    0.3 0;
            0    0   0    0.5 0.2];
        


% Set up the models we are comparing in the Matsuda case
[r, ~] = size(A);
E1 = ones(r, r);
E2 = E1;
E2(edge(1), edge(2)) = 0;
E2(edge(2), edge(1)) = 0;

% Set up parameters for Eichler
M = 30;
% here we find C_h by more accurate scheme
% than in direct_spectrum
%%%%%%%%%%%%%[C_h, h2, h4]= intcostap(percent);

%get weights
w            = quadratic_weights(M);
[C_w2, C_w4] = quadratic_constants();

% Define percent (of cosine I think)
percent=0.2;
deltat=1;

% Define the weight sequence functions
Z         = @(m) cos(pi * (-m/2:m/2)/m);  % Matsuda
W_Matsuda = @(m) Z(m)/sum(Z(m));
W_Eichler = @(m) quadratic_weights(m/2);  % Eichler

% We now find total_statistics for the edge edge for each algorithm
for t = 1:total_statistics
    disp(t)
    
    %Create new data
    Data.generateVAR1(A);
    
    % Using SURE, get the Matsuda and Eichler m values
    I = Data.Periodogram();

    % Perform SURE
    m_Matsuda_List = zeros(r,1);
    m_Eichler_List = zeros(r,1);
    for i = 1:r
        m_Matsuda_List(i) = SURE_opt(I(i,i,:), W_Matsuda);
        m_Eichler_List(i) = SURE_opt(I(i,i,:), W_Eichler);
        fprintf('\nMatsuda m: %d', m_Matsuda_List(i));
        fprintf('\nEichler m: %d\n', m_Eichler_List(i));
    end
    m_Matsuda = 2*round(mean(m_Matsuda_List)/2);
    m_Eichler = 2*round(mean(m_Eichler_List)/2);
    
    %Define weight sequence making sure we use the correct C and D
    Data.m = m_Matsuda;
    Data.getWeights;
    
    %We store the Weighted Periodogram now so that it doesn't have to be
    %re-calculated
    Data.getWtdPdgm;
    
    %The test matrix for the fully saturated model is simply the weighted
    %periodogram so we store this now
    Data.TestMatrix = Data.F;
    
    %We also store the inverse of the weighted periodogram for each
    %frequency value as otherwise this will be done multiple times in
    %TestMatrix.m
    for i = 1:Data.n/2
        Data.InverseTestMatrix(:,:,i) = inv(Data.TestMatrix(:,:,i));
    end
    
    % Get Matsuda statistic
    [Z, Divergence] = TestStat(Data, E1, E2);
    T_matsuda(t) = real(Z);
    
    
    % Now Eichler
    T_eichler(t) = EichlerStat(Data.X, m_Eichler/2, percent, deltat, edge);
    
    
end

[h,p] = kstest(T_matsuda)
[h,p] = kstest(T_eichler)

[h,p] = jbtest(T_matsuda)
[h,p] = jbtest(T_eichler)

[h,p] = lillietest(T_matsuda)
[h,p] = lillietest(T_eichler)