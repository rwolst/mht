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

fprintf('\nHow many statsitcs do you want to use for the normality test (200)?\n')
total_statistics = input('');

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
percent = 0.2;
deltat  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



T_matsuda = zeros(total_statistics, 1);
T_eichler = zeros(total_statistics, 1);

% We now find total_statistics for the edge edge for each algorithm
for t = 1:total_statistics
    disp(t)
    
    % 1: Matsuda
    % 2: Eichler
    
    %Create new data
    Data = Model;
    Data.n = n;
    Data.generateVAR1(A);
    for stat = 1:2
        % Clear any dynamic variables as they need to be recalculated
        Data.clearDynamicVariables()
        
        %Define weight sequence making sure we use the correct C and D
        if stat == 1
            Data.m = m0_Matsuda;
            Data.getWeights();
        elseif stat == 2
            Data.m = m0_Eichler;
            Data.w = quadratic_weights(Data.m);
        end
    
        %We store the Weighted Periodogram now so that it doesn't have to be
        %re-calculated
        if stat == 1
            Data.getWtdPdgm();
        elseif stat == 2
            % Calculate the tapered direct spectrum
            [Sd, C_h] = direct_spectrum(Data.X, percent, deltat);

            % Now the weighted direct spectrum
            Data.F = Data.calculateWtdPdgm(Sd, Data.w);
            
            QQ = sumcoh_ab(Data.F);
        end
    
        % The test matrix for the fully saturated model is simply the weighted
        % periodogram so we store this now
        Data.TestMatrix = Data.F;
    
        %We also store the inverse of the weighted periodogram for each
        %frequency value as otherwise this will be done multiple times in
        %TestMatrix.m
        for i = 1:Data.n/2
            Data.InverseTestMatrix(:,:,i) = inv(Data.TestMatrix(:,:,i));
        end
        
        if stat == 1  
            % Get Matsuda statistic
            [Z, Divergence] = TestStat(Data, E1, E2);
            T_matsuda(t) = real(Z);
        elseif stat == 2               
            T_eichler(t) = EichlerStat(Data.n, m0_Eichler, edge, C_h, QQ); 
        end
    end
end


[~,p_mat] = kstest(T_matsuda);
[~,p] = kstest(T_eichler);
fprintf(sprintf('\n\nKS normality test p-values    : Matsuda %g \t Eichler %g', p_mat, p))

[~,p_mat] = jbtest(T_matsuda);
[~,p] = jbtest(T_eichler);
fprintf(sprintf('\nJB normality test p-values    : Matsuda %g \t Eichler %g', p_mat, p))

[~,p_mat] = lillietest(T_matsuda);
[~,p] = lillietest(T_eichler);
fprintf(sprintf('\nLillie normality test p-values: Matsuda %g \t Eichler %g\n', p_mat, p))

scatter(T_eichler, T_matsuda)

        
    