function [T,E, Divergence] = MultHypTest(Data, m, alpha, controlType, statType)
%This function uses Matsuda's test statistic to determine missing edges in
%a graphical model, but does not iterate the procedure i.e. if there are
%three edges that when removed independently produce a test statistic below
%the critical value, we remove these edges and keep the others.

    % controlType  :  FWER or FDR or FDR_dep (if tests are negatively correlated)
    % statType     :  Matsuda or Eichler
    
    Data.r = length(Data.X(:,1));  % Dimension of our graph
    
    
    % Inititalise divergence as it sometimes isn't used
    Divergence = [];
    
    % For now pre-define the Eichler variables
    percent = 0.2;
    deltat  = 1;
    
    if nargin <4
        alpha = 0.1;
    end
    % the significance level of our test

    E1 = ones(Data.r,Data.r);
    E2 = ones(Data.r,Data.r,Data.r*(Data.r-1)/2);
    T = zeros(Data.r*(Data.r-1)/2,3);

    % We store the Weighted Periodogram now so that it doesn't have to be
    % re-calculated
    tic
    if isempty(Data.F)
        if strcmp(statType, 'Eichler')
            % We calculate the tapered weighted periodogram
            % Calculate the tapered direct spectrum
            [Sd, C_h] = direct_spectrum(Data.X, percent, deltat);

            % Now the weighted direct spectrum
            Data.F = Data.calculateWtdPdgm(Sd, Data.w);
            
            % Now find the unstandardized Eichler statistic
            QQ = sumcoh_ab(Data.F);
        else
            Data.getWtdPdgm;
        end
        'Pdgm Caclculated'
    else
        'Data.F not [] so assuming Pdgm already calculated'
    end
    fprintf('\nTime taken to calculate weighted periodogram: ')
    toc
    
    % The test matrix for the fully saturated model is simply the weighted
    % periodogram so we store this now
    Data.TestMatrix = Data.F;
    
    % We also store the inverse of the weighted periodogram for each
    % frequency value as otherwise this will be done multiple times in
    % TestMatrix.m
    tic
    for i = 1:Data.n/2
        Data.InverseTestMatrix(:,:,i) = inv(Data.TestMatrix(:,:,i));
    end
    
    counter = 0;
    for i = 1:(Data.r-1)
        for j = (i+1):Data.r 
            counter = counter + 1;
            E2(:,:,counter) = E1;
            E2(i, j, counter) = 0;
            E2(j, i, counter) = 0;
        end
    end      
    fprintf('\nTime taken to calculate inverse weighted periodogram: ')
    toc
    
    % Now put the maximum EKL level into our Data object
    Data.MaxCritLevel = norminv(1 - (2*alpha/(Data.r*(Data.r-1))),0,1);
    
    % We make a copy of data that each worker can share to save memory
    tic
    DataW = WorkerObjWrapper(Data);
    fprintf('\nTime taken to create wrapper object: ')
    toc
    
    
    tic
    if strcmp(statType, 'Matsuda')
        % Calculate Matsuda statistics
        parfor i = 1:Data.r*(Data.r-1)/2
            [Z, Divergence(i,:)] = TestStat(DataW.Value, E1, E2(:,:,i));
            T(i,3) = real(Z);
        end
    elseif strcmp(statType, 'Eichler')
        % Calculate the Eichler statistics
        for i = 1:Data.r*(Data.r-1)/2
            [r,c] = find(E2(:,:,i) == 0, 1);  % Row and column where we have no edge
            if r > c
                edge = [c, r];
            else
                edge = [r, c];
            end
            Z = EichlerStat(Data.n, m, edge, C_h, QQ);
            T(i,3) = real(Z);
        end
    end
    fprintf('\nTime taken to calculate statistics: ')
    toc
    
    % Now perform the multiple hypothesis test
    tic
    % Add the entry row and column to T
    counter = 0;
    for i = 1:(Data.r-1)
        for j = (i+1):Data.r  
            counter = counter + 1;
            T(counter,1) = i;
            T(counter,2) = j;
        end
    end    
    
    % Now perform our procedure
    assert(~any(isnan(T(:,3))))
    P = sortrows(T,3);
    L = length(T(:,1));  % Number of tests
    if strcmp(controlType, 'FWER') == 1
        exit = 0;
        l = 0;
        split = L + 1;
        while (exit == 0 && l < L)
            l = l + 1;
            CritLevel = norminv(1 - (alpha/(L-l+1)),0,1);
            %%CritLevel = norminv((1 - alpha)^(1/(L-l+1)),0,1);
            if P(L+1-l,3) < CritLevel
                exit = 1;
                split = l;
            end
        end
    elseif strcmp(controlType, 'FDR') == 1
        % We want to find the largest k such that 1 - normcdf(P(k)) > (k/L) * alpha
        % I.e. p-value(k) > (k/L)*alpha
        isGreater = 1 - normcdf(P(:,3)) > ((1:L)' / L) * alpha;
        
        % We now want the largest index i in isGreater such that
        % isGreater(i) == 1
        k = find(isGreater == 1, 1, 'Last');
        
        % Finally get our split value by accounting for the case when all
        % edges are significant
        if isempty(k)
            split = L + 1;
        else
            split = L-k+1;
        end
    elseif strcmp(controlType, 'FDR_dep') == 1
        % We want to find the largest k such that normpdf(P(k)) > (k/(c(L)*L)) * alpha
        % Where c(L) = sum(1./(1:L))
        c_L = sum(1./(1:L));
        isGreater = 1 - normcdf(P(:,3)) > ((1:L)' / (c_L * L)) * alpha;
        
        % We now want the largest index i in isGreater such that
        % isGreater(i) == 1
        k = find(isGreater == 1, 1, 'Last');
        
        % Finally get our split value by accounting for the case when all
        % edges are significant
        if isempty(k)
            split = L + 1;
        else
            split = L-k+1;
        end
    end

    % Now create our final model
    E = ones(Data.r,Data.r);

    if split ~= L + 1
        for i = 1:(L-split+1)
            E(P(i,1),P(i,2)) = 0;
            E(P(i,2),P(i,1)) = 0;
        end
    end
    
    fprintf('\nTime taken to choose graphical model using stepdown procedure: ')
    toc
    
    
end



