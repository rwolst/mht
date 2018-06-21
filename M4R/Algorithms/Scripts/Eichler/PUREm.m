    % Get the m bandwidth values using PURE
    % First get the pilot periodogram
    I = Data.Periodogram();
    
    % Now get the weight sequence function
    Z = @(m) cos(pi * (-m/2:m/2)/m);
    W = @(m) Z(m)/sum(Z(m));

    % We now perform PURE on each chaneel of our periodogram and average to
    % get true m
    m_PURE = zeros(1, Data.r);
    
    for i = 1:Data.r
        % Perform SURE
        m_pilot = SURE_opt(I(i,i,:), W);

        % Get weighted periodogram with given m_pilot
        Data.m = m_pilot;
        Data.getWeights();
        Data.getWtdPdgm();
        F = Data.F;
        
        % Perform PURE to get our true m
        m_PURE(i) = PURE_opt(I(i,i,:), F(i,i,:), W);
        
        fprintf('\nFor channel %d we have optimal PURE m = %g\n\n', i, m_PURE(i))
    end
    
    % Average the m and get the weights
    Data.m = 2 * round(mean(m_PURE)/2);  % Periodogram smoothing bandwidth
    Data.getWeights;
    fprintf('\nUsing m value %g\n', Data.m)
    
    % Clear the dynamic variables so we re-calculate the weighted
    % periodogram
    Data.clearDynamicVariables();