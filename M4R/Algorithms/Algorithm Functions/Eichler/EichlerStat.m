function [ T ] = EichlerStat( N, m, edge, C_h, QQ )
    % We use this function, given our variables, to return the statistic from
    % Eichler
    % N       : N is total time indices.
    % m       : The m/2:m/2 bandwidth (needs to be converted to Eichler)
    % percent : The percent of cosine taper we use (common value 0.2)
    % deltat  : (common value 1)
    % edge    : The edge to test for example edge = [1, 2]
    % C_h     : A constant returned from calculating the tapered periodogram
    % QQ      : Unstandardized Eichler stat
    
    % Convert bandwidth
    M = m/2;
    assert(ceil(M) == floor(M))
    
    % The smoothing weights can be hard coded but also calculated by
    % removing the comment
    C_w2 = 1.199973539737920;
    C_w4 = 0.867559112862386;
    % [C_w2, C_w4] = quadratic_constants();
    

    % We can now compute Eichler's statistic Q_T
    % Form standardized Q_T
    % Note B_T=2M/N so M_T=1/B_T=N/(2M)
    M_T = N./(2*M);
    T   = N.*QQ(edge(1),edge(2))-M_T.*(C_h.*C_w2);
    T   = T./(C_h.*sqrt(M_T.*2.*C_w4));
end

