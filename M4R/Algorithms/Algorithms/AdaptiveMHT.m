function [T,E, Divergence] = AdaptiveMHT(Data, alpha)
%This function uses Matsuda's test statistic to determine missing edges in
%a graphical model, but does not iterate the procedure i.e. if there are
%three edges that when removed independently produce a test statistic below
%the critical value, we remove these edges and keep the others. On top of
%that it uses the Guo adaptive Holm method for FWER control
    
    Data.r = length(Data.X(:,1));
    %r is the dimension of our graph
    
    if nargin <2
        alpha = 0.1;
    end
    %the significance level of our test

    E1 = ones(Data.r,Data.r);
    T = zeros(Data.r*(Data.r-1)/2,3);

    %We store the Weighted Periodogram now so that it doesn't have to be
    %re-calculated
    if isempty(Data.F)
        Data.getWtdPdgm;
        'Pdgm Caclculated'
    else
        'Data.F not [] so assuming Pdgm already calculated'
    end
    
    counter = 0;
    
    %The test matrix for the fully saturated model is simply the weighted
    %periodogram so we store this now
    Data.TestMatrix = Data.F;
    for i = 1:(Data.r-1)
        for j = (i+1):Data.r
            counter = counter + 1;
            
            E2 = E1;
            E2(i, j) = 0;
            E2(j, i) = 0;
            
            T(counter, 1) = i;
            T(counter, 2) = j;
            [Z, Divergence(counter,:)] = TestStat(Data, E1, E2);
            T(counter, 3) = real(Z);
            
        end
    end
    

    %Convert stats to P-values
    P = T;
    P(:,3) = [P(:,3)>=0].*(1 - normcdf(P(:,3),0,1)) + [P(:,3)<0].*normcdf(P(:,3),0,1);
    P = sortrows(P,3);
    
    
    %Find estimate for proportion of true null hypotheses
    lambda = 0.2; %A parameter to be chosen
    
    R = sum(P(:,3)<=lambda);
    M = Data.r * (Data.r - 1)/2;
    thEst = (M - R + 1)/(M*(1-lambda));
    
    m(1) = M*thEst;
    for i = 2:M+1
        m(i) = sum(P(:,3)>alpha/m(i-1));
    end
    
    k = 1;
    m(M+2) = 0;
    for i = 2:M+1
        if m(i) >= m(i+1) && m(i)<=m(1)
            k = i;
        end
    end
    
    split = 0;
    for i = 1:R
        if P(i,3) <= alpha/m(k)
            split = i;
        end
    end
    
    
    E = ones(Data.r,Data.r);
    
    if split ~= M
        for i = 1:(M-split)
            E(P(M-i+1,1),P(M-i+1,2)) = 0;
            E(P(M-i+1,2),P(M-i+1,1)) = 0;
        end
    end
end

