function [T, Divergence] = TestStat(Data,E1,E2)
%This returns the test statistic as a function of the models E1 and E2. The
%statistic is asymptotically normally distributed.

    M1 = 0;
    M2 = 0;
    
    %We count the number of missing edges in the models E1 and E2 by
    %checking for zeros in their defining matrices. We then call these M1
    %and M2
    for i = 1:Data.r
        for j = 1:i
            if E1(i,j) == 0
                M1 = M1 + 1;
            end
            
            if E2(i,j) == 0
                M2 = M2 + 1;
            end
        end
    end
    
    %We make a point of noting if only one edge has been removed from the
    %saturated graph, as this allows us to make a shortcut in building our
    %testmatrix,
    if (M1 == 0) && (M2 == 1)
        single = 1;
    else
        single = 0;
    end
    
    % Now put the maximum EKL level into our Data object
    Data.MaxEKLLevel = Data.MaxCritLevel/((Data.m*Data.n/(Data.D*(M2-M1)))^(0.5)) + (Data.C*(M2-M1)/Data.m);
    
    %When testing with the True Spectrum, we set C&D to 0&1
    %C=0;
    %D=1;
    [Z, Divergence] = KLDModel(Data,E1,E2, single);
    T = (Z -(Data.C*(M2-M1)/Data.m))*(Data.m*Data.n/(Data.D*(M2-M1)))^(0.5);


end

