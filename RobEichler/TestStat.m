function T = TestStat( X,E1,E2,F,m,C,D )
% References:
% Y. Matsuda (2006) "A test statistic for graphical modelling of
% multivariate time series", Biometrika 93, 399-409.
% R. Wolstenholme and A.T. Walden (2013) "A multiple hypothesis test
% approach to graphical modelling of multivariate time series.''
%

% This function returns the test statistic as a function of the 
% models defined by matrices E1 and E2. Here E1 and E2 are square 
% matrices with a 1 entry indicating an edge exists and a 0 entry 
% indicating an edge doesn't exist.
% Thestatistic is asymptotically normally distributed.

%
%
% inputs:
%
% X(r,n):   the data matrix, r-vectors, n of them  
% E1(r,r), E2(r,r): as above
% F(r,r,(n/2)+1): Spectral matrix, r x r x number of frequencies
% m:       number of weights is m+1 (often m/2 is called M in literature
%          but Matsuda uses m as here)
% C, D:    Matsuda coefficients which depend on function u where w_k=u(k/m)
%
% outputs:
%
% T: value of Matsuda's  statistic, eqn (11) of Matsuda(2006)
%
    n = length(X(1,:));
    r = length(X(:,1));
    M1 = 0;
    M2 = 0;    
% We count the number of missing edges in the models E1 and E2 by
% checking for zeros in their defining matrices. We then call these M1
% and M2
    for i = 1:r
        for j = 1:i
            if E1(i,j) == 0
                M1 = M1 + 1;
            end
            
            if E2(i,j) == 0
                M2 = M2 + 1;
            end
        end
    end
% We make a point of noting if only one edge has been removed from the
% saturated graph, as this allows us to make a shortcut in building our
% testmatrix,
    if (M1 == 0) && (M2 == 1)
        single = 1;
    else
        single = 0;
    end
% now compute Matsuda's formula for the statistic
% here KLDModel is what he calls eKL
    T = (KLDModel(X,E1,E2,F, single)-(C*(M2-M1)/m))*(m*n/(D*(M2-M1)))^(0.5);
end

