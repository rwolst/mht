function z = KLDModel(X, E1, E2, F, single )
% References:
% Y. Matsuda (2006) "A test statistic for graphical modelling of
% multivariate time series", Biometrika 93, 399-409.
% R. Wolstenholme and A.T. Walden (2013) "A multiple hypothesis test
% approach to graphical modelling of multivariate time series.''
%

%   Finds the Kullback-Leibler divergence between two graphical models and
%   their spectral density matrices at the Fourier frequencies.
%   At each Fourier frequency, the spectral density matrix of
%   the multivariate time series X is estimated using a weighted 
%   periodogram with weights w. 
%   The spectral density matrices corresponding to the graphical 
%   models E1 and E2 (which are matrices representing connections 
%   between (j,k) by a 1 or a 0) are then found. Finally the 
%   Kullback-Liebler divergence at the Fourier frequencies is found
%   and returned as output z.
%
% inputs:
%
% X(r,n):   the data matrix, r-vectors, n of them  
% E1(r,r), E2(r,r): as above
% F(r,r,(n/2)+1): Spectral matrix, r x r x number of frequencies
% single: 1 if we only remove one edge to get our test matrix (as in MHT)
%
% output:
%
% z: estimated Kullback-Leibler divergence eKL, Matsuda eqn (7)
%
    n = length(X(1,:));
    Sum = 0;
    
   
    for j = 1:(n/2)
        %For an AR(1) process with co-efficient matrix A, the sdf is given
        %by: F = inv((eye(5) - A*exp(-2i*pi*(j+1)/n))^2);
        
        
        T1 = TestMatrix(E1,F(:,:,j),single);
        T2 = TestMatrix(E2,F(:,:,j),single);
        
        %TestMatrix Check
        %inv(T2)
        %F(:,:,j) - T2
        
        Sum = Sum + KLDTerm(T1,T2);
    end
    
    z = (1/n)*Sum;          
end

