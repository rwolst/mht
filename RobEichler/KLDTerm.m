function  e = KLDTerm(T1, T2)
% References:
% Y. Matsuda (2006) "A test statistic for graphical modelling of
% multivariate time series", Biometrika 93, 399-409.
% R. Wolstenholme and A.T. Walden (2013) "A multiple hypothesis test
% approach to graphical modelling of multivariate time series.''
%

% Finds the Kullback-Leibler divergence term corresponding to the spectral 
% density matrices T1 and T2 *at particular Fourier frequency* 

% find dimension of square matrix
    n = length(T1(:,1));
% now contributing term to Kullback-Leibler sum
    e = trace(T1/T2)-n - log(det(T1/T2));
end

