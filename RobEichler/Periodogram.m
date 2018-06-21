function Sp = Periodogram(X, deltat)% zero mean data!!
% References:
% Y. Matsuda (2006) "A test statistic for graphical modelling of
% multivariate time series", Biometrika 93, 399-409.
% R. Wolstenholme and A.T. Walden (2013) "A multiple hypothesis test
% approach to graphical modelling of multivariate time series.''
%

%   Returns the periodogram values for a r-vector-valued time series 
%   with n samples. Values are computed at Fourier frequencies 
%   0,1/(n deltat), 2/(n \deltat) up to 1/(2 deltat) 
%   ie at the +ve Fourier frequencies
%
%
% input: 
%
% X(r,n) a matrix with n columns of r-vector-valued zero meandata
% deltat: sampling interval of data
%
% output:
%
% Sp(r,r, (n/2)+1): periodogram spectral matrix at the +ve Fourier freqs
%
    n = length(X(1,:));
% we want to FT along rows
    W = transpose(fft(transpose(X)));
    r = length(W(:,1));

%
    Sp(1:r, 1:r, 1:(n/2)+1) = zeros; % ?????????
%   We now find the complex transpose of the column in W corresponding
%   to the frequencies lambda(j)
    
    for j = 1:(n/2)+1
        for i = 1:r
            for k = 1:r
                Sp(i,k,j) = W(i,j)*(W(k,j)')./n;
            end    
        end
    end
    Sp=Sp*deltat;
    
%   Periodograms are symmetric for real-valued time series so we only
%   consider the first n/2 values.

end

