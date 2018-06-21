function [Sd C_h] = direct_spectrum(X, p, deltat)% zero mean data!!

%   Returns the single-taper direct spectrum estimate for
%   values for a r-vector-valued time series 
%   with n samples. Values are computed at Fourier frequencies 
%   0,1/(n deltat), 2/(n \deltat) up to 1/(2 deltat) 
%   ie at the +ve Fourier frequencies
%
%   taper is a 100p% cosine taper

% input: 
%
% X(r,n) a matrix with r rows of n-length zero mean data
% p  defines taper as 100p% cosine
% deltat: sampling interval of data
%
% output:
%
% Sp(r,r, (n/2)+1): direct  spectral matrix at the +ve Fourier freqs
%
    [r, n]=size(X);
% compute cosine taper  (normalized to unity)  
    h=cosine_taper(n, p);
    C_h=n.*sum(h.^4);
% taper each series (row)
    for j=1:r
        X(j,1:n)=h(1:n).*X(j,1:n);
    end
% we want to FT along rows
    W = transpose(fft(transpose(X)));

%
    Sd(1:r, 1:r, 1:(n/2)+1) = zeros; % ?????????
%   We now find the complex transpose of the column in W corresponding
%   to the frequencies lambda(j)
    
    for j = 1:(n/2)+1
        for i = 1:r
            for k = 1:r
                Sd(i,k,j) = W(i,j)*(W(k,j)');
            end    
        end
    end
    Sd=Sd*deltat;
    
%   Direct spectrum is  symmetric for real-valued time series so we only
%   consider the positive values.

end

