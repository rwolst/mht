function S = WtdPdgm(X, w, deltat)
% input 
%
% X(r,n) a matrix with n columns of r-vector-valued data
% deltat: sampling interval of data

%This function outputs a weighted Periodogram for Fourier frequency
%2*(pi)*j/n and sequence of weights w(k)
%   For a sequence of weights w(k) for k = 0 to m we construct a weighted 
%   Periodogram as detailed in Matsuda 2006 equation 5.
    n = length(X(1,:));
    r = length(X(:,1));
    M = (length(w) - 1)/2;
    Norm = 0;
% Spectral matrix r x r x number of frequencies
    S(1:r,1:r,1:(n/2)+1) = zeros;
% compute sum of weights for standardization
    for i = -M:M
        Norm = Norm + w(i + M + 1);
    end
% the periodogram computed next is r x r x (n/2)+1 ie only pos freqs    
    Sp = Periodogram(X, deltat);
    
%   Note,  the number of frequencies  is equal to the number of 
%   observations and the spectrum is symmetric about 0 and modulo
%   indexing can be used to reflect about zero and Nyquist
    
    for j = 1:(n/2)+1
        for i = -M:M
            period = mod(j-1+i,n); 
            if (mod(period,n) > n/2 - 1)
                S(:,:,j) = S(:,:,j) + w(i + M + 1)*Sp(:,:,n-period+1);
            else
                S(:,:,j) = S(:,:,j) + w(i + M + 1)*Sp(:,:,period+1);
            end
        end
    end
    
    S = S/Norm;
    
 end

