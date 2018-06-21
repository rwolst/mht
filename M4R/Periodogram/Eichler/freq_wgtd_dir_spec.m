function S = freq_wgtd_dir_spec(Sd, w)
% input 
%
% p gives 100p% cosine taper in the direct spectrum estimate
% deltat: sampling interval of data


% S(1:r,1:r,1:nreq) is the input direct spectrum estimate using tapering
%                   nreq=(n/2)+1
% w(-M:M) is the weighting matrix

[r,r,nreq]=size(Sd);
n=2*(nreq-1);
% Spectral matrix r x r x number of frequencies
    S(1:r,1:r,1:nreq) = zeros;
% weighting is done from -M to M (length(w)= 2M+1)
    M = (length(w) - 1)/2;
% check sum of weights =1
%Norm=sum(w);
    
%   Note,  the number of frequencies  is equal to the number of 
%   observations and the spectrum is symmetric about 0 and modulo
%   indexing can be used to reflect about zero and Nyquist
    
    for j = 1:nreq
        for i = -M:M
            period = mod(j-1+i,n); 
            if (mod(period,n) > n/2 - 1)
                S(:,:,j) = S(:,:,j) + w(i + M + 1)*conj(Sd(:,:,n-period+1));
            else
                S(:,:,j) = S(:,:,j) + w(i + M + 1)*Sd(:,:,period+1);
            end
        end
    end
    
%    S = S/Norm;
    
 end

