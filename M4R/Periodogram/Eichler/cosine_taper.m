function h=cosine_taper(N,p)
%
% calculates a cosine taper of length N of various percentages p*100
%
  
    taper(1:N) =ones;
    M = floor(p*N/2);
    if (M > 0)
        range = 1:M;
        taper(range) =  0.5 * (1 - cos(pi*range/(M + 1)));
        taper((N+1-M):N) = taper(M:-1:1);
    end
    C = sqrt(sum(taper.^2));
    h=taper/C;
  
