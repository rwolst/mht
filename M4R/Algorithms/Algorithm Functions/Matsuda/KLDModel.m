function [z, Divergence] = KLDModel( Data, E1, E2 ,single)
%Finds the Kulback-Liebler divergence between two graphical models and
%their spectral density matrices at the Fourier frequencies
%   At each Fourier frequency, we estimate the spectral density matrix of
%   the multivariate time series X using a weighted periodogram and weights
%   w. We then
%   work out the spectral density matrices corresponding to the graphical 
%   models E1 and E2
%   (which are matrices representing connections between (a,b) by a 1 or a
%   0). Finally we work out the Kullback-Liebler divergence at the Fourier
%   frequencies and return this as our answer.

    Sum = 0;
    %F = WtdPdgm(X,w);
    
     %We can test that our test works by using the true spectrum instead of a weighted periodogram
     %F = TrueVARPeriodogram([0.2 0 0.3 0 0.3; 0.3 -0.2 0 0 0; 0.2 0 0.3 0 0;0.2 0.3 0 0.3 0;0.2 0 0.2 0.2 0.2],n);
    
    Divergence = zeros(1,Data.n/2);
    for j = Data.StatLoopIndex        
        %For an AR(1) process with co-efficient matrix A, the sdf is given
        %by: F = inv((eye(5) - A*exp(-2i*pi*(j+1)/n))^2);
        
        
        G1 = Data.TestMatrix(:,:,j);
        G2 = TestMatrix(E2,Data.TestMatrix(:,:,j),single,Data.InverseTestMatrix(:,:,j));
        
        %TestMatrix Check
        %inv(G2)
        %F(:,:,j) - G2
        
        Divergence(j) = KLDTerm(G1,G2);
        Sum = Sum + Divergence(j);
        
        if Data.StopStatEarly
            if Sum > Data.n * Data.MaxEKLLevel
                % We stop early as statistic is certain to be rejected
                break
            end
        end
    end
    
    z = (1/Data.n)*Sum;    
        
end

