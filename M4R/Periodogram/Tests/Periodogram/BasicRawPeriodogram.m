function [ I ] = BasicRawPeriodogram( X )
% The periodogram calculated in the most basic way (i.e. no FFT etc.)
    N = length(X(1,:));
    r = length(X(:,1));
    I = zeros(r,r,N);
    
    % Slow FT the data X
    for j = 1:N
        % j is used to get our frequency components
        frequency = 2*pi*(j - 1)/N;
        
        % Now perform the FT for the given frequency
        total = 0;
        for k = 1:N
            total = total + X(:,k) * exp(- frequency * sqrt(-1) * k);
        end
        
        I(:,:,j) = (1/N)*(total*total');
    end

end

