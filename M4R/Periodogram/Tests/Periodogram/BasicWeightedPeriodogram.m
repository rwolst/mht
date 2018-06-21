function [ J ] = BasicWeightedPeriodogram( X, w )
% The weighted periodogram calculated in the most basic way i.e with for
% loops

    % First get the raw periodogram
    I         = BasicRawPeriodogram(X);
    [r, ~, n] = size(I);
    I_index   = 0:1:n-1;  % The frequency index of the values in I
    
    % Now smooth the periodogram using our weights
    % To do this, we first extend the periodogram between -n/2:n-1 noting
    % that it is only unique between 0:n-1 so we can apply a filter to this
    % extended object to efficiently calculate the weighted periodogram
    
    % Extend periodogram
    I_ext       = zeros(r, r, 3*n/2);
    I_ext_index = -n:1:n;
    I_ext       = I(:,:,mod(I_ext_index, n)  + 1);  % Interpret extended I from I by symmetry
    
    % For index = 0 ... n-1, smooth the periodogram using the weights
    J = zeros(r, r, n/2 + 1);  % We only need it upt to n/2 + 1 as above this it is periodic
    for i = 1 : n/2 + 1
        % Get the index for I_ext that we need to use
        minIndex = find(I_ext_index == I_index(i)) - n/2;
        maxIndex = find(I_ext_index == I_index(i)) + n/2;
        
        % Now smooth the correct values in I_ext with the weights
        J(:,:,i) = sum(bsxfun(@times, w, I_ext(:,:,minIndex : maxIndex)),3);
    end
end