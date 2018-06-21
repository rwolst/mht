function [ m, fval ] = PURE_opt( I, F_pilot, W )
% Finds optimal bandwidth m for a given periodogram and weight function
% Input:
%   I        - Raw single channel periodogram from 0:n/2
%   F_pilot  - A pilot weighted periodogram calculated from the SURE optimal bandwidth
%   W        - weight sequence function

    % Define objective function
    obj_f = @(m0) PURE_obj_f(I, F_pilot, W, m0);
    
    % Get size of periodogram
    [~, ~, temp] = size(I);
    n = 2*(temp - 1);
    
    % Use optimiser to minimise
    [m0, fval] = ga(obj_f, 1, [], [], [], [], [1], [n/2], [], [1]);
    m = 2*m0;
end

