function [ R ] = SURE_obj_f( I, W, m0 )
% The objective function for the SURE bandwidth selection method
% Input:
%   I  - Raw single channel periodogram
%   m0 - half the bandwidth (so we can constrain ga to integers)
%   W  - weight sequence function

    % First get weight sequence based on half bandwidth m0
    m = 2*m0;
    w = W(m);
    
    % Now get the frequency smoothed periodogram
    F = Model.calculateWtdPdgm(I, w);
    n = 2*(size(I, 3)-1);

    % We can now calculate the SURE risk estimate
    % How different our smoothed periodogram is from the orginal.
    pdgm_diff         = sum((I(1 : n/2) - F(1 : n/2)).^2);  % Lee says to n/2 but why not n/2 + 1?
    % How much of the actual raw periodogram point we use in each smoothing i.e. the central weight. In a cross-validation case this is 0.
    pdgm_point_weight = ((1 - 2*w(m/2+1))/2) * sum(I(1 : n/2).^2);
    
    % We want to keep pdgm_diff as low as possible (indicating close estimate) while having the pdgm_point_weight as high as possible (indicating smoothness)
    R = pdgm_diff - pdgm_point_weight;
end


