function [ R ] = PURE_obj_f( I, F_pilot, W, m0 )
% The objective function for the PURE bandwidth selection method
% Input:
%   I        - Raw single channel periodogram
%   F_pilot  - A pilot weighted periodogram calculated from the SURE optimal bandwidth
%   m0       - half the bandwidth (so we can constrain ga to integers)
%   W        - weight sequence function

    % First get weight sequence based on half bandwidth m0
    m = 2*m0;
    w = W(m);
    
    % Now get the frequency smoothed periodogram
    Data = Model;
    Data.w = w;
    Data.Pdgm = I;
    Data.getWtdPdgm();
    F = Data.F;
    n = Data.n;
    
    % We can now calculate the PURE risk estimate
    % How different our smoothed periodogram is from the pilot.
    pdgm_diff         = sum((F_pilot(1 : n/2) - F(1 : n/2)).^2);  % Lee says to n/2 but why not n/2 + 1?
    
    % We now manipulate our periodogram smoothing to get the second term
    Data1 = Model;
    Data1.Pdgm = F_pilot.^2;
    Data1.w = w.^2;
    Data1.getWtdPdgm();
    
    % How smooth our new weighted periodogram is compared to the pilot?
    pdgm_smooth = sum(Data1.F(1 : n/2));
    
    % The estimated risk
    R = pdgm_diff + pdgm_smooth;
    
end

