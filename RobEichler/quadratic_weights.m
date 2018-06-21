function [w C_w2 C_w4]=quadratic_weights(M)
%
% these are the Bartlett-Priestley quadratic weights as in Eichler (07,08);
% if rescaled by dividing by 2M you get approx a sum of unity.
% we prefer to rescale exactly
%

for k=-M:M
    f=k./(2.*M); % actually is (k/N)/(2M/N) ie f/B_N in Eichler
    
% Bartlett/Priestley window on [-1/2,1/2] is the following function
    w(k+ M+1)=(3/2).*(1-(2.*f).^2);
end
% sum(w) % ~= 2M
%w=w./sum(w);
w=w./(2*M);
%
% now we compute the constants 
% C_w2= \int k^2(u) du  and C_w4= \int k^4(u) du
%
% where k(u) is the inverse FT of the Bartlett/Priestley window
% given above. This should have k(0)=1 and according to Priestley p447
% it is given by
% k(u)=(3/(pi u)^2)*( [sin(pi u)/(pi u)]-cos(pi u) )
C_w2=quadv(@BartlettPriestley_sq,-20,20); % value  is 1.2
C_w4=quadv(@BartlettPriestley_p4,-10,10); % value is 0.8676