function [ C_w2, C_w4 ] = quadratic_constants( )
    % Now we compute the constants as in Eichler (07,08)
    % C_w2= \int k^2(u) du  and C_w4= \int k^4(u) du
    %
    % where k(u) is the inverse FT of the Bartlett/Priestley window
    % This should have k(0)=1 and according to Priestley p447
    % it is given by
    % k(u)=(3/(pi u)^2)*( [sin(pi u)/(pi u)]-cos(pi u) )
    C_w2=quadv(@BartlettPriestley_sq,-20,20); % value  is 1.2
    C_w4=quadv(@BartlettPriestley_p4,-10,10); % value is 0.8676
end

