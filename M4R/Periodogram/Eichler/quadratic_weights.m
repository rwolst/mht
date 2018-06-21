function [w] = quadratic_weights(m)
% These are the Bartlett-Priestley quadratic weights as in Eichler (07,08);
% if rescaled by dividing by 2M you get approx a sum of unity.
% we prefer to rescale exactly

    % First convert our m to the Matsuda M
    M = m/2;
    assert(ceil(M) == floor(M))

    k = -M:M;
    f = k./(2*M); % actually is (k/N)/(2M/N) ie f/B_N in Eichler

    % Bartlett/Priestley window on [-1/2,1/2] is the following function
    w(k + M + 1) = (3/2)*(1-(2*f).^2);

    % Normalise
    w = w/(2*M);
end

