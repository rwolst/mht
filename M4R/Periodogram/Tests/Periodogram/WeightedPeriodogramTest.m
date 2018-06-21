% We test our periodogram function works correctly for
% 1) Single channel time series (white noise, AR(1))
% 2) Multi-channel time series
% by checking it can accuractely estimate the true spectral density

% Single channel test
% White noise
n = 32;
m = 4;
X = normrnd(zeros(1,n), 1);

% The spectral density chould be constant across all frequencies i.e. the
% periodogram should be relatively constant across all frequencies
Data = Model(X);
Data.m = m;
Data.getWeights();
Data.getWtdPdgm();

J1 = Data.F;
% We need to make the weights go from 1:n+1 as opposed to 1:m+1, by
% filling in with zeros
w = zeros(1, 1, n+1);
w(1, 1, n/2 - m/2 + 1:1:n/2 + m/2 + 1) = Data.w;

J2 = BasicWeightedPeriodogram(Data.X, w);

% Now test I1 and I2 are equal
% Note I1 runs from frequencies 0:n/2 and I2 runs from 0:n-1
assert(all(abs(J1 - J2).^2 < 1e-5))

% Multi-channel test
% White noise
n = 32;
X = normrnd(zeros(2,n), 1);

% The spectral density chould be constant across all frequencies i.e. the
% periodogram should be relatively constant across all frequencies
Data = Model(X);
Data.m = m;
Data.getWeights();
Data.getWtdPdgm();

J1 = Data.F;
% We need to make the weights go from 1:n+1 as opposed to 1:m+1, by
% filling in with zeros
w = zeros(1, 1, n+1);
w(1, 1, n/2 - m/2 + 1:1:n/2 + m/2 + 1) = Data.w;

J2 = BasicWeightedPeriodogram(Data.X, w);

% Now test I1 and I2 are equal
% Note I1 runs from frequencies 0:n/2 and I2 runs from 0:n-1
assert(all(all(all(abs(J1 - J2).^2 < 1e-5))))

disp('Test Passed!')