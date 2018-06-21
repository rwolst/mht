% We test our periodogram function works correctly for
% 1) Single channel time series (white noise, AR(1))
% 2) Multi-channel time series
% by checking it can accuractely estimate the true spectral density

% Single channel test
% White noise
n = 32;
X = normrnd(zeros(1,n), 1);

% The spectral density chould be constant across all frequencies i.e. the
% periodogram should be relatively constant across all frequencies
Data = Model(X);
I1 = Data.Periodogram();
I2 = BasicRawPeriodogram(Data.X);

% Now test I1 and I2 are equal
% Note I1 runs from frequencies 0:n/2 and I2 runs from 0:n-1
assert(all(abs(I1 - I2(:,:,1:n/2+1)).^2 < 1e-5))

% AR(1)

% Multi-channel test
% White noise
n = 32;
X = normrnd(zeros(2,n), 1);

% The spectral density chould be constant across all frequencies i.e. the
% periodogram should be relatively constant across all frequencies
Data = Model(X);
I1 = Data.Periodogram();
I2 = BasicRawPeriodogram(Data.X);

% Now test I1 and I2 are equal
% Note I1 runs from frequencies 0:n/2 and I2 runs from 0:n-1
assert(all(all(all(abs(I1 - I2(:,:,1:n/2+1)).^2 < 1e-5))))

disp('Test Passed!')