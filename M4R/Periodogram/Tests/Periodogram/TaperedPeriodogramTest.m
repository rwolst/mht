% Define X as white noise
n = 1024;
m = 32;
X = rand(10,1024);
percent = 0.2;
deltat  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Andrew's way
% Get data dimensions
[~, N] = size(X);

% Get the smoothing weights
w            = quadratic_weights(m);
[C_w2, C_w4] = quadratic_constants();

% Calculate the tapered direct spectrum
[Sd, C_h] = direct_spectrum(X, percent, deltat);

% Now the weighted direct spectrum
S  = freq_wgtd_dir_spec(Sd, w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our way
Data = Model(X);
Data.w = quadratic_weights(m);

% Calculate the tapered direct spectrum
[Sd, C_h] = direct_spectrum(Data.X, percent, deltat);

% Now the weighted direct spectrum
Data.F = Data.calculateWtdPdgm(Sd, Data.w);

% Perform check
assert(all(all(all(abs(Data.F - S).^2 < 1e-5))))

disp('Test Passed!')