% We test our PURE bandwidth finder gives a logical answer for a cosine
% weight function
fprintf('-------------------------------------------')
fprintf('\nTesting the PURE optimal bandwidth method\n')
fprintf('-------------------------------------------')
fprintf('\nMost important thing is that the weighted and true periodograms are similar')
fprintf('\n\nWould you like to test\n1) The Matsuda weight function\n2) The Eichler weight function\n')
x = input('');

% Define weight generating function
if x == 1
    Z = @(m) cos(pi * (-m/2:m/2)/m);  % Matsuda
    W = @(m) Z(m)/sum(Z(m));
elseif x == 2
    W = @(m) quadratic_weights(m/2);  % Eichler
end

% White noise
fprintf('\n----------------------')
fprintf('\n| White noise process |')
fprintf('\n----------------------')
fprintf('\nExpecting the smoothing to be large as the true spectrum is flat (although not always necessarilly the case)\n')
n = 4096;
X = normrnd(zeros(1,n), 1);

% True white noise spectrum
f = @(v) 1;
S = f(1:n/2 + 1);

% We should get an optimal bandwidth that makes sense
I = Model.calculatePdgm(X, n);

% Perform SURE
[m, fval] = SURE_opt(I, W);


% Get the pilot periodogram with optimal m
F1 = Model.calculateWtdPdgm(I, W(m));

% Compare with true
fprintf('\n\nFor SURE n = %g, we have optimal m = %g and true MSE = %g\n', n, m, sum((squeeze(S) - squeeze(F1)).^2))

% Perform PURE
[m, fval] = PURE_opt(I, F1, W);

% Compare with true
F = Model.calculateWtdPdgm(I, W(m));

fprintf('\n\nFor PURE n = %g, we have optimal m = %g and true MSE = %g\n', n, m, sum((squeeze(S) - squeeze(F)).^2))

% ... and plot
% Plot true and weighted on same graph
figure
title('White noise process, SURE (r) vs PURE (b) vs true (g)')
hold on
plot(squeeze(F1), 'r')
plot(squeeze(F),  'b')
plot(squeeze(S),  'g')
hold off

% AR(1)
fprintf('\n----------------')
fprintf('\n| AR(1) process |')
fprintf('\n----------------')
fprintf('\nExpecting the smoothing to result in a spectrum close to the true spectrum\n')
n = 4096;
A = 0.5;
X = Model.VAR1(A, n);
I = Model.calculatePdgm(X, n);

% True AR(1) spectrum
f = @(v) (1./(1-2*A*cos(2*pi*v/n) + A^2)).';
S = f(1:n/2 + 1);

% Perform SURE
[m, fval] = SURE_opt(I, W);

% Get the pilot periodogram with optimal m
F1 = Model.calculateWtdPdgm(I, W(m));

% Compare with true
fprintf('\n\nFor SURE n = %g, we have optimal m = %g and true MSE = %g\n', n, m, sum((squeeze(S) - squeeze(F1)).^2))

% Perform PURE
[m, fval] = PURE_opt(I, F1, W);

% Compare with true
F = Model.calculateWtdPdgm(I, W(m));

fprintf('\n\nFor PURE n = %g, we have optimal m = %g and true MSE = %g\n', n, m, sum((squeeze(S) - squeeze(F)).^2))

% ... and plot
% Plot true and weighted on same graph
figure
title('AR(1) process, SURE (r) vs PURE (b) vs true (g)')
hold on
plot(squeeze(F1), 'r')
plot(squeeze(F),  'b')
plot(squeeze(S),  'g')
hold off