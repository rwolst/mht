% We test our SURE bandwidth finder gives a logical answer for a cosine
% weight function
fprintf('-------------------------------------------')
fprintf('\nTesting the SURE optimal bandwidth method\n')
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

% We should get an optimal bandwidth that makes sense
I = Model.calculatePdgm(X, size(X,2));

% Perform SURE
[m, fval] = SURE_opt(I, W);
fprintf('\n\nFor n = %g, we have optimal m = %g\n', n, m)

% Finally plot the periodogram with optimal m
F = Model.calculateWtdPdgm(I, W(m));

% Plot true and weighted on same graph
figure
title('White noise process, weighted vs true')
hold on
plot(1:n/2 + 1, 1, 'r')
plot(squeeze(F))
hold off

figure
title('White noise process, raw')
hold on
plot(squeeze(I))
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

% Perform SURE
[m, fval] = SURE_opt(I, W);
fprintf('\n\nFor n = %g, we have optimal m = %g\n', n, m)

% Now get the true best bandwidth by putting in true spectrum
f = @(v) 1./(1-2*A*cos(2*pi*v/n) + A^2);

% Finally plot the periodogram with optimal m's
F = Model.calculateWtdPdgm(I, W(m));

% Plot true and weighted on same graph
figure
title('AR(1) process, weighted vs true')
hold on
plot(1:n/2 + 1, f(1:n/2 + 1), 'r')
plot(squeeze(F))
hold off

% figure
% hold on
% plot(1:n/2 + 1, f(1:n/2 + 1), 'r')
% plot(squeeze(Data2.F))
% hold off

figure
title('AR(1) process, raw')
hold on
plot(squeeze(I))
hold off