% This script uses a VAR process with a non-identity error covariance
% matrix. We then run the 

clear

%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%
k = 0; % Allows us to control inverse covariance matrix
series = 5;
series_plot = [1,2]; % The series we want to see the cross-spectrum for
n = 2048;
m = 256;
sigma = eye(series); % The covariance of the VAR errors
sigma(1,2)=k; sigma(2,1)=k;
counter = 0;
step = 0.01;


% Select one of the VAR models below by uncommenting

%A(:,:,4) = [0.6 0 0 0 0 ; 0 0.5 0 0.6 0; 0 0 0.8 0 0; 0 0 0 0.5 0; 0 0 -0.2 0 0.7];
%A(:,:,3) = [0 0.65 0 0 0; 0 -0.3 0 0 0; 0 0 -0.7 0 0; 0 0 0.9 0 0.4; 0 0 0 0 -0.5];
%A(:,:,2) = [0 0 0 0 0 ; 0 0 0 0 0 ; 0 0 0 0 -0.1; 0 0 0 0 0; 0 0 0 0 0];
%A(:,:,1) = [0 0 0 0 0; 0 0 -0.3 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
x = 0.0;
A(:,:,1) = [0.2 0 0 0 0.3; 
            0.3 -0.2 x 0 0; 
            0.2 x 0.3 0 0;
            0.2 0.3 0 0.3 0;
            0.2 0 0.2 0.5 0.2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = Model;
Data.n = n;
Data.m = m;

Data.getWeights;
X = VAR1(A,n,zeros(series,1),sigma);
Data.X = X;

[T,E, Divergence] = MultHypTest(Data, m, 0.05);


