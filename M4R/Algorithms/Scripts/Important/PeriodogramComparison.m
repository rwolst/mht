clear

n = 1024;
m = 64;
sigma = 1;
counter = 0;
step = 0.01;


%A(:,:,4) = [0.6 0 0 0 0 ; 0 0.5 0 0.6 0; 0 0 0.8 0 0; 0 0 0 0.5 0; 0 0 -0.2 0 0.7];
%A(:,:,3) = [0 0.65 0 0 0; 0 -0.3 0 0 0; 0 0 -0.7 0 0; 0 0 0.9 0 0.4; 0 0 0 0 -0.5];
%A(:,:,2) = [0 0 0 0 0 ; 0 0 0 0 0 ; 0 0 0 0 -0.1; 0 0 0 0 0; 0 0 0 0 0];
%A(:,:,1) = [0 0 0 0 0; 0 0 -0.3 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];

x = 0.0;
A(:,:,1) = [0.2 0 0.3 0 0.3; 0.3 -0.2 x 0 0; 0.2 x 0.3 0 0;0.2 0.3 0 0.3 0;0.2 0 0.2 0.2 0.2];

X = VAR1(A,n);

p = length(A(1,1,:));
r = length(X(:,1));
S = zeros(r, r, 0.5/step + 1);

data = Model(X);
data.m = m;
data.getWeights();
data.getWtdPdgm();
I = data.F;

% %Work out the weightings for the periodogram
% for i = -m/2:m/2
%     w(i + m/2 + 1) = cos(i*pi/m);
% end
% 
% %Calculate periodogram and actual spectrum
% I = WtdPdgm(X,w);

 for f = 0:step:0.5
    counter = counter + 1;
    G = eye(r);
    
    for j = 1:p
        G = G - A(:,:,p-j+1)*exp(2i*pi*f*j);
    end
    
    S(:,:,counter) = sigma * inv(G) * inv(G');
    
 end
 
 for series1 = 1:5
    for series2 = 1:5

        %Prepare the periodogram plot with a and b
        a = zeros(n,1);
        b = zeros(n,1);

        for j = 1:n/2
            a(j) = (-0.5) + (j-1)/n;
            a(n/2 + j) = j/n;
            b(n/2 - j + 1) = I(series1,series2,j);
            b(n/2 + j) = I(series1,series2,j);
        end


        %Prepare the true spectrum plot with c and d
        c = zeros(counter,1);
        d = zeros(counter,1);



        for j = 1:counter
            c(j) = (-0.5) + (j-1)/(2*counter);
            d(j) = S(series1,series2,counter - j + 1);

            c(counter + j) = j/(2*counter);
            d(counter + j) = S(series1,series2,j);

        end

        figure;       
        plot(a,b,c,d)
        title(sprintf('Graph of Series %d and Series %d', series1, series2))
    end
 end


%Green true spectrum, blue periodogram