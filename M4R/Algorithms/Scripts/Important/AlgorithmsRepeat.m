
for i = 4:12
    n = 2^i;
    m = n/8;

    %Enter Model
    x = 0.0;
    A(:,:,1) = [0.2 0 0.3 0 0.3; 0.3 -0.2 x 0 0; 0.2 x 0.3 0 0;0.2 0.3 0 0.3 0;0.2 0 0.2 0.2 0.2];


    %Define weight sequence making sure we use the correct C and D
    for k = -m/2:m/2
        w(k + m/2 + 1) = cos(pi*k/m);
    end
    C = 0.617;
    D = 0.446;

    %Create data
    X = VAR1(A,n);

    %Choose the algortihm to run
    %E=MatsudaTest(X,w,m,C,D);
    [T,E] = MultHypTest(X,w,m,C,D);

    Z(:,i-3)=T(:,3);
end