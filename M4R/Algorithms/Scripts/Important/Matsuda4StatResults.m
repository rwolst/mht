%%%%%%%%%%%%%%%%%%%%%%%%Variable Definitions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1024;
m = 64;
for k = -m/2:m/2
    w(k + m/2 + 1) = cos(pi*k/m);
end
x = 0;
C = 0.617;
D = 0.446;
A(:,:,1) = 0.1*[2 0 -1 0 -5; 4 -2 0 2 0; -2 0 3 0 1; 3 1 0 3 0; 0 0 0 5 2];
TotalRepeats = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;

%%%%%%%%%%%%%%My Model%%%%%%%%%%%%%%%%
a(1,:) = [2,3];
a(2,:) = [2,5];
a(3,:) = [3,4];
a(4,:) = [3,5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%New%%%%%%%%%%%%%%%%%%%%
for i = 1:TotalRepeats
    count = count + 1
    %X = VAR1(A,n,[0 0 0 0],[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1]);
    X = VAR1(A,n);
    E = ones(5);
    
    for j = 1:4
        E1 = E;
        
        E1(a(j,1),a(j,2)) = 0;
        E1(a(j,2),a(j,1)) = 0;
        
        T(1,j,i) = SpecificModelTest(X, w, m, C, D, E, E1);
    end
    
    %Use max for backwards stepwise selection, min for forwards
    [~,I1] = max(T(1,:,i));
    
    E(a(I1,1),a(I1,2)) = 0;
    E(a(I1,2),a(I1,1)) = 0;
    
    for j = 1:4
        E1 = E;
        
        E1(a(j,1),a(j,2)) = 0;
        E1(a(j,2),a(j,1)) = 0;
        
        if (j==I1)
            T(2,j,i) = NaN;
        else
            T(2,j,i) = SpecificModelTest(X, w, m, C, D, E, E1);
        end
    end
    
    [~,I2] = max(T(2,:,i));
    
    E(a(I2,1),a(I2,2)) = 0;
    E(a(I2,2),a(I2,1)) = 0;
  
    for j = 1:4
        E1 = E;
        
        E1(a(j,1),a(j,2)) = 0;
        E1(a(j,2),a(j,1)) = 0;
        
        if (j==I1)||(j==I2) 
            T(3,j,i) = NaN;
        else
            T(3,j,i) = SpecificModelTest(X, w, m, C, D, E, E1);
        end
    end
    
    [~,I3] = max(T(3,:,i));
    
    E(a(I3,1),a(I3,2)) = 0;
    E(a(I3,2),a(I3,1)) = 0;
    
    for j = 1:4
        E1 = E;
        
        E1(a(j,1),a(j,2)) = 0;
        E1(a(j,2),a(j,1)) = 0;
        
        if (j==I1)||(j==I2)||(j==I3)
            T(4,j,i) = NaN;
        else
            T(4,j,i) = SpecificModelTest(X, w, m, C, D, E, E1);
        end
    end
end
