function M = SparseStatVar1(p)
%This simulates a sparse stationary VARp(1) matrix

%We first simulate sparse matrix A according to rule i+j = 1 mod 5

    %The leading diagonal cannot be 0
    for k = 1:p
        A(k,k) = normrnd(0,1);
    end

    for i = 1:(p-1)
        for j = (i+1):p
            if mod((i+j),5) == 1
                A(i,j) = normrnd(0,1);
                A(j,i) = normrnd(0,1);
            else
                A(i,j) = 0;
                A(j,i) = 0;
            end
        end
    end

    %We now send A to M where M is necessarilly stationary
    [V,D1] = eig(A);
    D2 = zeros(p,p);
    
    for i = 1:p
        if abs(D1(i,i)) > 1
            D2(i,i) = 1/(D1(i,i));
        else
            D2(i,i) = D1(i,i);
        end
    end
    
    M = V*D2*inv(V);
    
    for i = 1:p
        for j = 1:p
            if abs(M(i,j)) < 10^(-10)
                M(i,j) = 0;
            end
        end
    end
end

