count = 0;
testType = 'FDR';
for jeezy = 1:3
for p = 10:2:50
    count = count + 1;
    
    Data = Model;
    %Enter sample size and bandwidth
    Data.n = 2048;
    Data.m = 512;
    alpha = 0.05;

    'p'
    p
    
    dim{count} = p;
    
    M{count} = real(SparseStatVar1(p));

    Data.getWeights;
    Data.generateVAR1(M{count})

    time{count} = cputime;
    [T,E{count}, ~] = MultHypTest(Data,Data.m,alpha, testType);
    time{count} = cputime - time{count};
%     time2{count} = cputime;
%     [T2,E2{count}, Divergence2] = AdaptiveMHT(Data,0.1);
%     time2{count} = cputime - time2{count};
    
    %We can now test E against what it should be

    for k = 1:p
        A{count}(k,k) = 1;
    end

    miss{count} = 0;
    for i = 1:p-1
        for j = (i+1):p
            if (abs(M{count}(:,i)'*M{count}(:,j)) > 0.0001) || (abs(M{count}(i,j) > 0.0001)) || (abs(M{count}(j,i) > 0.0001))
                A{count}(i,j) = 1;
                A{count}(j,i) = 1;
                miss{count} = miss{count} + 1;
            else
                A{count}(i,j) = 0;
                A{count}(j,i) = 0;
            end
        end
    end

    'E-A'
    B{count} = E{count}-A{count}
%     B2{count} = E2{count}-A{count}

    typeI{count} = 0;
    typeII{count} = 0;
%     typeI2{count} = 0;
%     typeII2{count} = 0;
    for i = 1:(p-1)
        for j = (i+1):p
            if B{count}(i,j) == 1
                typeI{count} = typeI{count}+1; %If there is a I then we have included a wrong edge i.e. type I error
            end
            if B{count}(i,j) == -1
                typeII{count} = typeII{count}+1;
            end
%             if B2{count}(i,j) == 1
%                 typeI2{count} = typeI2{count}+1; %If there is a I then we have included a wrong edge i.e. type II error
%             end
%             if B2{count}(i,j) == -1
%                 typeII2{count} = typeII2{count}+1;
%             end

        end
    end
end
end
    
