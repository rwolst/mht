function [typeI, typeII] = LargePVAR(varargin)
    % This analyses our algorithm for large p, if the run flag is true,
    % data is generated and if the analyse flag is true, data is analysed for
    % type I and II errors
    
    % If you want to plot this function to compare typeI and typeII errors,
    % use the small script where the dividing numbers are the true number
    % of edges and non-edges for n= 2048, m=512, p = 150
%     for i = 1:15
%         [typeI(i), typeII(i)] = LargePVAR(10^(-i));
%     end
%     plot(typeI/(11175-3975), typeII/3975)

    run = false;
    analyse = true;
    adaptive = false;

    nVarargin = length(varargin);
    if nVarargin == 1
        alpha = varargin{1};
    else
        alpha = 0.05;
    end
    
    n = 2048;
    m = 512;
    p = 150;
    lambda = 0.0000001; %A parameter to be chosen which is an estimate for proportion of true null hypotheses   

    if run == true
        Data= Model;

        Data.n = n;
        Data.m = m;
        Data.getWeights;


        start = tic;
        p
        M = SparseStatVar1(p);
        Data.generateVAR1(M);

        [T,E,Divergence] = MultHypTest(Data, m , alpha, 'parallel');
        %E = MatsudaTest(Data,alpha);
        elapsed(p) = toc(start);

        %We now save the model parameters
        save(strcat('LargePVAR_MandT_p_', int2str(p), '_n_', int2str(n), '_m_', int2str(m)), 'M', 'T', 'Data');
    end

    if analyse == true
        load(strcat('LargePVAR_MandT_p_', int2str(p), '_n_', int2str(n), '_m_', int2str(m)));
        
        Data = Model;
        Data.r = p;
        if adaptive == true
            E = AdaptiveTest(Data, T, alpha, lambda);
        else
            E = NonAdaptiveTest(Data, T, alpha);
        end


        % Now find type I and II error
        %.................................%

        [typeI, typeII] = GetErrors(M, E, p);
    end
end

function E = NonAdaptiveTest(Data, T, alpha)
    P = T;
%     P(:,3) = [P(:,3)>=0].*(1 - normcdf(P(:,3),0,1)) + [P(:,3)<0].*normcdf(P(:,3),0,1);
    P = sortrows(P,3);

    exit = 0;
    l = 0;
    L = length(T(:,1));
    split = L + 1;
    
    while (exit == 0 && l < L)
        l = l + 1;
        CritLevel = norminv(1 - (alpha/(L-l+1)),0,1);
        %%CritLevel = norminv((1 - alpha)^(1/(L-l+1)),0,1);
        if P(L+1-l,3) < CritLevel
            exit = 1;
            split = l;
        end
    end
    
    E = ones(Data.r,Data.r);
    
    if split ~= L + 1
        for i = 1:(L-split+1)
            E(P(i,1),P(i,2)) = 0;
            E(P(i,2),P(i,1)) = 0;
        end
    end
end

function E = AdaptiveTest(Data, T, alpha, lambda)
    % Find adaptive critical levels
    P = T;
    P(:,3) = [P(:,3)>=0].*(1 - normcdf(P(:,3),0,1)) + [P(:,3)<0].*normcdf(P(:,3),0,1);
    P = sortrows(P,3);
    
    
    %Find estimate for proportion of true null hypotheses   
    R = sum(P(:,3)<=lambda);
    dim = Data.r * (Data.r - 1)/2;
    thEst = (dim - R + 1)/(dim*(1-lambda));
    
    m(1) = dim*thEst;
    for i = 2:dim+1
        m(i) = sum(P(:,3)>alpha/m(i-1));
    end
    
    k = 1;
    m(dim+2) = 0;
    for i = 2:dim+1
        if m(i) >= m(i+1) && m(i)<=m(1)
            k = i;
        end
    end
    
    split = 0;
    for i = 1:R
        if P(i,3) <= alpha/m(k)
            split = i;
        end
    end
    
    
    E = ones(Data.r,Data.r);
    
    if split ~= dim
        for i = 1:(dim-split)
            E(P(dim-i+1,1),P(dim-i+1,2)) = 0;
            E(P(dim-i+1,2),P(dim-i+1,1)) = 0;
        end
    end
end

function [typeI, typeII] = GetErrors(M, E, p)
    % We build a to be the true adjacancy matrix for VAR(1) process defined
    % by M
    A = zeros(p,p);
    for k = 1:p
        A(k,k) = 1;
    end

    miss = 0;
    for i = 1:p-1
        for j = (i+1):p
            if (abs(M(:,i)'*M(:,j)) > 0.0001) || (abs(M(i,j) > 0.0001)) || (abs(M(j,i) > 0.0001))
                A(i,j) = 1;
                A(j,i) = 1;
                miss = miss + 1;
            else
                A(i,j) = 0;
                A(j,i) = 0;
            end
        end
    end

    B = E-A;

    typeI = 0;
    typeII = 0;
    
    % Now find the errors, note that type I error is when an non-edge in
    % the true model is included in our model i.e. our test failed to
    % accept it as a non-edge. The highest number of type I errors possible
    % is p*(p-1)/2 - miss and highest number of typeII errors is miss.
    for i = 1:(p-1)
        for j = (i+1):p
            if B(i,j) == 1
                typeI = typeI+1; %If there is a 1 then we have included a wrong edge in E i.e. type I error
            end
            if B(i,j) == -1
                typeII = typeII+1; % If there is a -1 we didn't include an edge that should be i.e. type II error
            end
        end
    end    
end




