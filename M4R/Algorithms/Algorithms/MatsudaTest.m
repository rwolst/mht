function [T,E] = MatsudaTest(Data,alpha)
%This function returns the graphs in an array E at each step of Matsuda's
%algorithm.
%Note that while the code for this works, it is far from optimal as many
%variables get recalculated at different steps.
    
    if nargin < 2
        alpha = 0.05;
    end
    %the significance level of our test
    
    marker = 0;
    Data.getWtdPdgm;
    
    T = [1, 2, 1000];
    %Simply initialising T where the value 100 will almost certainly be
    %larger than the calculated test statistics.
    
    E(:,:,1) = ones(Data.r,Data.r);
    %We define E[1] as a completely saturated graph.
    
    
    k = 1;
    
    while (marker == 0 && k <= Data.r*(Data.r-1)/2)
        %We initially set our values for the test matrix to be compared
        %against
        if k == 1
            Data.TestMatrix = Data.F;
        elseif k == 2
            %In this case can set the single variable to 1
            for j = 1:(Data.n/2)
                Data.TestMatrix(:,:,j)  = TestMatrix(E(:,:,k),Data.F(:,:,j),1);
            end
        else
            %In this case we must iterate to find answer as single is 0
            for j = 1:(Data.n/2)
                Data.TestMatrix(:,:,j)  = TestMatrix(E(:,:,k),Data.F(:,:,j),0);
            end
        end
        CritLevel = norminv((1-alpha)^(1/(Data.r*(Data.r-1)/2 - k + 1)),0,1);
        T = [1, 2, 100];
        %Simply initialising T (the vector that stores the minimum 
        %statistic and edge where it occured) where the value 100 will 
        %almost certainly be
        %larger than the calculated test statistics.
        
        for i = 1:(Data.r-1)
            i;
            for j = (i+1):Data.r
                j;
                E(:,:,k+1) = E(:,:,k);
                
                %Calculate the test statistic setting non-zero edges to
                %zero.
                if E(i, j, k+1) ~= 0
                    E(i,j,k+1) = 0;
                    E(j,i,k+1) = 0;
                    
                    Stat = TestStat(Data, E(:,:,k), E(:,:,k+1));
                    
                    %Set the minimum test statistic to T and include the
                    %corresponding edge (i,j).
                    if Stat < T(1,3)
                        T(1,1) = i;
                        T(1,2) = j;
                        T(1,3) = Stat;
                    end   
                end
            end
        end
                
        if real(T(1,3)) < CritLevel
            E(:,:,k+1) = E(:,:,k);
            E(T(1,1),T(1,2),k+1) = 0;
            E(T(1,2),T(1,1),k+1) = 0;
        else
            E(:,:,k+1) = zeros(Data.r,Data.r);
            marker = 1;
        end
        
        k = k + 1;
    end
    
end

