%('C:\Users\rjw08\Dropbox\MATLAB\M4R\Plots and Variables\Variables\PhD\MHTStatsn4096m256CVLL100stats.mat')

l = length(T(1,1,:));
step = 0.0005;

for loop = 0.00:step:0.9
    %We concentrate more of the alpha values near to zero for the Matsuda
    alpha = loop^5;
    %alpha = loop;
    
    FWER = 0;
    Power = 0;
    
    %We have to change the below MHT function based on the model (Matsuda
    %one is OK)
    if length(T(:,1,1)) == 10
        [A,~,~] = MHTStatEvaluation(T, alpha);
    else
        [A,~,~] = MatsudaStatEvaluation(T,alpha);
    end


    for i = 1:l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Incorrect Power Definition
        if A(i,3) == 1
            Power = Power + 1;
        end
        
        if A(i,1)==1
            Power = Power + 1;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Correct Power Definition
%         if A(i,3) == 1 && A(i,1)==1
%             Power = Power + 1;
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (A(i,2)==1)
            FWER = FWER + 1;
        end
    end
    
    B(round(loop/step)+1)=FWER/l;
    
    %Incorrect
    %C(round(loop/step)+1)=Power/(2*l);
    %Correct
    C(round(loop/step)+1)=Power/(l);
end

for i = 1:((0.9/step)+1)
    D(i) = step*(i-1);
end