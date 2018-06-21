function G = TestMatrix( E,y,single,varargin)
%This outputs matrix G(a,b) = y(a,b) where (a,b) element of edge set E and
%inv(G)(a,b) = 0 otherwise
%   E is a matrix s.t. E(a,b) = E(b,a) = 1 if there is a connection between
%   to vertices (a,b) and E(a,b) = 0 otherwise. This function sets G(a,b) =
%   y(a,b) subject to the constraint that an edge exists between (a,b)

%   We allow the ability to pass the inverse of matrix y in as an extra
%   parameter.
    
    k = length(E(:,1));
    limit = k*5;
    
    if (length(varargin) == 1)
        inverse_y = varargin{1};
    else
        inverse_y = inv(y);
    end

    
    %Define set F of (a,b)'s who aren't in the model.
%     M = 0;
%     for i = 1:(k-1)
%         for j = (i+1):k
%             if E(i,j) == 0
%                 M = M + 1;
%                 F(M,:) = [i, j];        
%             end
%         end
%     end
    [i,j] = find(E==0);
    F = [i,j];
    F = F(F(:,1)<F(:,2),:);
    M = length(F);
    %M is now the number of missing edges in our model.
           
    
    G = y;
    
    if single == 1
        %This differentiates whether the model E we are using is the full
        %model (M == 0) or the model with one edge missing (M == 1)
        if M ~= 0
            i = F(1,1);
            j = F(1,2);
            H = inverse_y;

            G(i,j) = y(i,j) + (H(i,j)/(H(i,i)*H(j,j) - (H(i,j)*H(j,i))));
            G(j,i) = y(j,i) + (H(j,i)/(H(i,i)*H(j,j) - (H(i,j)*H(j,i))));
        end
    else
        %Below is the recursion used in Matsuda (15)
        if M ~= 0
            Temp = zeros(k,k);

            for n = 1:limit
                l = mod(n,M);
                if l == 0
                    l = M;
                end
                
                H = inv(G);
                for i = 1:k
                    for j = 1:k
                        if i >= j
                            if [j, i] == F(l, :)                                
                                Temp(i,j) = G(i,j) + H(i,j)/(H(i,i)*H(j,j) - H(i,j)*H(j,i));
                            else
                                Temp(i,j) = G(i,j);
                            end
                        else
                            if [i, j] == F(l, :)
                                Temp(i,j) = G(i,j) + H(i,j)/(H(i,i)*H(j,j) - H(i,j)*H(j,i));
                            else
                                Temp(i,j) = G(i,j);
                            end
                        end
                    end
                end
                G = Temp;
            end    
        end
    end
end

