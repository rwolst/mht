function T = TestMatrix(E, S, single )
% References:
% Y. Matsuda (2006) "A test statistic for graphical modelling of
% multivariate time series", Biometrika 93, 399-409.
% R. Wolstenholme and A.T. Walden (2013) "A multiple hypothesis test
% approach to graphical modelling of multivariate time series.''
%
%  This outputs matrix T(a,b) = S(a,b) where (a,b) element of edge set E and
%  [inv(T)]_(a,b) = 0 otherwise. This is the Wermuth-Scheidt (1977)
%  algorithm.
%
% inputs:
%
%   E(r,r): is a matrix s.t. E(a,b) = E(b,a) = 1 if there is a connection between
%   to vertices (a,b) and E(a,b) = 0 otherwise. 
%   S(r,r): spectral matrix
%   single: 1 if we only remove one edge to get our test matrix 
%
%  output:
%   T: as above
%
    
    r = length(E(:,1));
    limit = r*5;

    
%  Define set FF of (a,b)'s which aren't in the model.
    M = 0;
    for i = 1:(r-1)
        for j = (i+1):r
            if E(i,j) == 0
                M = M + 1;
                FF(M,:) = [i, j];        
            end
        end
    end
%  M is now the number of missing edges in our model.
           
    
    T = S;
    
    if single == 1
% this differentiates between whether we are using the saturated model
% M==0 or the model with one edge missing M==1
        if M ~= 0 
            i = FF(1,1);
            j = FF(1,2);
            H = inv(S);

            T(i,j) = S(i,j) + (H(i,j)/(H(i,i)*H(j,j) - (H(i,j)*H(j,i))));
            T(j,i) = S(j,i) + (H(j,i)/(H(i,i)*H(j,j) - (H(i,j)*H(j,i))));
        end
    else
%  Below is the recursion used in Matsuda, eqn(15)
        if M ~= 0
            Temp = zeros(r,r);

            for n = 1:limit
                l = mod(n,M);
                if l == 0
                    l = M;
                end
                
                H = inv(T);
                for i = 1:r
                    for j = 1:r
                        if i >= j
                            if [j, i] == FF(l, :)                                
                                Temp(i,j) = T(i,j) + H(i,j)/(H(i,i)*H(j,j) - H(i,j)*H(j,i));
                            else
                                Temp(i,j) = T(i,j);
                            end
                        else
                            if [i, j] == FF(l, :)                              
                                Temp(i,j) = T(i,j) + H(i,j)/(H(i,i)*H(j,j) - H(i,j)*H(j,i));
                            else
                                Temp(i,j) = T(i,j);
                            end
                        end
                    end
                end
                T = Temp;
            end    
        end
    end
end

