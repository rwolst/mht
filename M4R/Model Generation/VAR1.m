function X = VAR1( A,n, varargin )
%This function creates a VAR(p) process for 1...n observations, with
%co-efficient matrix A where A(:,:,p) is the (t-1)th co-efficient and 
%A(:,:,1) is the (p-1)th. varargin is the mu and sigma of the multivariate
%normal distribution if it is specified.
%   We first take X(:,1)= zeros and do an intial burn of 1000 iterations to
%   remove system transients. We then record n observations and return them
%   in the matrix X.

    nVarargs = length(varargin);
    
    %So I don't have to re-type
    %A = [0.2 0 0.3 0 0.3; 0.3 -0.2 0 0 0; 0.2 0 0.3 0 0;0.2 0.3 0 0.3 0;0.2 0 0.2 0.2 0.2];
    Burn = 1000;
    p = length(A(1,1,:));
    r = length(A(1,:,1));   
    Y = zeros(r,p);
    e = zeros(r,1);
    X = zeros(r,n);
    
    for i = 1:Burn
        
        %Generate normal distributed error vector
        if nVarargs == 0
            e = mvnrnd(zeros(r,1),eye(r))';
        else
            e = mvnrnd(varargin{1},varargin{2})';
        end
        
        %The VAR1 equation
        
        if p>1
            Temp = 0;
            for k = 1:p
                Temp = Temp + A(:,:,k)*Y(:,k);
            end
            
            for j = 2:p
                Y(:,j-1) = Y(:,j);
            end       
            Y(:,p) = Temp + e;
        else
            Y(:,1) = A(:,:,1)*Y(:,1) + e; 
        end
        

    end
    
    
    for i = 1:p
        X(:,i) = Y(:,i); 
    end
    
    
    for i = p+1:n
       
        if nVarargs == 0
            e = mvnrnd(zeros(r,1),eye(r))';
        else
            e = mvnrnd(varargin{1},varargin{2})';
        end
        
        if p>1
            Temp = 0;
            for k = i-p:i-1
                Temp = Temp + A(:,:,k-(i-p)+1)*X(:,k);
            end
                 
            X(:,i) = Temp + e;
        else
            X(:,i) = A(:,:,1)*X(:,i-1) + e; 
        end
        
    end

end

