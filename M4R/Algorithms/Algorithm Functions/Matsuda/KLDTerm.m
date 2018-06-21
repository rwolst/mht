function  e = KLDTerm(X,Y)
%Function finds the Fourier frequency term in the Kullback-Liebler 
%divergence corresponding to the spectral density matrices X and Y

    n = length(X(:,1));
    temp = X/Y;
    e = trace(temp - eye(n)) - log(det(temp));
    

end

