classdef Model<handle
    %This is the class containing all the data needed to run our tests
    
    properties
        n  % Observation Length
        m  % Smoothing Length
        r  % Dimension
        X  % Data
        F  % Weighted Periodogram
        Pdgm  % Raw periodogram
        w  % Weights
        C  % Constant
        D  % Constant
        DFT_points; % Allows us to estimate our DFT for a given length
        TestMatrix  % Dynamical Storage of Model TestMatrix
        InverseTestMatrix  % Dynamical Storage of inverse of Model TestMatrix
        StopStatEarly  % This is a flag that will stop calculating a test statistic once it is past the critical level
        MaxCritLevel  % The maximum Crit level for the tests in terms of the test statistic (if we ever hit above this we stop our calculation)
        MaxEKLLevel  % The maximum EKL level for the tests in terms of the test statistic (if we ever hit above this we stop our calculation)
        StatLoopIndex  % If we are stopping stat early, it may be worth looping over frequencies in a different way to simply 1...n/2
        
    end
    
    methods    
        function obj = clearDynamicVariables(obj)
            obj.F = [];
            obj.TestMatrix = [];
            obj.InverseTestMatrix = [];
        end
        
        function obj = Model(X, DFT_points)
            % Construct class given input X
            if nargin == 1
                obj.X = X;
                obj.n = length(X(1,:));
                obj.r = length(X(:,1));
                obj.DFT_points = obj.n;
            elseif nargin == 2
                obj.X = X;
                obj.n = length(X(1,:));
                obj.r = length(X(:,1));
                obj.DFT_points = DFT_points;
            end
            
            obj.StopStatEarly = true;
            obj.MaxCritLevel = inf;
            obj.MaxEKLLevel = inf;
        end
        
        function obj = initialiseLoopIndex(obj)
            % The ordering in which we loop over our weighted periodogram
            % when calculating the eKL
            obj.StatLoopIndex = 1:obj.n/2;
        end
        
        function obj = getWeights(obj)
            if isempty(obj.m)
                sprintf('First Enter a Value for m')
            else
                obj.w = cos(pi * (-obj.m/2:obj.m/2)/obj.m);
                obj.w = obj.w/sum(obj.w);
                
                obj.C = 0.617;
                obj.D = 0.446;
                obj.initialiseLoopIndex();
            end
        end
        
        function obj = generateVAR1(obj, A, n)
            % n is an optional parameter but o/w we need obj.n
            if nargin < 2
                print ('Enter Value for A')
            elseif nargin == 2
                if isempty(obj.n)
                    sprintf('Need a Value for n')
                end
            elseif nargin == 3
                obj.n = n;
                obj.DFT_points = n;
            end
            obj.X = obj.VAR1(A,obj.n);
            obj.r = length(obj.X(:,1));
            
        end
        
        function Pdgm = Periodogram(obj)
            % Returns the Periodogram value for a specific frequency lambda(j) from the
            % matrix X
            %   This discrete Fourier transforms the matrix X containing r dimnsional
            %   vectors of n time series observations, into the matrix W. We then use
            %   the fact that the Periodogram I = (W)(W*) the complex transpose.
            %   In this function we assume the Fourier frequencies running from 0 to
            %   2*pi i.e. 2*pi*j/n in this range.
            % The output is in the range 1:obj.n/2 + 1 noting obj.n/2 + 2:obj.n + 1 = 1:obj.n/2
            if isempty(obj.X)
                sprintf('We must have data values in X')
            else
                if size(obj.DFT_points) == 0
                    % Must define the number of DFT points so set to n
                    obj.DFT_points = obj.n;
                end
                   
                Pdgm = obj.calculatePdgm(obj.X, obj.DFT_points);
            end
        end      
                                       
        function obj = getWtdPdgm(obj)
            % This function outputs a weighted Periodogram for Fourier frequency
            % 2*(pi)*j/n and sequence of weights w(k)
            %   For a sequence of weights w(k) for k = 0 to m we construct a weighted 
            %   Periodogram as detailed in Matsuda 2006 equation 5.
            
            % If periodogram not already stored, we create a new one
            if isempty(obj.Pdgm)
                obj.Pdgm = obj.Periodogram();             
            end
            [obj.r, ~, temp] = size(obj.Pdgm);
            obj.n = 2*(temp-1);

            % Now use our static method to get the weighted periodogram
            obj.F = obj.calculateWtdPdgm(obj.Pdgm, obj.w);
        end
    end
    
    methods(Static)
        function W = DFT(X, DFT_points)
            W = (DFT_points)^(-0.5)*(transpose(fft(transpose(X), DFT_points)));
        end
        
        function Pdgm = calculatePdgm(X, DFT_points)
            % Returns the Periodogram value for a specific frequency lambda(j) from the
            % matrix X
            %   This discrete Fourier transforms the matrix X containing r dimnsional
            %   vectors of n time series observations, into the matrix W. We then use
            %   the fact that the Periodogram I = (W)(W*) the complex transpose.
            %   In this function we assume the Fourier frequencies running from 0 to
            %   2*pi i.e. 2*pi*j/n in this range.
            % The output is in the range 1:obj.n/2 + 1 noting obj.n/2 + 2:obj.n + 1 = 1:obj.n/2             
            FT = Model.DFT(X, DFT_points);
            gFT = FT(:, 1 : DFT_points/2 + 1);
            Pdgm = bsxfun(@times,permute(gFT,[1 3 2]),permute(conj(gFT),[3 1 2]));  % This permutes array dimensions so we can bsxfun them efficiently
        end   
        
        function F = calculateWtdPdgm(Pdgm, w)
            % This function outputs a weighted Periodogram for Fourier frequency
            % 2*(pi)*j/n and sequence of weights w(k)
            %   For a sequence of weights w(k) for k = 0 to m we construct a weighted 
            %   Periodogram as detailed in Matsuda 2006 equation 5.
            
            temp       = size(Pdgm);
            r          = temp(1);
            DFT_points = 2*(temp(3) - 1);
            m          = length(w) - 1;
            
            % Take advantage of the Matlab filter function
            % Create new Periodogram with m/2 values wrapped around in front and 
            % behind it (it seems like there is redundancy here)
            WrapPdgm                                         = zeros(r, r, DFT_points/2 + 1 + m);
            WrapPdgm(:, :, m/2 + 1 : DFT_points/2 + m/2 + 1) = Pdgm;
            % Wrapping is equivalent to taking the complex conjugate of the
            % matrix i.e. W(f) -> W(-f). As out matrices are Hermitian,
            % this is equivalent to taking the transpose which can be done
            % with permute
            WrapPdgm(:, :, 1 : m/2)                          = permute(flip(Pdgm(:, :, 2 : m/2 + 1), 3), [2 1 3]);  
            WrapPdgm(:, :, DFT_points/2 + m/2 + 2 : end)     = permute(flip(Pdgm(:, :, DFT_points/2 - m/2 + 1 : end - 1), 3), [2 1 3]);

            % Perform filtering
            filt_data = filter(w, 1 ,reshape(WrapPdgm , r*r, []), [], 2);
            F = reshape(filt_data(:, m + 1 : end), r, r,[]);  
        end
        
        function X = VAR1(A, n)
            % A simple inclusion of our model generation function into this
            % class
            X = VAR1(A, n);
        end
    end
    
end

