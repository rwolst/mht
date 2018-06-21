r = 10;  % Dimension
n = 2048;  % Time points
m = 512;  % Bandwidth of smoothing
one_line_pdgm = true;
one_line_filter = true;
use_GPU = false;
pdgm_method = 2;
filter_method = 2;
run_old_method = true;


% Generate some random rxn data 
if use_GPU
    X = rand(r,n,'gpuArray');
else  
   X = rand(r, n); 
end

% Generate normalised weights according to a cos window
w = cos(pi * (-m/2:m/2)/m);
w = w/sum(w);

tic;
% Generate non-smoothed Periodogram
FT = (n)^(-0.5)*(ctranspose(fft(ctranspose(X))));

if one_line_pdgm == false
    if use_GPU
        Pdgm = zeros(r, r, n/2 + 1, 'gpuArray');
    else
        Pdgm = zeros(r, r, n/2 + 1);
    end

    for j = 1:n/2 + 1
        Pdgm(:,:,j) = FT(:,j)*FT(:,j)';
    end
else
    if pdgm_method == 1
        Pdgm_cell = cellfun(@(x) x * x', mat2cell(FT(:, 1 : n/2 + 1), [r], ones(n/2 + 1, 1)), 'UniformOutput', false);
        Pdgm = reshape(cell2mat(Pdgm_cell),r,r,[]);
    elseif pdgm_method == 2   
        gFT = FT(:,1:n/2 + 1);
        %// Perform non-smoothed Periodogram, thus removing the first loop
        Pdgm = bsxfun(@times,permute(gFT,[1 3 2]),permute(conj(gFT),[3 1 2]));
    end
end

% Finally smooth with our weights
if use_GPU
    SmPdgm = zeros(r,r,n/2+1,'gpuArray');
else
    SmPdgm = zeros(r,r,n/2+1);
end


% Take advantage of the GPU filter function
% Create new Periodogram J with m/2 values wrapped around in front and 
% behind it (it seems like there is redundancy here)
if use_GPU
    WrapPdgm = zeros(r,r,n/2 + 1 + m, 'gpuArray');
else
    WrapPdgm = zeros(r,r,n/2 + 1 + m);
end

WrapPdgm(:,:,m/2+1:n/2+m/2+1) = Pdgm;
WrapPdgm(:,:,1:m/2) = flip(Pdgm(:,:,2:m/2+1),3);
WrapPdgm(:,:,n/2+m/2+2:end) = flip(Pdgm(:,:,n/2-m/2+1:end-1),3);

% Perform filtering
if one_line_filter == false
    for i = 1:r
        for j = 1:r
            temp = filter(w, [1], squeeze(WrapPdgm(i,j,:)));
            SmPdgm(i,j,:) = temp(m+1:end);
        end
    end
else
    if filter_method == 1
        temp = filter(w, 1, WrapPdgm, [], 3);
        SmPdgm = temp(:, :, m + 1 : end);        
    elseif filter_method == 2
        filt_data = filter(w,1,reshape(WrapPdgm,r*r,[]),[],2);
        SmPdgm = reshape(filt_data(:,m+1:end),r,r,[]);        
    end
end

fprintf('\nTime for the Accelerated Smoothing ')
toc;


if run_old_method
    % Now check F with true weighted Periodogram
    if use_GPU
        X = gather(X);
    end
    Data = Model(X);
    Data.m = m;
    Data.getWeights();

    tic;
    Data.getWtdPdgm();
    fprintf('\nTime for the Non-Accelerated Smoothing ')
    toc;

    fprintf('\nWe should see 1 if both periodograms are equal\n')
    disp(all(all(all((SmPdgm - Data.F) < 0.001))))
end