load('C:\Users\rjw08\Documents\MATLAB\M4R\Plots and Variables\Variables\PhD\RealData\Mospatientspos.mat')



Epoch = 1;
EpochData(:,:) = Mospatientspos(Epoch,:,:);

Data = Model(EpochData')

Data.m = Data.n/16;

Data.getWeights()

[T,E,Divergence] = MultHypTest(Data);
