load('C:\Users\rjw08\Documents\MATLAB\M4R\HighDimSpacep30n2048m256.mat')
%load('C:\Users\rjw08\Documents\MATLAB\M4R\HighDimSpacen2048m128.mat')

Dimension = cell2mat(dim);


%Miss is actually the number of edges in the model and is a bit of a
%misnomer
Edges = cell2mat(miss);

NullEdges = (Dimension.^2 - Dimension)/2 - Edges;

%The non adaptive ratios
TypeIRatio = cell2mat(typeI)./NullEdges;
TypeIIRatio = cell2mat(typeII)./Edges;

%The adaptive ratios
TypeIRatioAdaptive = cell2mat(typeI2)./NullEdges;
TypeIIRatioAdaptive = cell2mat(typeII2)./Edges;