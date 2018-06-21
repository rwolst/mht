addpath(genpath('.'))

for i = 8:8
    matlabpool('local',i)
    SpeedAnalysisDimension
    matlabpool close force local
end

quit()
