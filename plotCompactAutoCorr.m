clear
close all


%% set parameters
strain = 'daf22_npr1'; % {'daf22_npr1','daf22','npr1','N2'}
minBlobSize = 150;
minTraj = 25;

%% load file list and randomly select a movie
filenames = importdata(['datalist/' strain '_list.txt']);
numFiles = length(filenames);
fileCtr = randi(numFiles,1);

%% load data
filename = filenames{fileCtr};
trajData = h5read(filename,'/trajectories_data');
blobFeats = h5read(filename,'/blob_features');
frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
compactness = blobFeats.compactness;

%% filter worms
% get all worm indices available
uniqueWorms = unique(trajData.worm_index_joined);
% exclude those below a certain blob size threshold
uniqueWorms = unique(trajData.worm_index_joined(blobFeats.area>=minBlobSize));
% exclude those with very short trajectories
for wormCtr = 1:length(uniqueWorms)
    if nnz(trajData.worm_index_joined == uniqueWorms(wormCtr))<minTraj
        uniqueWorms(wormCtr) = NaN;
    end
end
uniqueWorms = uniqueWorms(uniqueWorms>0);

%% randomly sample 20 worm trajectories for auto correlation plots
sampleWorms = datasample(uniqueWorms,20,'Replace',false);
samplePlot = figure;
% for each sample trajectory
for wormCtr = 1:length(sampleWorms)
    % get compactness values for the trajectory
    compactWorm = compactness(trajData.worm_index_joined==sampleWorms(wormCtr));
    % check that frames are monotonically increasing
    if ~all(diff(trajData.frame_number(trajData.worm_index_joined==sampleWorms(wormCtr)))>0)
        warning ('frames missing from this trajectory')
    end
    % plot autocorrelation of compactness for that trajectory
    set(0,'CurrentFigure',samplePlot)
    subplot(4,5,wormCtr) = acf(compactWorm,length(compactWorm)-1);
end