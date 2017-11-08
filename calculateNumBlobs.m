clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',15,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% set parameters
strains = {'N2','npr1','daf22','daf22_npr1'}; % {'N2','npr1','daf22','daf22_npr1'}
saveResults = false;

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    legendList{strainCtr} = strain;
    filenames = importdata(['datalist/' strain '_list.txt']);
    
    %% initialise
    numFiles = length(filenames);

    %% go through individual movies
    for fileCtr = 1:numFiles
        
        %% load data
        filename = filenames{fileCtr};
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        skelData = h5read(filename,'/skeleton');
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        numFrames = max(trajData.frame_number);
        
        %% obtain features, filtering out single worms
        multiWormLogInd = logical(~trajData.is_good_skel);
        uniqueBlob = NaN(1,numFrames);
        
        %% go through each frame, get the number of blobs for each frame
        for frameCtr = 1:max(trajData.frame_number)
            thisFrameLogInd = trajData.frame_number == trajData.frame_number(frameCtr);
            uniqueBlob(frameCtr) = numel(unique(trajData.worm_index_joined(thisFrameLogInd & multiWormLogInd)));
        end
    end
end
            