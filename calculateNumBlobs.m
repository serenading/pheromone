clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

%% set parameters
strains = {'N2','npr1','daf22','daf22_npr1'}; % {'N2','npr1','daf22','daf22_npr1'}
saveResults = true;

%% initialise
numBlobFig = figure; hold on

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    legendList{strainCtr} = strain;
    filenames = importdata(['datalist/' strain '_list.txt']);
    
    %% initialise
    numFiles = length(filenames);
    uniqueBlob.(strains{strainCtr}) = NaN(numFiles,90000); % 90000 frames is 1hr at 25fps

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
        
        %% go through each frame, get the number of blobs for each frame
        for frameCtr = 1:numFrames
            thisFrameLogInd = trajData.frame_number == trajData.frame_number(frameCtr);
            uniqueBlob.(strains{strainCtr})(fileCtr,frameCtr) = numel(unique(trajData.worm_index_joined(thisFrameLogInd & multiWormLogInd)));
        end
    end
    
    %% smooth data and average across movies
    uniqueBlob.(strains{strainCtr}) = smoothdata(uniqueBlob.(strains{strainCtr}),2,'movmedian',frameRate*3);
    uniqueBlobMean.(strains{strainCtr}) = nanmean(uniqueBlob.(strains{strainCtr}),1);
    uniqueBlobStd.(strains{strainCtr}) = nanstd(uniqueBlob.(strains{strainCtr}),0,1);
    
    %% plot for strain
    strainFig = figure; hold on
    for fileCtr = 1:numFiles
        plot(uniqueBlob.(strains{strainCtr})(fileCtr,:))
    end
    title(strains{strainCtr})
    
    %% plot for overall fig
    set(0,'CurrentFigure',numBlobFig)
    time_x = (1:90000)/frameRate/60; % time in incriment of 3 seconds on x-axis
    plot(time_x,uniqueBlobMean.(strains{strainCtr})(1:90000))
end

%% format plot and save
set(0,'CurrentFigure',numBlobFig)
legend(legendList)
xlabel('time (min)')
ylabel('number of clusters')
figurename = 'figures/numBlobs';
if saveResults
    exportfig(numBlobFig,[figurename '.eps'],exportOptions)
    system(['epstopdf ' figurename '.eps']);
    system(['rm ' figurename '.eps']);
    % save data too
    save('figures/numBlobs.mat','uniqueBlob','uniqueBlobMean','uniqueBlobStd');
end
 