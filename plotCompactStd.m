clear
close all


%% set parameters
strains = {'daf22_npr1','daf22','npr1','N2'}; % {'daf22_npr1','daf22','npr1','N2'}
minBlobSize = 150;
minTraj = 25;

compactnessStandDevFig = figure; hold on
Legend = cell(length(strains),1);

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    Legend{strainCtr} = strain;
    filenames = importdata(['datalist/' strain '_list.txt']);
    numFiles = length(filenames);
    
    %% go through individual movies
    standDevPool = cell(numFiles,1);
    for fileCtr = 1:numFiles
        %% load data
        filename = filenames{fileCtr};
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        compactness = blobFeats.compactness;
        % get all worm indices available
        uniqueWorms = unique(trajData.worm_index_joined);
        %% exclude worms below a certain blob size threshold
        excludedWorms = unique(trajData.worm_index_joined(blobFeats.area<minBlobSize)); 
        filteredWorms = setdiff(uniqueWorms,excludedWorms);
        uniqueWorms = filteredWorms;
        %% go through each worm
        standDev = NaN(length(uniqueWorms));
        for wormCtr = 1:length(uniqueWorms)
            % exclude worms below min traj length threshold
            if nnz(trajData.worm_index_joined == uniqueWorms(wormCtr))>=minTraj
            compactWorm = compactness(trajData.worm_index_joined == uniqueWorms(wormCtr));
            standDev(wormCtr) = std(compactWorm);
            end
        end
        standDev(isnan(standDev))=[];
        % pool cross different movies
        standDevPool{fileCtr} = standDev;
    end
    standDevPool = horzcat(standDevPool{:});
    set(0,'CurrentFigure',compactnessStandDevFig)
    histogram(standDevPool,'Normalization','pdf','DisplayStyle','stairs')
end
set(0,'CurrentFigure',compactnessStandDevFig)
legend(Legend)
xlabel('standard deviation of compactness')
ylabel('P')
set(compactnessStandDevFig,'PaperUnits','centimeters')
figurename = 'figures/compactnessStd';
savefig(compactnessStandDevFig,[figurename '.fig'])
load('exportOptions.mat')
exportfig(compactnessStandDevFig,[figurename '.eps'],exportOptions)