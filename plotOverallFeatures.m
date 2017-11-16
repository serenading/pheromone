clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',15,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',15,...
    'LineWidth',3);

%% set parameters
strains = {'npr1','daf22_npr1','N2','daf22',}; % {'N2','npr1','daf22','daf22_npr1'}
numSampleSkel = 500; % number of skeletons (per file) to sample in order to determine overall skeleton lengths for normalisation
areaCutOff = 8;
perimeterCutOff = 3.5;

saveResults = true;

%% initialise
swLengthFig = figure; hold on
swWidthFig = figure; hold on
swPerimeterFig = figure; hold on
swAreaFig = figure; hold on
perimeterFig = figure; hold on
areaFig = figure; hold on
perimeterThresholdFig = figure;
areaThresholdFig = figure;

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    legendList{strainCtr} = strain;
    filenames = importdata(['datalist/' strain '_list.txt']);
    
    %% initialise
    numFiles = length(filenames);
    perimeter = cell(numFiles,1);
    area = cell(numFiles,1);
    perimeterThres.(strains{strainCtr}) = NaN(1,numFiles);
    areaThres.(strains{strainCtr}) = NaN(1,numFiles);
    swLengths.(strains{strainCtr}) =  NaN(numFiles,numSampleSkel); 
    swWidths.(strains{strainCtr}) = NaN(numFiles,numSampleSkel);
    swPerimeters.(strains{strainCtr}) = NaN(numFiles,numSampleSkel);
    swAreas.(strains{strainCtr}) = NaN(numFiles,numSampleSkel);

    %% go through individual movies
    for fileCtr = 1:numFiles
        
        %% load data
        filename = filenames{fileCtr};
        trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        skelData = h5read(filename,'/skeleton');
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        
        %% obtain features, filtering out single worms
        multiWormLogInd = logical(~trajData.is_good_skel);
        perimeter{fileCtr} = blobFeats.perimeter(multiWormLogInd);
        area{fileCtr} = blobFeats.area(multiWormLogInd);
        
        %% calculate worm skeleton length (as daf-22 containing animals appear smaller) to normalise features against later
        % load xy coordinates
        xcoords = squeeze(skelData(1,:,:));
        ycoords = squeeze(skelData(2,:,:));
        % filter for single worms
        singleWormLogInd = logical(trajData.is_good_skel);
        xcoords = xcoords(:,singleWormLogInd);
        ycoords = ycoords(:,singleWormLogInd);
        singleWormArea = blobFeats.area(singleWormLogInd);
        singleWormPerimeter = blobFeats.perimeter(singleWormLogInd);
        % extract single worm features and calculate skeleton length
        [~,sampleSkelIdx] = datasample(1:size(xcoords,2),numSampleSkel,'Replace',false); % sample 500 random single worm skeletons
        xcoords = xcoords(:,sampleSkelIdx);
        ycoords = ycoords(:,sampleSkelIdx); 
        singleWormArea = singleWormArea(sampleSkelIdx);
        singleWormPerimeter = singleWormPerimeter(sampleSkelIdx);
        for skelCtr = 1:numSampleSkel
            skel_xcoords = xcoords(:,skelCtr);
            skel_ycoords = ycoords(:,skelCtr);
            dx = skel_xcoords(2:end)-skel_xcoords(1:end-1);
            dy = skel_ycoords(2:end)-skel_ycoords(1:end-1);
            dz = sqrt(dx.^2 + dy.^2);
            swLengths.(strains{strainCtr})(fileCtr,skelCtr) = sum(dz);
            swAreas.(strains{strainCtr})(fileCtr,skelCtr) = singleWormArea(skelCtr);
            swPerimeters.(strains{strainCtr})(fileCtr,skelCtr) = singleWormPerimeter(skelCtr);
            swWidths.(strains{strainCtr})(fileCtr,skelCtr) = swAreas.(strains{strainCtr})(fileCtr,skelCtr)...
                /swLengths.(strains{strainCtr})(fileCtr,skelCtr);
        end
        
        %% normalise area and perimeter from this movie with sw features from this movie; store value for threshold box plot later
        perimeterThres.(strains{strainCtr})(fileCtr) = numel(find(perimeter{fileCtr}/median(swPerimeters.(strains{strainCtr})(fileCtr,:))>perimeterCutOff))/numel(perimeter{fileCtr});
        areaThres.(strains{strainCtr})(fileCtr) = numel(find(area{fileCtr}/median(swAreas.(strains{strainCtr})(fileCtr,:))>areaCutOff))/numel(area{fileCtr});
    end
    
    %% pool data across movies
    swLength = median(swLengths.(strains{strainCtr})(:));
    swWidth = median(swWidths.(strains{strainCtr})(:));
    swArea = median(swAreas.(strains{strainCtr})(:));
    swPerimeter = median(swPerimeters.(strains{strainCtr})(:));
    perimeter = vertcat(perimeter{:});
    area = vertcat(area{:});
    
   %% use worm skeleton lengths to normalise blob features
    perimeter = perimeter/swPerimeter;
    area = area/swArea;

    %% plot figures
    set(0,'CurrentFigure',swLengthFig)
    histogram(swLengths.(strains{strainCtr}),'Normalization','pdf','DisplayStyle','stairs')
    set(0,'CurrentFigure',swWidthFig)
    histogram(swWidths.(strains{strainCtr}),'Normalization','pdf','DisplayStyle','stairs')
    set(0,'CurrentFigure',swAreaFig)
    histogram(swAreas.(strains{strainCtr}),'Normalization','pdf','DisplayStyle','stairs')
    set(0,'CurrentFigure',swPerimeterFig)
    histogram(swPerimeters.(strains{strainCtr}),'Normalization','pdf','DisplayStyle','stairs')
    
    set(0,'CurrentFigure',perimeterFig)
    histogram(perimeter,'Normalization','pdf','DisplayStyle','stairs')
    set(0,'CurrentFigure',areaFig)
    histogram(area,'Normalization','pdf','DisplayStyle','stairs')
end


%% format and save figures
if numel(strains) == 4
    boxplotGroup = [ones(size(perimeterThres.(strains{1}))) 2*ones(size(perimeterThres.(strains{2}))) 3*ones(size(perimeterThres.(strains{3}))) 4*ones(size(perimeterThres.(strains{4})))];
    %
    set(0,'CurrentFigure',perimeterThresholdFig)
    boxplot([perimeterThres.(strains{1}),perimeterThres.(strains{2}),perimeterThres.(strains{3}),perimeterThres.(strains{4})],boxplotGroup,'notch','on','labels',legendList);
    title('perimeter')
    set(perimeterThresholdFig,'PaperUnits','centimeters')
    figurename = ['figures/perimeterThreshold_' num2str(perimeterCutOff)];
    if saveResults
        exportfig(perimeterThresholdFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
    %
    set(0,'CurrentFigure',areaThresholdFig)
    boxplot([areaThres.(strains{1}),areaThres.(strains{2}),areaThres.(strains{3}),areaThres.(strains{4})],boxplotGroup,'notch','on','labels',legendList);
    title('area')
    set(areaThresholdFig,'PaperUnits','centimeters')
    figurename = ['figures/areaThreshold_' num2str(areaCutOff)];
    if saveResults
        exportfig(areaThresholdFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
else
    warning('fix strain list for threshold box plots')
end

%% change legend format for the double mutant strain

% if strcmp(legendList{3},'daf22_npr1')
%     legendList{3} = 'daf22\_npr1'; % add back slash so n doesn't become subscript
% else
%     warning('need to rename daf22_npr1 to avoid subscript appearance in legend')
% end
% 
% 
% set(0,'CurrentFigure',swLengthFig)
% legend(legendList)
% xlabel('single worm length')
% ylabel('probability')
% set(swLengthFig,'PaperUnits','centimeters')
% figurename = 'figures/swLength';
% if saveResults
%     exportfig(swLengthFig,[figurename '.eps'],exportOptions)
%     system(['epstopdf ' figurename '.eps']);
%     system(['rm ' figurename '.eps']);
% end
% 
% set(0,'CurrentFigure',swWidthFig)
% legend(legendList)
% xlabel('single worm width')
% ylabel('probability')
% set(swWidthFig,'PaperUnits','centimeters')
% figurename = 'figures/swWidth';
% if saveResults
%     exportfig(swWidthFig,[figurename '.eps'],exportOptions)
%     system(['epstopdf ' figurename '.eps']);
%     system(['rm ' figurename '.eps']);
% end
% 
% set(0,'CurrentFigure',swAreaFig)
% legend(legendList)
% xlabel('single worm area')
% ylabel('probability')
% set(swAreaFig,'PaperUnits','centimeters')
% figurename = 'figures/swArea';
% if saveResults
%     exportfig(swWidthFig,[figurename '.eps'],exportOptions)
%     system(['epstopdf ' figurename '.eps']);
%     system(['rm ' figurename '.eps']);
% end
% 
% set(0,'CurrentFigure',swPerimeterFig)
% legend(legendList)
% xlabel('single worm perimeter')
% ylabel('probability')
% set(swPerimeterFig,'PaperUnits','centimeters')
% figurename = 'figures/swPerimeter';
% if saveResults
%     exportfig(swPerimeterFig,[figurename '.eps'],exportOptions)
%     system(['epstopdf ' figurename '.eps']);
%     system(['rm ' figurename '.eps']);
% end
% 
% 
% set(0,'CurrentFigure',perimeterFig)
% legend(legendList)
% xlabel('perimeter')
% ylabel('probability')
% xlim([0 20])
% figurename = 'figures/perimeter';
% if saveResults
%     exportfig(perimeterFig,[figurename '.eps'],exportOptions)
%     system(['epstopdf ' figurename '.eps']);
%     system(['rm ' figurename '.eps']);
% end
% 
% 
% set(0,'CurrentFigure',areaFig)
% legend(legendList)
% xlabel('area')
% ylabel('probability')
% xlim([0 20])
% figurename = 'figures/area';
% if saveResults
%     exportfig(areaFig,[figurename '.eps'],exportOptions)
%     system(['epstopdf ' figurename '.eps']);
%     system(['rm ' figurename '.eps']);
% end