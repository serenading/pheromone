% A. Script visualises sweeping using pixel data from selected frames, using up to 5 hours of long sweeping brightfield recordings (Leah's 1 patch dataset)
% step 1: read image frame
% step 2: generate binary image based on intensity threshold to pick out worm pixels from background
% step 3: dilate image to "connect" loose pharynxes and apply area thresholding to binary image to pick up clusters
% step 4: draw clusters over time (optional: plot centroid)
% step 5: plot food contour on top of the image (optional)
% B,C. Script calculates cluster centroid speed and plots them as timeseries and box plot.
% D. Script plots MSD of cluster centroids
% E. Script generates a down-sampled avi video, stringing together segments over multiple 1-hour recordings


close all
clear

dataset = 'long'; % 'short' or 'long'
strains = {'daf22_npr1','daf22','npr1','N2'};
sampleEveryNSec = 120;  % in seconds
dafblobAreaThreshold = 2000;
nondafBlobAreaThreshold = 3000; % single worm area ~ 500
plotVisualisation = true;
saveCentroidValues = false;
plotCentroidSpeeds = false;
plotMeanSquaredDisplacement = false;
makeDownSampledVideo = false;

if strcmp(dataset,'long')
    maxSeg = 15; % maximum number of 1-hour recordings
elseif strcmp(dataset,'short')
    maxSeg = 2; % maximum number of 1-hour recordings
end

if plotVisualisation
    plotCentroid = false;
    plotFoodContour = true;
end

if plotCentroidSpeeds || plotMeanSquaredDisplacement
    smoothWindow = 20; % number of sampled frames for smoothing i.e. smoothWindow = 20x sampleEveryNSec = 30 means smoothing over 600 seconds
    plotFileList = [1,2,3,4,6]; % the replicates to plot for overall cluster speed - ignore rep 5 with two converging clusters
end

if plotMeanSquaredDisplacement
    initialTimeStepToUse = 60; % use file matching this timeStep to start with
end

pixelsize = 10; % microns per pixel. can read by pixelsize = double(h5readatt(skelFilename,'/trajectories_data','microns_per_pixel'));
frameRate = 25; % frames per second. can read by frameRate = double(h5readatt(skelFilename,'/plate_worms','expected_fps'));

exportOptions = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

exportOptions2 = struct('Format','EPS2',...
    'Color','rgb',...
    'Width',25,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1); % for individual rep timeseries only

addpath('auxiliary/')

for strainCtr = 1:length(strains)
    if strcmp(dataset,'long')
        [annotationNum,annotationFilenames,~] = xlsread('datalist/pheromoneLong.xlsx',strainCtr,'A1:E200','basic');
        % xy coordinates and radius of food contour obtained by hand annotation using VGG
        foodCtnCoords_xyr.daf22_npr1 = [1009,935,406;1023,993,429;975,1083,427;998,981,392;...
            1080,1068,409;1175,923,422;1077,993,419;1132,980,395];
        foodCtnCoords_xyr.daf22 = [1118,1060,412;1134,924,397;1059,1017,413;1149,912,409;...
            1019,1057,387;1058,1039,397;1090,904,396;1092,1042,399];
        foodCtnCoords_xyr.npr1 = [1166,926,421;993,876,422;951,783,403;900,875,383;...
            963,915,428;980,1115,399;1085,1112,411;1204,858,396];
        foodCtnCoords_xyr.N2 = [1097,1008,406;996,1077,425;1053,948,402;1212,936,398;...
            1032,1043,427;1182,1038,396;940,1027,426;1247,864,401];
    elseif strcmp(dataset,'short')
        [annotationNum,annotationFilenames,~] = xlsread('datalist/pheromoneTwoHours.xlsx',strainCtr,'A1:E12','basic');
        % xy coordinates and radius of food contour obtained by hand annotation using VGG
        if strcmp(strains{strainCtr},'daf22_npr1')
            foodCtnCoords_xyr = [1015,986,362;1006,963,367;1272,1013,371;942,988,363;1062,961,359];
        end
    end

    % go through each recording replicate
    for fileCtr = 1:max(annotationNum(:,1))
        if plotVisualisation
            clusterVisFig = figure; hold on
        end
        if makeDownSampledVideo
            video = VideoWriter(['pheromone_' dataset '_' strains{strainCtr} '_rep' num2str(fileCtr) '_timeStep' num2str(sampleEveryNSec) '.avi']); % create the video object
            video.FrameRate = 15;% set the frame rate
            open(video); % open the file for writing
        end
        totalFrames = 0;
        totalSegs = 0;
        for segCtr = 1:max(annotationNum(:,2)) % go through each hour of the recording replicate
            fileIdx = find(annotationNum(:,1) == fileCtr & annotationNum(:,2) == segCtr);
            firstFrame{segCtr} = annotationNum(fileIdx,4)+1; % +1 to adjust for python 0 indexing
            lastFrame{segCtr} = annotationNum(fileIdx,5)+1;
            if lastFrame{segCtr} - firstFrame{segCtr} > 0 % if this recording has any valid frames
                totalFrames = totalFrames+lastFrame{segCtr}-firstFrame{segCtr}+1;
                totalSegs = totalSegs+1;
                filename{segCtr} = annotationFilenames{fileIdx};
            end
        end
        totalSampleFrames = ceil(totalFrames/sampleEveryNSec/25);
        % initialise
        plotColors = parula(totalSampleFrames);
        clusterCentroidCoords{fileCtr} = NaN(2,totalSampleFrames);
        clusterSolidity{fileCtr} = NaN(1,totalSampleFrames);
        clusterArea{fileCtr} = NaN(1,totalSampleFrames);
        cumFrame = 0; % keep track of cumulative frames across replicate segments
        leftoverFrames = 0; % keep track of leftover frames at the end of one segment that combines with the start of the next segment
        
        for segCtr = 1:totalSegs
            % load data
            fileInfo = h5info(filename{segCtr});
            dims = fileInfo.Datasets(2).Dataspace.Size;
            if leftoverFrames>0
                assert(firstFrame{segCtr} ==1); % if there are leftover frames from the previous segment, then this segment must start from the very first frame
                firstFrame{segCtr} = sampleEveryNSec*frameRate-leftoverFrames+firstFrame{segCtr};
            end
            movieFrames = firstFrame{segCtr}:sampleEveryNSec*frameRate:...
                floor((lastFrame{segCtr}-firstFrame{segCtr}+1)/sampleEveryNSec/frameRate)*sampleEveryNSec*frameRate+firstFrame{segCtr};
            for frameCtr = 1:numel(movieFrames)
                imageFrame = h5read(filename{segCtr},'/mask',[1,1,movieFrames(frameCtr)],[dims(1),dims(2),1]);
                if makeDownSampledVideo
                    timeStamp = round((cumFrame+frameCtr-1)*sampleEveryNSec/60); % timestamp in minutes
                    imageFrame = AddTextToImage(imageFrame,['t = ' num2str(timeStamp) ' min'],[100 100],[255,255,255],'Arial',50);
                    writeVideo(video,imageFrame); %write the image to file
                end
                % generate binary segmentation based on black/white contrast
                binaryImage = imageFrame>0 & imageFrame<70;
                binaryImage = imfill(binaryImage, 'holes');
                % filter by blob size
                blobMeasurements = regionprops(binaryImage, 'Area','Centroid','Solidity');
                blobCentroidsCoords = reshape([blobMeasurements.Centroid],[2, numel([blobMeasurements.Centroid])/2]);
                if strainCtr >2
                    blobLogInd = [blobMeasurements.Area] > nondafBlobAreaThreshold; % apply blob area threshold values
                else
                    blobLogInd = [blobMeasurements.Area] > dafblobAreaThreshold;
                end
                blobLogInd = blobLogInd & blobCentroidsCoords(2,:) > 250; % get rid of the annoying box at the edge
                % restrict to blobs near the food patch centre (within 500 pixels or 5 mm)
                if strcmp(dataset,'short') & strcmp(strains{strainCtr},'daf22_npr1')
                   blobLogInd = blobLogInd & blobCentroidsCoords(1,:)<foodCtnCoords_xyr(fileCtr,1)+500 & blobCentroidsCoords(1,:)>foodCtnCoords_xyr(fileCtr,1)-500;
                   blobLogInd = blobLogInd & blobCentroidsCoords(2,:)<foodCtnCoords_xyr(fileCtr,2)+500 & blobCentroidsCoords(2,:)>foodCtnCoords_xyr(fileCtr,2)-500;
                elseif strcmp(dataset,'long')
                    blobLogInd = blobLogInd & blobCentroidsCoords(1,:)<foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,1)+500 & blobCentroidsCoords(1,:)>foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,1)-500;
                    blobLogInd = blobLogInd & blobCentroidsCoords(2,:)<foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,2)+500 & blobCentroidsCoords(2,:)>foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,2)-500;
                end
                blobBoundaries = bwboundaries(binaryImage,8,'noholes');
                if plotVisualisation
                    % plot individual blob boundaries that meet area threshold requirements
                    set(0,'CurrentFigure',clusterVisFig)
                    for blobCtr = 1:numel(blobLogInd)
                        if blobLogInd(blobCtr)
                            fill(blobBoundaries{blobCtr}(:,1)*pixelsize/1000,blobBoundaries{blobCtr}(:,2)*pixelsize/1000,plotColors(cumFrame+frameCtr,:),'edgecolor','none')
                            alpha 0.5
                        end
                    end
                end
                % get centroids and solidity
                if nnz(blobLogInd)>0
                    [~,maxAreaIdx] = max([blobMeasurements.Area].*blobLogInd);
                    clusterCentroidCoords{fileCtr}(1,cumFrame+frameCtr) = blobCentroidsCoords(1,maxAreaIdx);
                    clusterCentroidCoords{fileCtr}(2,cumFrame+frameCtr) = blobCentroidsCoords(2,maxAreaIdx);
                    area = [blobMeasurements.Area];
                    solidity = [blobMeasurements.Solidity];
                    clusterSolidity{fileCtr}(cumFrame+frameCtr) = solidity(maxAreaIdx);
                    clusterArea{fileCtr}(cumFrame+frameCtr) = area(maxAreaIdx);
                    if plotVisualisation & plotCentroid
                        set(0,'CurrentFigure',clusterVisFig)
                        plot(clusterCentroidCoords{fileCtr}(2,cumFrame+frameCtr)*pixelsize/1000,clusterCentroidCoords{fileCtr}(1,cumFrame+frameCtr)*pixelsize/1000,'k--x')
                    end
                end
            end
            cumFrame = cumFrame+numel(movieFrames);
            leftoverFrames = lastFrame{segCtr} - max(movieFrames);
        end
        
        if plotVisualisation
            set(0,'CurrentFigure',clusterVisFig)
            axis equal
            colorbar
            caxis([0 ceil(totalFrames/25/60)])
            cb = colorbar; cb.Label.String = 'minutes';
            if strcmp(dataset,'short') & strcmp(strains{strainCtr},'daf22_npr1')
                xmax = round(foodCtnCoords_xyr(fileCtr,2)*pixelsize/1000+5);
                xmin = round(foodCtnCoords_xyr(fileCtr,2)*pixelsize/1000-5);
                ymax = round(foodCtnCoords_xyr(fileCtr,1)*pixelsize/1000+5);
                ymin = round(foodCtnCoords_xyr(fileCtr,1)*pixelsize/1000-5);
            elseif strcmp(dataset,'long')
                xmax = round(foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,2)*pixelsize/1000+5);
                xmin = round(foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,2)*pixelsize/1000-5);
                ymax = round(foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,1)*pixelsize/1000+5);
                ymin = round(foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,1)*pixelsize/1000-5);
            end
            xlim([xmin xmax])
            ylim([ymin ymax])
            xticks(xmin:2:xmax)
            yticks(ymin:2:ymax)
            xlabel('x (mm)')
            ylabel('y (mm)')
            if plotFoodContour
                if strcmp(dataset,'short') &strcmp(strains{strainCtr},'daf22_npr1')
                    viscircles([foodCtnCoords_xyr(fileCtr,2),foodCtnCoords_xyr(fileCtr,1)]*pixelsize/1000,foodCtnCoords_xyr(fileCtr,3)*pixelsize/1000,'Color','k','LineStyle','--','LineWidth',1);
                elseif strcmp(dataset,'long')
                    viscircles([foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,2),foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,1)]*pixelsize/1000,foodCtnCoords_xyr.(strains{strainCtr})(fileCtr,3)*pixelsize/1000,'Color','k','LineStyle','--','LineWidth',1);
                end
            end
            % export figure
            if plotCentroid
                if plotFoodContour
                    if strainCtr <=2
                        figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroidFood_blobArea' num2str(dafblobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFPhe_' dataset];
                    else
                        figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroidFood_blobArea' num2str(nondafBlobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFPhe_' dataset];
                    end
                else
                    if strainCtr <=2
                        figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroid_blobArea' num2str(dafblobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFPhe_' dataset];
                    else
                        figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelCentroid_blobArea' num2str(nondafBlobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFPhe_' dataset];
                    end
                end
            elseif plotFoodContour
                if strainCtr<=2
                    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelFood_blobArea' num2str(dafblobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFPhe_' dataset];
                else
                    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixelFood_blobArea' num2str(nondafBlobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFPhe_' dataset];
                end
            else
                if strainCtr<=2
                    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixel_blobArea' num2str(dafblobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFPhe_' dataset];
                else
                    figurename = ['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_rep' num2str(fileCtr) '_blobsOverTimePixel_blobArea' num2str(nondafBlobAreaThreshold) '_timeStep' num2str(sampleEveryNSec) '_dataBFPhe_' dataset];
                end
            end
            exportfig(clusterVisFig,[figurename '.eps'],exportOptions)
        end
        
        % calculate centroid speed (in microns per minute)
        clusterCentroidSpeed{fileCtr} = NaN(1,totalSampleFrames);
        for frameCtr = 1:totalSampleFrames-10
            if ~isnan(clusterCentroidCoords{fileCtr}(1,frameCtr))
                for stepCtr = 1:10 % in case the next sample frame has no cluster, go up to 10 time steps away
                    if ~isnan(clusterCentroidCoords{fileCtr}(1,frameCtr+stepCtr))
                        break
                    end
                end
                clusterCentroidSpeed{fileCtr}(frameCtr) = sqrt((clusterCentroidCoords{fileCtr}(1,frameCtr+stepCtr)-clusterCentroidCoords{fileCtr}(1,frameCtr))^2 +...
                    (clusterCentroidCoords{fileCtr}(2,frameCtr+stepCtr)-clusterCentroidCoords{fileCtr}(2,frameCtr))^2)...
                    /stepCtr*pixelsize*60/sampleEveryNSec;
                if clusterCentroidSpeed{fileCtr}(frameCtr)>500
                    clusterCentroidSpeed{fileCtr}(frameCtr) = NaN;
                end
            end
        end
        if makeDownSampledVideo
            close(video); %close the file
        end
    end
    
    % save cluster centroid coordinates and speeds
    
    if saveCentroidValues
        save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidSpeed_dataBFPhe_' dataset '_' strains{strainCtr} '_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterCentroidSpeed');
        save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidCoords_dataBFPhe_' dataset '_' strains{strainCtr} '_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterCentroidCoords');
        save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidSolidity_dataBFPhe_' dataset '_' strains{strainCtr} '_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterSolidity');
        save(['figures/sweeping/' strains{strainCtr} '_clusterCentroidArea_dataBFPhe_' dataset '_' strains{strainCtr} '_timeStep' num2str(sampleEveryNSec) '.mat'],'clusterArea');
    end
end

%% plot median speeds for different replicates (npr-1 or daf-22;npr1 only)
if plotCentroidSpeeds
    load(['figures/sweeping/npr1_clusterCentroidSpeed_dataBFPhe_' dataset '_' strains{strainCtr} '_timeStep' num2str(sampleEveryNSec) '.mat'])
    smoothSpeeds = NaN(numel(plotFileList),maxSeg*3600/sampleEveryNSec); % initialise
    recordingColors = distinguishable_colors(numel(plotFileList));
    legends = cell(1,numel(plotFileList));
    
    clusterCentroidSpeedFig = figure; hold on % show each replicate individually
    poolRepFig = figure; % shaded error bars showing average across each specified replicate
    smoothedBoxPlotFig = figure; % show each replicate as a box plot
    unsmoothedBlotPlotFig = figure; % show each replicate as a box plot
    speedVSolidityFig = figure; % plots speed against solidity
    
    for fileCtr = 1:length(plotFileList)
        fileIdx = plotFileList(fileCtr);
        totalFrames = 0;
        totalSegs = 0;
        for segCtr = 1:5 % go through each hour of the recording replicate (5 hours maximum)
            recIdx = find(annotationNum(:,1) == fileCtr & annotationNum(:,2) == segCtr);
            firstFrame{segCtr} = annotationNum(recIdx,4)+1; % +1 to adjust for python 0 indexing
            lastFrame{segCtr} = annotationNum(recIdx,5)+1;
            filename{segCtr} = annotationFilenames{recIdx};
            if lastFrame{segCtr} - firstFrame{segCtr} > 0 % if this recording has any valid frames
                totalFrames = totalFrames+lastFrame{segCtr}-firstFrame{segCtr}+1;
                totalSegs = totalSegs+1;
            end
        end
        recordingsPlotX = 1:sampleEveryNSec/60:ceil(totalFrames/25/60);
        set(0,'CurrentFigure',clusterCentroidSpeedFig)
        if numel(recordingsPlotX) == numel(clusterCentroidSpeed{fileIdx})
            plot(recordingsPlotX,smoothdata(clusterCentroidSpeed{fileIdx},'movmedian',smoothWindow,'omitnan'),'Color',recordingColors(fileCtr,:))
        elseif numel(recordingsPlotX) < numel(clusterCentroidSpeed{fileIdx})
            plot(recordingsPlotX,smoothdata(clusterCentroidSpeed{fileIdx}(1:numel(recordingsPlotX)),'movmedian',smoothWindow,'omitnan'),'Color',recordingColors(fileCtr,:))
        elseif numel(recordingsPlotX) > numel(clusterCentroidSpeed{fileIdx})
            plot(recordingsPlotX(1:numel(clusterCentroidSpeed{fileIdx})),smoothdata(clusterCentroidSpeed{fileIdx},'movmedian',smoothWindow,'omitnan'),'Color',recordingColors(fileCtr,:))
        end
        legends{fileCtr} = num2str(fileIdx);
        %figure; plot(smoothdata(clusterCentroidSpeed{fileIdx},'movmedian',smoothWindow,'omitnan'),clusterSolidity{fileIdx},'.')
    end
    
    set(0,'CurrentFigure',clusterCentroidSpeedFig)
    xlabel('Time (min)'), ylabel('Cluster Speed (microns/min)')
    legend(legends,'Location','eastoutside')
    xlim([0 60*maxSeg])
    set(gca,'Xtick',[0:30:60*maxSeg])
    figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBFPhe_' dataset '_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_individualReps'];
    exportfig(clusterCentroidSpeedFig,[figurename '.eps'],exportOptions2)
    
    set(0,'CurrentFigure',poolRepFig)
    recordingsPlotX = 1:sampleEveryNSec/60:60*maxSeg;
    shadedErrorBar(recordingsPlotX,nanmedian(smoothSpeeds,1),nanstd(smoothSpeeds),'k');
    xlabel('Time (min)'), ylabel('Cluster Speed (microns/min)')
    xlim([0 250])
    set(gca,'Xtick',[0:30:60*maxSeg])
    figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBFPhe_' dataset '_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_pooledReps'];
    exportfig(poolRepFig,[figurename '.eps'],exportOptions)
    
    set(0,'CurrentFigure',smoothedBoxPlotFig)
    boxGroups = [];
    for fileCtr = 1:length(plotFileList)
        fileIdx = plotFileList(fileCtr);
        boxGroups = [boxGroups fileIdx*ones(1,3600*maxSeg/sampleEveryNSec)];
    end
    boxplot(smoothSpeeds(:),boxGroups(:))
    set(0,'CurrentFigure',smoothedBoxPlotFig)
    xlabel('Replicate'), ylabel('Cluster Speed (microns/min)')
    set(gca,'XTickLabel',legends)
    figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBFPhe_' dataset '_timeStep' num2str(sampleEveryNSec) '_smoothWindow' num2str(smoothWindow) '_smoothedBoxPlot'];
    exportfig(smoothedBoxPlotFig,[figurename '.eps'],exportOptions)
    
    set(0,'CurrentFigure',unsmoothedBlotPlotFig)
    unsmoothedSpeeds = [];
    unsmoothedGroups = [];
    for fileCtr = 1:length(plotFileList)
        fileIdx = plotFileList(fileCtr);
        unsmoothedSpeeds = [unsmoothedSpeeds clusterCentroidSpeed{fileIdx}];
        unsmoothedGroups = [unsmoothedGroups fileIdx*ones(1,numel(clusterCentroidSpeed{fileIdx}))];
    end
    boxplot(unsmoothedSpeeds(:),unsmoothedGroups(:))
    xlabel('Replicate'), ylabel('Cluster Speed (microns/min)')
    figurename = ['figures/sweeping/npr1_clusterCentroidSpeed_dataBFPhe_' dataset '_timeStep' num2str(sampleEveryNSec) '_unsmoothedBlotPlot'];
    exportfig(unsmoothedBlotPlotFig,[figurename '.eps'],exportOptions)
    
end

if plotMeanSquaredDisplacement
    % read xy coordinates from specified initial timeStep file
    load(['/Users/sding/Documents/pheromone/figures/sweeping/npr1_clusterCentroidCoords_dataBFPhe_' dataset '_timeStep' num2str(initialTimeStepToUse) '.mat'])
    % initialise
    msdFig = figure; hold on
    legends = cell(1,numel(plotFileList));
    
    for fileCtr = 1:length(plotFileList)
        fileIdx = plotFileList(fileCtr);
        % smooth xy coordinates
        smoothClusterCentroidCoords = smoothdata(clusterCentroidCoords{fileIdx},1,'movmedian',smoothWindow,'omitnan');
        % interpolate over NaN xy coordinates
        xcoords = naninterp(smoothClusterCentroidCoords(1,:))';
        ycoords = naninterp(smoothClusterCentroidCoords(2,:))';
        data=[xcoords ycoords];
        % initialise
        nData = size(data,1); %number of data points
        numberOfDeltaT = floor(nData/4); %for MSD, dt should be up to 1/4 of number of data points (Saxton, M. J. (1997). Single-particle tracking: The distribution of diffusion coeffcients. Biophys. J. 72, 1744?1753.);
        msd = zeros(numberOfDeltaT,3); %We'll store [mean, std, n]
        % calculate msd for all deltaT's
        for dt = 1:numberOfDeltaT
            deltaCoords = data(1+dt:end,1:2) - data(1:end-dt,1:2);
            squaredDisplacement = sum(deltaCoords.^2,2); % dx^2+dy^2
            msd(dt,1) = mean(squaredDisplacement)*pixelsize/1e6; % average in cm
            msd(dt,2) = std(squaredDisplacement)*pixelsize/1e6; % std in cm
            msd(dt,3) = length(squaredDisplacement); % n
        end
        % plot MSD
        plot((1:size(msd,1))*initialTimeStepToUse/60,msd(:,1))
        % generate legend
        legends{fileCtr} = num2str(fileIdx);
    end
    % format MSD plot
    xlabel('tau (min)'); ylabel('MSD (cm)')
    legend(legends,'Location','eastoutside')
    figurename = ['figures/sweeping/npr1_clusterMSD_dataBFPhe_' dataset '_initialTimeStep' num2str(initialTimeStepToUse)];
    exportfig(msdFig,[figurename '.eps'],exportOptions)
end