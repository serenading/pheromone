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
strains = {'daf22_npr1','daf22','npr1','N2'}; % {'daf22_npr1','daf22','npr1','N2'}
compactnessFig = figure; hold on
perimeterFig = figure; hold on
areaFig = figure; hold on
quirkinessFig = figure; hold on
compactnessTest = cell(length(strains),1);
perimeterTest = cell(length(strains),1);
areaTest = cell(length(strains),1);
quirkinessTest = cell(length(strains),1);
compactnessCat = cell(length(strains),1);
perimeterCat = cell(length(strains),1);
areaCat = cell(length(strains),1);
quirkinessCat = cell(length(strains),1);

%% go through strains, densities, movies
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    filenames = importdata(['datalist/' strain '_list.txt']);
    numFiles = length(filenames);
    compactness = cell(numFiles,1);
    perimeter = cell(numFiles,1);
    area = cell(numFiles,1);
    quirkiness = cell(numFiles,1);
    %% go through individual movies
    for fileCtr = 1:numFiles
        %% load data
        filename = filenames{fileCtr}
        %trajData = h5read(filename,'/trajectories_data');
        blobFeats = h5read(filename,'/blob_features');
        %skelData = h5read(filename,'/skeleton');
        frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
        compactness{fileCtr} = blobFeats.compactness;
        perimeter{fileCtr} = blobFeats.perimeter;
        area{fileCtr} = blobFeats.area;
        quirkiness{fileCtr} = blobFeats.quirkiness;
    end
    compactness = vertcat(compactness{:});
    perimeter = vertcat(perimeter{:});
    area = vertcat(area{:});
    quirkiness = vertcat(quirkiness{:});
    %% save values for Kruskal-Wallis nonparametric test
    compactnessTest{strainCtr} = compactness;
    perimeterTest{strainCtr} = perimeter;
    areaTest{strainCtr} = area;
    quirkinessTest{strainCtr} = quirkiness;
    %% plot figures
    set(0,'CurrentFigure',compactnessFig)
    histogram(compactness,'Normalization','pdf','DisplayStyle','stairs')
    set(0,'CurrentFigure',perimeterFig)
    histogram(perimeter,'Normalization','pdf','DisplayStyle','stairs')
    set(0,'CurrentFigure',areaFig)
    histogram(area,'Normalization','pdf','DisplayStyle','stairs')
    set(0,'CurrentFigure',quirkinessFig)
    histogram(quirkiness,'Normalization','pdf','DisplayStyle','stairs')
end

%% Kruskal-Wallis nonparametric test
compactnessTestPool = vertcat(compactnessTest{:});
perimeterTestPool = vertcat(perimeterTest{:});
areaTestPool = vertcat(areaTest{:});
quirkinessTestPool = vertcat(quirkinessTest{:});
compactnessTestCat = NaN(length(compactnessTestPool),1);
perimeterTestCat = NaN(length(perimeterTestPool),1);
areaTestCat = NaN(length(areaTestPool),1);
quirkinessTestCat = NaN(length(quirkinessTestPool),1);
positionCtr = 1;
% generate category labels the same length as the saved values
for strainCtr = 1:length(strains)
    strain = strains{strainCtr};
    compactnessTestCat(positionCtr:(positionCtr+length(compactnessTest{strainCtr})-1)) = strainCtr;
    perimeterTestCat(positionCtr:(positionCtr+length(compactnessTest{strainCtr})-1)) = strainCtr;
    areaTestCat(positionCtr:(positionCtr+length(compactnessTest{strainCtr})-1)) = strainCtr;
    quirkinessTestCat(positionCtr:(positionCtr+length(compactnessTest{strainCtr})-1)) = strainCtr;
    positionCtr = positionCtr + length(compactnessTest{strainCtr});
end
% the actual test
compactnessp = kruskalwallis(compactnessTestPool,compactnessTestCat)
perimeterp = kruskalwallis(perimeterTestPool,perimeterTestCat)
areap = kruskalwallis(areaTestPool,areaTestCat)
quirkinessp = kruskalwallis(quirkinessTestPool,quirkinessTestCat)

%% format and save figures
set(0,'CurrentFigure',compactnessFig)
legend('daf22\_npr1','daf22','npr1','N2')
xlabel('compactness')
ylabel('probability')
xlim([0 1])
set(compactnessFig,'PaperUnits','centimeters')
figurename = 'figures/compactness';
%savefig(compactnessFig,[figurename '.fig'])
%exportfig(compactnessFig,[figurename '.eps'],exportOptions)
%
set(0,'CurrentFigure',perimeterFig)
legend('daf22\_npr1','daf22','npr1','N2')
xlabel('perimeter')
ylabel('probability')
xlim([0 1500])
figurename = 'figures/perimeter';
%savefig(perimeterFig,[figurename '.fig'])
%exportfig(perimeterFig,[figurename '.eps'],exportOptions)
%
set(0,'CurrentFigure',areaFig)
legend('daf22\_npr1','daf22','npr1','N2')
xlabel('area')
ylabel('probability')
xlim([0 2000])
figurename = 'figures/area';
%savefig(areaFig,[figurename '.fig'])
%exportfig(areaFig,[figurename '.eps'],exportOptions)
%
set(0,'CurrentFigure',quirkinessFig)
legend('daf22\_npr1','daf22','npr1','N2')
xlabel('quirkiness')
ylabel('probability')
xlim([0 1])
figurename = 'figures/quirkiness';
%savefig(quirkinessFig,[figurename '.fig'])
%exportfig(quirkinessFig,[figurename '.eps'],exportOptions)