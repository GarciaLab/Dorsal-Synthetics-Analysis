function [x, y] = plotGreenBoxWithData

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

numBins = 20;

fiducialTime = 6;

for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

enhancerStruct = combinedCompiledProjects_allEnhancers(...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));


% bin nuclei
%this function calculates the Dorsal fluorescence at some arbitrary time
%in nc12 and adds it to the struct in a 'DorsalFluoArbitraryTime' field
enhancerStruct = DorsalFluoArbitraryTime(enhancerStruct,fiducialTime);
nucleiFluorescence = [enhancerStruct.DorsalFluoArbitraryTime];
nucleiFluorescence = [enhancerStruct.dorsalFluoFeature];

binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(enhancerStruct)
    enhancerStruct(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end

coveredBins = unique([enhancerStruct.dorsalFluoBin2]);
binValues = binValues(coveredBins);
%store everything in these arrays
mean_fraction_acrossNuclei_perBin = [];
mean_fraction_acrossEmbryos_perBin = [];
se_fraction_acrossEmbryos_perBin = [];

OnNucleiPerBin = [];
OffNucleiPerBin = [];

dataSets = unique({enhancerStruct.dataSet});
prefixes = unique({enhancerStruct.prefix});
a = [];
c = [];
for p = 1:length(prefixes)
    
    mean_fraction_acrossNuclei_perBin = [];
    mean_timeOn_acrossNuclei_perBin = [];
    se_timeOn_acrossNuclei_perBin = [];
    
    for b = 1:length(coveredBins)
        
        binID = coveredBins(b);
        binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin2]== binID & strcmpi({enhancerStruct.prefix}, prefixes{p}) );
        
        activeNuc_Bin = length([binStruct.particleTimeOn]);
        inactiveNuc_Bin = length(binStruct) - activeNuc_Bin;
        % take the means and errors across nuclei first because it's easy
        mean_fraction_acrossNuclei_perBin(b) = activeNuc_Bin/length(binStruct);
        mean_timeOn_acrossNuclei_perBin(b) = nanmean([binStruct.particleTimeOn]);
        se_timeOn_acrossNuclei_perBin(b) = nanstd([binStruct.particleTimeOn])./sqrt(activeNuc_Bin);
    end
    a = [a, mean_fraction_acrossNuclei_perBin];
    c = [c, mean_timeOn_acrossNuclei_perBin];
end

x = c;
y = a;


%clean data
y(x<2) = [];
x(x < 2) = [];
y(x>8) = [];
x(x>8) = [];
y(isnan(x)) = [];
x(isnan(x)) = [];
x(isnan(y)) = [];
y(isnan(y)) = [];

P = [x;y]';

nBins = [15, 5];

displayFigures = false;

if displayFigures
    h = binscatter(x, y, nBins)
    h.ShowEmptyBins = 'on';
    xlabel('mean turn on')
    ylabel('fraction active')
    ax = gca;
    hold on
    [k,	~] = convhull([x;y]');
    % plot(x(k),y(k))
    hold on
    % fill(x(k),y(k), 'g', 'FaceAlpha', .3, 'LineWidth', 3, 'EdgeColor', 'g')
    % hull = alphaShape(x_hull, y_hull,Inf,'HoleThreshold',1E30 );
    % colormap brewermap(20,'Blues')
    % g = histogram2(x, y,nBins, 'DisplayStyle','tile','ShowEmptyBins','on')
    
    size(h.Values)
    h.Values
    h.XBinEdges
    h.YBinEdges
    
    % [row,col,v] = find(h.Values > 10)
    
end

%%


