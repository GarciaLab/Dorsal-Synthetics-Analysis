function [FractionsPerEmbryo,TimeOnsPerEmbryo] =  getPerEmbryoDataForFits(includeString,excludeString,numBins)
% DataTypePattern shluld be a string corresponding to a tab in
% dataStatus.xls

% this function is based on plotgreendata_slim_SA. Check that code first of
% you want to add functionalities such as getting the per embryo values of
% other metrics or getting per nucleus metrics.

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

fiducialTime = 6; %mins

for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

% define rules for what data we'll consider in the analysis and plots
dataTypeIncludeRule = contains({combinedCompiledProjects_allEnhancers.dataSet},includeString);
dataTypeExcludeRule = ~contains(lower({combinedCompiledProjects_allEnhancers.dataSet}),lower(excludeString));
dataTypeRule = dataTypeIncludeRule & dataTypeExcludeRule;
enhancerStruct = combinedCompiledProjects_allEnhancers(dataTypeRule &...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));

% % just gonna grab data for one enhancer
% enhancerStruct = combinedCompiledProjects_allEnhancers(...
%     [combinedCompiledProjects_allEnhancers.cycle]==12 &...
%     ({combinedCompiledProjects_allEnhancers.dataSet} == DataTypePattern1 |...
%       {combinedCompiledProjects_allEnhancers.dataSet} == DataTypePattern2  )&...
%     ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));

% bin nuclei
binLimits = linspace(0,4000,numBins);
for n = 1:length(enhancerStruct)
    FluoFeature = enhancerStruct(n).dorsalFluoFeature;
    bin = find(FluoFeature-binLimits<0,1,'first')-1;
    enhancerStruct(n).dorsalFluoBin3 = bin;
end

% these bins are represented in the data
coveredBins = unique([enhancerStruct.dorsalFluoBin3]);
binValues = binLimits(coveredBins);

% define some filters
minEmbryosPerBin = 3;
minNucleiPerEmbryoPerBin = 1;
minOnset = 2; % (min) earliest possible spot detection time to be counted
maxOnset = 8; %(min) latest possible spot detection time to be counted

%store everything in these arrays
%mean_maxFluo_acrossNuclei_perBin = [];
%se_maxFluo_acrossNuclei_perBin = [];
%mean_accFluo_acrossNuclei_perBin = [];
%se_accFluo_acrossNuclei_perBin = [];
%mean_fraction_acrossNuclei_perBin = [];
%
% mean_maxFluo_acrossEmbryos_perBin = [];
% se_maxFluo_acrossEmbryos_perBin = [];
% mean_accFluo_acrossEmbryos_perBin = [];
% se_accFluo_acrossEmbryos_perBin = [];
% mean_fraction_acrossEmbryos_perBin = [];
% se_fraction_acrossEmbryos_perBin = [];

OnNucleiPerBin = [];
OffNucleiPerBin = [];

FractionsPerEmbryo = nan(50,length(coveredBins));
TimeOnsPerEmbryo = nan(50,length(coveredBins));
% MaxFluoPerEmbryo = nan(50,length(coveredBins));
% AccFluoPerEmbryo = nan(50,length(coveredBins));
% DurationPerEmbryo = nan(50,length(coveredBins));

% now make one struct per bin
for b = 1:length(coveredBins)
    binID = coveredBins(b);
    binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin3]== binID);
    activeNuc_Bin = length([binStruct.particleTimeOn]);
    inactiveNuc_Bin = length(binStruct) - activeNuc_Bin;
    % take the means and errors across nuclei first because it's easy      
%     mean_maxFluo_acrossNuclei_perBin(b) = nanmean([binStruct.particleFluo95]);
%     se_maxFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleFluo95])./sqrt(activeNuc_Bin);
%     mean_accFluo_acrossNuclei_perBin(b) =  nanmean([binStruct.particleAccumulatedFluo]);
%     se_accFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleAccumulatedFluo])./sqrt(activeNuc_Bin);
%    mean_fraction_acrossNuclei_perBin(b) = activeNuc_Bin/length(binStruct);
    %filter spurious time ons due to errors
%     particlesTimeOns = [binStruct.particleTimeOn];
%     particlesTimeOns = particlesTimeOns(particlesTimeOns>minOnset);
%     particlesTimeOns = particlesTimeOns(particlesTimeOns<maxOnset);
%     mean_timeOn_acrossNuclei_perBin(b) = nanmean(particlesTimeOns);
%     se_timeOn_acrossNuclei_perBin(b) = nanstd(particlesTimeOns)./sqrt(activeNuc_Bin);
%     mean_duration_acrossNuclei_perBin(b) = nanmean([binStruct.particleDuration]);
%     se_duration_acrossNuclei_perBin(b) = nanstd([binStruct.particleDuration])./sqrt(activeNuc_Bin);    
%     totalRNA_acrossNuclei_perBin(b) = nansum([binStruct.particleAccumulatedFluo]);
    
    % save this for bootstrapping later
%     OnNucleiPerBin(b) = activeNuc_Bin;
%     OffNucleiPerBin(b) = inactiveNuc_Bin;
    
    % now deal with the mean and errors across embryos
    [uniqueEmbryos, ~, J]=unique({binStruct.prefix});
    occurences = histc(J, 1:numel(uniqueEmbryos));
    numEmbryos = length(occurences);
    qualityEmbryos = occurences >= minNucleiPerEmbryoPerBin;
%     maxFluo_perEmbryo = [];
%     accFluo_perEmbryo = [];
    fraction_perEmbryo = []; 
    nuclei_perEmbryo = [];
    timeOn_perEmbryo = [];
%     duration_perEmbryo = [];
%     totalRNA_perEmbryo = [];
%     
    % and now make one struct per embryo per bin
    if numEmbryos >= minEmbryosPerBin     
        for e = 1:numEmbryos
            embryoPrefix = uniqueEmbryos{e};
%            if qualityEmbryos(e) %if this prefix has enough nuclei
            embryoStruct = binStruct(strcmpi({binStruct.prefix},embryoPrefix));
            nuclei_perEmbryo(e) = length(embryoStruct);
%             maxFluo_perEmbryo(e) = nanmean([embryoStruct.particleFluo95]);
%             accFluo_perEmbryo(e) =  nanmean([embryoStruct.particleAccumulatedFluo]);
            fraction_perEmbryo(e) = length([embryoStruct.particleTimeOn])/length(embryoStruct);
            % filter the onset times to 2<true<8 (min since anaphase)
            perEmbryoTimeOns = [embryoStruct.particleTimeOn];
            perEmbryoTimeOns = perEmbryoTimeOns(perEmbryoTimeOns>minOnset);
            perEmbryoTimeOns = perEmbryoTimeOns(perEmbryoTimeOns<maxOnset);
            timeOn_perEmbryo(e) = nanmean(perEmbryoTimeOns);
%             duration_perEmbryo(e) = nanmean([embryoStruct.particleDuration]);
%             totalRNA_perEmbryo(e) =  nansum([embryoStruct.particleAccumulatedFluo]).*fraction_perEmbryo(e);
%            end
        end
    end
    
    %add things taken across embryos to the arrays where we store stuff
%     mean_maxFluo_acrossEmbryos_perBin(b) = nanmean(maxFluo_perEmbryo);
%     se_maxFluo_acrossEmbryos_perBin(b) = nanstd(maxFluo_perEmbryo)./sqrt(numEmbryos);
%     mean_accFluo_acrossEmbryos_perBin(b) = nanmean(accFluo_perEmbryo);
%     se_accFluo_acrossEmbryos_perBin(b) = nanstd(accFluo_perEmbryo)./sqrt(numEmbryos);
    mean_fraction_acrossEmbryos_perBin(b) = nanmean(fraction_perEmbryo);
    se_fraction_acrossEmbryos_perBin(b) = std(fraction_perEmbryo)./sqrt(numEmbryos);
    mean_timeOn_acrossEmbryos_perBin(b) = nanmean(timeOn_perEmbryo);
    se_timeOn_acrossEmbryos_perBin(b) = nanstd(timeOn_perEmbryo)./sqrt(numEmbryos);
%     mean_duration_acrossEmbryos_perBin(b) = nanmean(duration_perEmbryo);
%     se_duration_acrossEmbryos_perBin(b) = nanstd(duration_perEmbryo)./sqrt(numEmbryos);
%     mean_totalRNA_acrossEmbryos_perBin(b) = nanmean(totalRNA_perEmbryo);
%     se__totalRNA_acrossEmbryos_perBin(b) = nanstd(totalRNA_perEmbryo)./sqrt(numEmbryos);
    
    %save the single embryo datum for later use
    FractionsPerEmbryo(1:length(fraction_perEmbryo),b) = fraction_perEmbryo;
    TimeOnsPerEmbryo(1:length(timeOn_perEmbryo),b) = timeOn_perEmbryo;
%     MaxFluoPerEmbryo(1:length(maxFluo_perEmbryo),b) = maxFluo_perEmbryo;
%     AccFluoPerEmbryo(1:length(accFluo_perEmbryo),b) = accFluo_perEmbryo;
%     DurationPerEmbryo(1:length(fraction_perEmbryo),b) = duration_perEmbryo;
    %clear  nuclei_per_Embryo maxFluo_perEmbryo accFluo_perEmbryo fraction_perEmbryo timeOn_perEmbryo
     
end
