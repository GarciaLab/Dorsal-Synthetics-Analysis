function [embryoRNA embryoRNAError] = averagesTake2(DataType,numBins,metric,fiducialTime,errorgroup,Color,ax)
% metric can be 'maxfluo', 'accumulatedfluo' or 'fraction'
% errorgroup is over what the error is taken, 'embryos' or 'nuclei'. For
% fraction the error over nuclei is bootstraped.

% load everything
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

% we'll use the struct called 'combinedCompiledProjects_allEnhancers', wich
% contains one entry per nucleus.
% First, make a smaller struct called 'enhancerStruct'containing only the nc12 nuclei for the enhancer
% specified by the 'DataTye' argument.
for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

enhancerStruct = combinedCompiledProjects_allEnhancers(contains({combinedCompiledProjects_allEnhancers.dataSet},DataType)  &...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));

% bin nuclei 
%this function calculates the Dorsal fluorescence at some arbitrary time
%in nc12 and adds it to the struct in a 'DorsalFluoArbitraryTime' field
enhancerStruct = DorsalFluoArbitraryTime(enhancerStruct,fiducialTime);
nucleiFluorescence = [enhancerStruct.DorsalFluoArbitraryTime];

% ****IMPORTANT!!!**** this line uses the standard dorsal fluorescence, not the
% arbitrary one at a given 'fiducial time'
%nucleiFluorescence = [enhancerStruct.dorsalFluoFeature];



binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(enhancerStruct)
    enhancerStruct(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end

%plot([enhancerStruct.dorsalFluoFeature],[enhancerStruct.dorsalFluoBin2],'o')
coveredBins = unique([enhancerStruct.dorsalFluoBin2]);
binValues = binValues(coveredBins);

% now make one struct per bin
% define some filters
minEmbryosPerBin = 3;
minNucleiPerEmbryoPerBin = 1;

%store everything in these arrays
mean_maxFluo_acrossNuclei_perBin = [];
se_maxFluo_acrossNuclei_perBin = [];
mean_accFluo_acrossNuclei_perBin = [];
se_accFluo_acrossNuclei_perBin = [];
mean_fraction_acrossNuclei_perBin = [];
%
mean_maxFluo_acrossEmbryos_perBin = [];
se_maxFluo_acrossEmbryos_perBin = [];
mean_accFluo_acrossEmbryos_perBin = [];
se_accFluo_acrossEmbryos_perBin = [];
mean_fraction_acrossEmbryos_perBin = [];
se_fraction_acrossEmbryos_perBin = [];

OnNucleiPerBin = [];
OffNucleiPerBin = [];
for b = 1:length(coveredBins)
    binID = coveredBins(b);
    binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin2]== binID);
    activeNuc_Bin = length([binStruct.particleTimeOn]);
    inactiveNuc_Bin = length(binStruct) - activeNuc_Bin;
    % take the means and errors across nuclei first because it's easy      
    mean_maxFluo_acrossNuclei_perBin(b) = nanmean([binStruct.particleFluo95]);
    se_maxFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleFluo95])./sqrt(activeNuc_Bin);
    mean_accFluo_acrossNuclei_perBin(b) =  nanmean([binStruct.particleAccumulatedFluo]);
    se_accFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleAccumulatedFluo])./sqrt(activeNuc_Bin);
    mean_fraction_acrossNuclei_perBin(b) = activeNuc_Bin/length(binStruct);
    mean_timeOn_acrossNuclei_perBin(b) = nanmean([binStruct.particleTimeOn]);
    se_timeOn_acrossNuclei_perBin(b) = nanstd([binStruct.particleTimeOn])./sqrt(activeNuc_Bin);
    mean_duration_acrossNuclei_perBin(b) = nanmean([binStruct.particleDuration]);
    se_duration_acrossNuclei_perBin(b) = nanstd([binStruct.particleDuration])./sqrt(activeNuc_Bin);    
    totalRNA_acrossNuclei_perBin(b) = nansum([binStruct.particleAccumulatedFluo]);
    
    % save this for bootstrapping later
    OnNucleiPerBin(b) = activeNuc_Bin;
    OffNucleiPerBin(b) = inactiveNuc_Bin;
    
    % now deal with the mean and errors across embryos
    [uniqueEmbryos, ~, J]=unique({binStruct.prefix});
    occurences = histc(J, 1:numel(uniqueEmbryos));
    numEmbryos = length(occurences);
    qualityEmbryos = occurences >= minNucleiPerEmbryoPerBin;
    maxFluo_perEmbryo = [];
    accFluo_perEmbryo = [];
    fraction_perEmbryo = []; 
    nuclei_perEmbryo = [];
    timeOn_perEmbryo = [];
    duration_perEmbryo = [];
    totalRNA_perEmbryo = [];
    
    if numEmbryos >= minEmbryosPerBin     
        for e = 1:numEmbryos
            embryoPrefix = uniqueEmbryos{e};
%            if qualityEmbryos(e) %if this prefix has enough nuclei
            embryoStruct = binStruct(strcmpi({binStruct.prefix},embryoPrefix));
            nuclei_perEmbryo(e) = length(embryoStruct);
            maxFluo_perEmbryo(e) = nanmean([embryoStruct.particleFluo95]);
            accFluo_perEmbryo(e) =  nanmean([embryoStruct.particleAccumulatedFluo]);
            fraction_perEmbryo(e) = length([embryoStruct.particleTimeOn])/length(embryoStruct);
            timeOn_perEmbryo(e) = nanmean([embryoStruct.particleTimeOn]);
            duration_perEmbryo(e) = nanmean([embryoStruct.particleDuration]);
            totalRNA_perEmbryo(e) =  nansum([embryoStruct.particleAccumulatedFluo]).*fraction_perEmbryo(e);
%            end
        end
    end
    
    %add things taken across embryos to the arrays where we store stuff
    mean_maxFluo_acrossEmbryos_perBin(b) = nanmean(maxFluo_perEmbryo);
    se_maxFluo_acrossEmbryos_perBin(b) = nanstd(maxFluo_perEmbryo)./sqrt(numEmbryos);
    mean_accFluo_acrossEmbryos_perBin(b) = nanmean(accFluo_perEmbryo);
    se_accFluo_acrossEmbryos_perBin(b) = nanstd(accFluo_perEmbryo)./sqrt(numEmbryos);
    mean_fraction_acrossEmbryos_perBin(b) = nanmean(fraction_perEmbryo);
    se_fraction_acrossEmbryos_perBin(b) = std(fraction_perEmbryo)./sqrt(numEmbryos);
    mean_timeOn_acrossEmbryos_perBin(b) = nanmean(timeOn_perEmbryo);
    se_timeOn_acrossEmbryos_perBin(b) = nanstd(timeOn_perEmbryo)./sqrt(numEmbryos);
    mean_duration_acrossEmbryos_perBin(b) = nanmean(duration_perEmbryo);
    se_duration_acrossEmbryos_perBin(b) = nanstd(duration_perEmbryo)./sqrt(numEmbryos);
    mean_totalRNA_acrossEmbryos_perBin(b) = nanmean(totalRNA_perEmbryo);
    se__totalRNA_acrossEmbryos_perBin(b) = nanstd(totalRNA_perEmbryo)./sqrt(numEmbryos);
    
    %clear  nuclei_per_Embryo maxFluo_perEmbryo accFluo_perEmbryo fraction_perEmbryo timeOn_perEmbryo
     
end


%% bootstrap the fraction active across nuclei
cObsP = @(x) (sum(x)/length(x)); %bootstrapped function: fraction
nSamples = 1000;
btsrp_error_fraction_perBin = [];
for b = 1:length(coveredBins)
    NOnObs =  OnNucleiPerBin(b);
    NOffObs = OffNucleiPerBin(b);
    NOnSample = [ones(1,NOnObs) zeros(1,NOffObs)]; %this is the original sample
    if NOnObs>1
        [bootObsOn,~] = bootstrp(nSamples, cObsP,NOnSample); % this is the bootstrapped sample
    elseif NOffObs ==0
        bootObsOn = 1;
    else
        bootObsOn = 0;
    end
    btsrp_error_fraction_perBin(b) = std(bootObsOn);
end

%% Integrate the mean total mRNA produced per nucleus (inactive and active ones) across DV
% to show the total mRNA produced per embryo

embryoRNA = nansum(mean_totalRNA_acrossEmbryos_perBin);
embryoRNAError = sqrt(nansum(se__totalRNA_acrossEmbryos_perBin).^2);




%% make figures

if strcmpi(metric,'maxfluo') 
    if strcmpi(errorgroup,'nuclei')
    errorbar(ax,binValues,mean_maxFluo_acrossNuclei_perBin,se_maxFluo_acrossNuclei_perBin,'ro-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
    errorbar(ax,binValues,mean_maxFluo_acrossEmbryos_perBin,se_maxFluo_acrossEmbryos_perBin,'ko-','CapSize',0,'LineWidth',1.5,...
       'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('maximum spot fluorescence')
%set(gca,'XScale','log')
% ylim([0 600])
ylim([0,600])
xlim([0 3800])
%legend('across nuclei','across embryos')


elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_accFluo_acrossNuclei_perBin,se_accFluo_acrossNuclei_perBin,'o-','CapSize',0,'LineWidth',1.5,...
     'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
     elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_accFluo_acrossEmbryos_perBin,se_accFluo_acrossEmbryos_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('accumulated fluorescence')
%set(gca,'XScale','log')
xlim([0 3800])
ylim([0 1200])


elseif contains(lower(metric),'fraction')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_fraction_acrossNuclei_perBin,btsrp_error_fraction_perBin,'o-','CapSize',0,'LineWidth',1.5,...
            'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_fraction_acrossEmbryos_perBin,se_fraction_acrossEmbryos_perBin,'ko-','CapSize',0,'LineWidth',1.5,...
            'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('fraction of active nuclei')
%set(gca,'XScale','log')
ylim([0 1.1])
xlim([0 3500])


elseif contains(lower(metric),'timeon')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_timeOn_acrossNuclei_perBin,se_timeOn_acrossNuclei_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_timeOn_acrossEmbryos_perBin,se_timeOn_acrossEmbryos_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('turn on time')
%set(gca,'XScale','log')
ylim([0 10])
xlim([0 3800])

elseif contains(lower(metric),'total')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,totalRNA_acrossNuclei_perBin,se_timeOn_acrossNuclei_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_totalRNA_acrossEmbryos_perBin,se__totalRNA_acrossEmbryos_perBin,'o-','CapSize',0,'LineWidth',1.5,...
            'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('total produced mRNA')
%ylim([0 1800])
xlim([0 3800])
%set(gca,'YScale','log')

elseif contains(lower(metric),'duration')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_duration_acrossNuclei_perBin,se_duration_acrossNuclei_perBin,'o-','CapSize',0,'LineWidth',1.5,...
        'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_duration_acrossEmbryos_perBin,se_duration_acrossEmbryos_perBin,'o-','CapSize',0,'LineWidth',1.5,...
            'Color',Color,'MarkerFaceColor',Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
xlabel('Dorsal concentration (AU)')
ylabel('spot duration (min)')
ylim([0 6])
xlim([0 3800])
%set(gca,'YScale','log')

    

end
%set(gca,'Fontsize',20);



end 
