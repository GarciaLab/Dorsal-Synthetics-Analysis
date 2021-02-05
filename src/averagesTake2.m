function averagesTake2(DataType,numBins,metric,fiducialTime,errorgroup,ax)
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

enhancerStruct = combinedCompiledProjects_allEnhancers(contains({combinedCompiledProjects_allEnhancers.dataSet},DataType)  &...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));

% bin nuclei 
%fiducialTime = 3; %minutes
enhancerStruct = DorsalFluoArbitraryTime(enhancerStruct,fiducialTime);

nucleiFluorescence = [enhancerStruct.DorsalFluoArbitraryTime];
nucleiFluorescence = [enhancerStruct.dorsalFluoFeature];
binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(enhancerStruct)
    enhancerStruct(n).dorsalFluoBin = binnedNuclearFluo(n);
end

coveredBins = unique([enhancerStruct.dorsalFluoBin]);
binValues = binValues(coveredBins);

% now make one struct per bin
% define some filters
minEmbryosPerBin = 1;
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
    binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin]== binID);
    activeNuc_Bin = length([binStruct.compiledParticle]);
    inactiveNuc_Bin = length(binStruct) - activeNuc_Bin;
    % take the means and errors across nuclei first because it's easy      
    mean_maxFluo_acrossNuclei_perBin(b) = nanmean([binStruct.particleFluo95]);
    se_maxFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleFluo95])./sqrt(activeNuc_Bin);
    mean_accFluo_acrossNuclei_perBin(b) =  nanmean([binStruct.particleAccumulatedFluo]);
    se_accFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleAccumulatedFluo])./sqrt(activeNuc_Bin);
    mean_fraction_acrossNuclei_perBin(b) = activeNuc_Bin/length(binStruct);
    mean_timeOn_acrossNuclei_perBin(b) = nanmean([binStruct.particleTimeOn]);
    se_timeOn_acrossNuclei_perBin(b) = nanstd([binStruct.particleTimeOn])./sqrt(activeNuc_Bin);
    
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
    nuclei_per_Embryo = [];
    timeOn_perEmbryo = [];
    
    if numEmbryos >= minEmbryosPerBin     
        for e = 1:numEmbryos
            embryoPrefix = uniqueEmbryos{e};
%            if qualityEmbryos(e) %if this prefix has enough nuclei
            embryoStruct = binStruct(strcmpi({binStruct.prefix},embryoPrefix));
            nuclei_per_Embryo(e) = length(embryoStruct);
            maxFluo_perEmbryo(e) = nanmean([embryoStruct.particleFluo95]);
            accFluo_perEmbryo(e) =  nanmean([embryoStruct.particleAccumulatedFluo]);
            fraction_perEmbryo(e) = length([embryoStruct.compiledParticle])/length(embryoStruct);
            timeOn_perEmbryo(e) = nanmean([embryoStruct.particleTimeOn]);
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

%% make figures

if strcmpi(metric,'maxfluo') 
    if strcmpi(errorgroup,'nuclei')
    errorbar(ax,binValues,mean_maxFluo_acrossNuclei_perBin,se_maxFluo_acrossNuclei_perBin,'r-','CapSize',0,'LineWidth',1.5...
        ,'MarkerFaceColor','w','MarkerSize',3)
    elseif strcmpi(errorgroup,'embryos')
    errorbar(ax,binValues,mean_maxFluo_acrossEmbryos_perBin,se_maxFluo_acrossEmbryos_perBin,'ko-','CapSize',0,'LineWidth',1.5,...
       'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',3)
    end
xlabel('Dorsal concentration (AU)')
ylabel('maximum spot fluorescence')
%set(gca,'XScale','log')
% ylim([0 600])
ylim([100,600])
xlim([0 4000])
%legend('across nuclei','across embryos')


elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_accFluo_acrossNuclei_perBin,se_accFluo_acrossNuclei_perBin,'r-','CapSize',0,'LineWidth',1.5,...
     'MarkerFaceColor','w','MarkerSize',3)
     elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_accFluo_acrossEmbryos_perBin,se_accFluo_acrossEmbryos_perBin,'k-','CapSize',0,'LineWidth',1.5,...
        'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',3)
    end
xlabel('Dorsal concentration (AU)')
ylabel('accumulated fluorescence')
%set(gca,'XScale','log')
xlim([0 4000])
ylim([0 1200])


elseif contains(lower(metric),'fraction')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_fraction_acrossNuclei_perBin,btsrp_error_fraction_perBin,'r-','CapSize',0,'LineWidth',1.5,...
            'MarkerFaceColor','r','MarkerEdgeColor','none')
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_fraction_acrossEmbryos_perBin,se_fraction_acrossEmbryos_perBin,'k-','CapSize',0,'LineWidth',1.5,...
            'MarkerFaceColor','k','MarkerEdgeColor','none')
    end
xlabel('Dorsal concentration (AU)')
ylabel('fraction of active nuclei')
%set(gca,'XScale','log')
ylim([0 1])
xlim([0 4000])


elseif contains(lower(metric),'timeon')
    if strcmpi(errorgroup,'nuclei')
        errorbar(ax,binValues,mean_timeOn_acrossNuclei_perBin,se_timeOn_acrossNuclei_perBin,'ro-','CapSize',0,'LineWidth',1.5,...
        'MarkerFaceColor','w','MarkerSize',3)
    elseif strcmpi(errorgroup,'embryos')
        errorbar(ax,binValues,mean_timeOn_acrossEmbryos_perBin,se_timeOn_acrossEmbryos_perBin,'k-','CapSize',0,'LineWidth',1.5,...
            'MarkerFaceColor','none','MarkerEdgeColor','none','MarkerSize',3)
    end
xlabel('Dorsal concentration (AU)')
ylabel('turn on time')
%set(gca,'XScale','log')
ylim([0 10])
xlim([0 4000])


end
%set(gca,'Fontsize',20);



end 
