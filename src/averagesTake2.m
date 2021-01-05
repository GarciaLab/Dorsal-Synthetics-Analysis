function averagesTake2(DataType,numBins,metric,fiducialTime,ax)
% metric can be maxfluo, accumulatedfluo or fraction

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


for b = 1:length(coveredBins)
    binID = coveredBins(b);
    binStruct = enhancerStruct([enhancerStruct.dorsalFluoBin]== binID);
    activeNuc_thisBin = length([binStruct.compiledParticle]);
    % take the means and errors across nuclei first because it's easy      
    mean_maxFluo_acrossNuclei_perBin(b) = nanmean([binStruct.particleFluo95]);
    se_maxFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleFluo95])./(activeNuc_thisBin-1);
    mean_accFluo_acrossNuclei_perBin(b) =  nanmean([binStruct.particleAccumulatedFluo]);
    se_accFluo_acrossNuclei_perBin(b) = nanstd([binStruct.particleAccumulatedFluo])./(length(binStruct)-1);
    mean_fraction_acrossNuclei_perBin(b) = activeNuc_thisBin/length(binStruct);
    %
      
    % now deal with the mean and errors across embryos
    [uniqueEmbryos, ~, J]=unique({binStruct.prefix});
    occurences = histc(J, 1:numel(uniqueEmbryos));
    numEmbryos = length(occurences);
    qualityEmbryos = occurences >= minNucleiPerEmbryoPerBin;
    maxFluo_perEmbryo = [];
    accFluo_perEmbryo = [];
    fraction_perEmbryo = []; 
    nuclei_per_Embryo = [];
    
    if numEmbryos >= minEmbryosPerBin     
        for e = 1:numEmbryos
            embryoPrefix = uniqueEmbryos{e};
            if qualityEmbryos(e) %if this prefix has enough nuclei
                embryoStruct = binStruct(strcmpi({binStruct.prefix},embryoPrefix));
                nuclei_per_Embryo(e) = length(embryoStruct);
                maxFluo_perEmbryo(e) = nanmean([embryoStruct.particleFluo95]);
                accFluo_perEmbryo(e) =  nanmean([embryoStruct.particleAccumulatedFluo]);
                fraction_perEmbryo(e) = length([embryoStruct.compiledParticle])/length(embryoStruct);
            end
        end
    end
    
    %add things taken across embryos to the arrays where we store stuff
    mean_maxFluo_acrossEmbryos_perBin(b) = nanmean(maxFluo_perEmbryo);
    se_maxFluo_acrossEmbryos_perBin(b) = nanstd(maxFluo_perEmbryo)./(numEmbryos-1);
    mean_accFluo_acrossEmbryos_perBin(b) = nanmean(accFluo_perEmbryo);
    se_accFluo_acrossEmbryos_perBin(b) = nanstd(accFluo_perEmbryo)./(numEmbryos-1);
    mean_fraction_acrossEmbryos_perBin(b) = nanmean(fraction_perEmbryo);
    se_fraction_acrossEmbryos_perBin(b) = std(fraction_perEmbryo)./(numEmbryos-1);
         
end

% make figures



if strcmpi(metric,'maxfluo') 
% figure %max fluo
% hold on
errorbar(binValues,mean_maxFluo_acrossNuclei_perBin,se_maxFluo_acrossNuclei_perBin,'r-','CapSize',0,'LineWidth',1.5...
    ,'MarkerFaceColor','w','MarkerSize',3)
%errorbar(ax,binValues,mean_maxFluo_acrossEmbryos_perBin,se_maxFluo_acrossEmbryos_perBin,'bo-','CapSize',0,'LineWidth',1.5,...
%    'MarkerFaceColor','w','MarkerSize',3)
%hold off
xlabel('Dorsal concentration (AU)')
ylabel('maximum spot fluorescence')
%set(gca,'XScale','log')
% ylim([0 600])
ylim([100,400])
xlim([0 4000])
%legend('across nuclei','across embryos')

elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
% figure %accumulated fluo
% hold on
errorbar(binValues,mean_accFluo_acrossNuclei_perBin,se_accFluo_acrossNuclei_perBin,'r-','CapSize',0,'LineWidth',1.5,...
     'MarkerFaceColor','w','MarkerSize',3)
%errorbar(ax,binValues,mean_accFluo_acrossEmbryos_perBin,se_accFluo_acrossEmbryos_perBin,'ro-','CapSize',0,'LineWidth',1.5,...
%    'MarkerFaceColor','w','MarkerSize',3)
%hold off
xlabel('Dorsal concentration (AU)')
ylabel('accumulated fluorescence')
%set(gca,'XScale','log')
xlim([0 4000])
ylim([0 1200])
%legend('across nuclei','across embryos')

elseif contains(lower(metric),'fraction')
% % figure;
% % hold on
%plot(binValues,mean_fraction_acrossNuclei_perBin,'ro-','LineWidth',1.5)
errorbar(ax,binValues,mean_fraction_acrossEmbryos_perBin,se_fraction_acrossEmbryos_perBin,'ko-','CapSize',0,'LineWidth',1.5)
%hold off
xlabel('Dorsal concentration (AU)')
ylabel('fraction of active nuclei')
%set(gca,'XScale','log')
%legend('across nuclei','across embryos')
ylim([0 1])
xlim([0 4000])

end
%set(gca,'Fontsize',20);



end 
