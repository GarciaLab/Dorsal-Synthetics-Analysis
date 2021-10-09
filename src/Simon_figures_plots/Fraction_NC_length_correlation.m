function Fraction_NC_length_correlation(numBins)

% we wish to know if there is a positive relationship between the length of
% the nuclear cycle and the fraction of nuclei that get a chance to
% transcribe



%% LOAD STUFF

% load everything
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])


%% FILTER DATA
% we'll use the struct called 'combinedCompiledProjects_allEnhancers', wich
% contains one entry per nucleus.
% First, make a smaller struct called 'enhancerStruct'containing only the nc12 nuclei for the enhancer
% specified by the 'DataTye' argument.
% some QC filtering
for i = 1:length(combinedCompiledProjects_allEnhancers)
    if isempty(combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature)
        combinedCompiledProjects_allEnhancers(i).dorsalFluoFeature = nan;
    end
end

dataTypeExcludeRule = ~contains(lower({combinedCompiledProjects_allEnhancers.dataSet}),'export');

filteredStruct = combinedCompiledProjects_allEnhancers(dataTypeExcludeRule & ...
    [combinedCompiledProjects_allEnhancers.cycle]==12 & ~isnan([combinedCompiledProjects_allEnhancers.dorsalFluoFeature]));

nucleiFluorescence = [filteredStruct.dorsalFluoFeature];

%% GET NC12 duration per dataset and add tu struct
[filteredStruct(:).nc12duration] = deal(nan);

Datasets = unique({filteredStruct.prefix});
for p = 1:length(Datasets)
    p/length(Datasets)
    prefix = Datasets{p}
    MitosisFrames = LiveExperiment(prefix).anaphaseFrames;
    FrameInfo = getFrameInfo(LiveExperiment(prefix));
    NC12start = MitosisFrames(end-2);
    NC12end = MitosisFrames(end-1);
    % calculate the duration of nc 12
    if NC12start ~= 0 && ~isnan(NC12end)
        NC12duration = FrameInfo(NC12end).Time - FrameInfo(NC12start).Time;
    else
        NC12duration = NaN;
    end    
    % add the duration to 'filteredStruct'
    [filteredStruct(strcmpi({filteredStruct.prefix},prefix)).nc12duration] = deal(NC12duration);
end

% keep only datasets where the start of NC12 was determined
hasNC12 = ~isnan([filteredStruct.nc12duration]);
filteredStruct = filteredStruct(hasNC12);
clearvars -except filteredStruct
%% bin nuclei according to Dl fluorescence
nucleiFluorescence = [filteredStruct.dorsalFluoFeature];
numBins = 19; % THIS SHOULD BE A FUNCTION VARIABLE!
binValues = linspace(0,3800,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(filteredStruct)
    filteredStruct(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end

coveredBins = unique([filteredStruct.dorsalFluoBin2]);
binValues = binValues(coveredBins);

% another way of binning
dorsalVals = linspace(0,3800,numBins);
binSize = diff(dorsalVals(1:2));
for n = 1:length(filteredStruct)
    FluoFeature = filteredStruct(n).dorsalFluoFeature;
    bin = find(FluoFeature-dorsalVals<0,1,'first')-1;
    filteredStruct(n).dorsalFluoBin3 = bin;
end
coveredBins = unique([filteredStruct.dorsalFluoBin3]);
binValues = dorsalVals(1:end-1) + (binSize/2);
clearvars -except filteredStruct
%% for each enhancer, each Dl fluo bin, plot fraction active of an embryo vs its NC12 duration
close all
dataTypes = unique({filteredStruct.dataSet}); % what enhancer
Correlations = [];
Palette = jet(length(dataTypes));
forLegend = {};
figure
hold on
for e = 1:length(dataTypes)
    Color = Palette(e,:);
    enhancer = dataTypes{e};
    forLegend{e} = enhancer;
    dataTypeIncludeRule = contains({filteredStruct.dataSet},enhancer);
    enhancerStruct = filteredStruct(dataTypeIncludeRule);  
    coveredBins = unique([enhancerStruct.dorsalFluoBin2]);
    durationDataForPlot = [];
    fractionDataForPlot = [];
    RvaluesForLegened = [];
    
    for b = coveredBins
        enhancerBinStruct = enhancerStruct([enhancerStruct.dorsalFluoBin2] == b);       
        coveredEmbryos = unique({enhancerBinStruct.prefix});
        fractionActivePerEmbryo = [];
        nc12DurationPerEmbryo = [];
        
        for p = 1:length(coveredEmbryos) % these are the embryos in this Dl bin of this enhnacer
            prefix = coveredEmbryos{p};
            enhancerBinEmbryoStruct = enhancerBinStruct(strcmpi({enhancerBinStruct.prefix},prefix));
            numParticles = length([enhancerBinEmbryoStruct.particleTimeOn]);
            numNuclei = length(enhancerBinEmbryoStruct);
            nc12DurationPerEmbryo = [nc12DurationPerEmbryo enhancerBinEmbryoStruct(1).nc12duration];
%             durationDataForPlot = [durationDataForPlot nc12DurationPerEmbryo];
            
            if numParticles == 0
                fractionActivePerEmbryo = [fractionActivePerEmbryo 0];
            else
                fractionActivePerEmbryo = [fractionActivePerEmbryo numParticles/numNuclei];
            end
%             fractionDataForPlot = [fractionDataForPlot fractionActivePerEmbryo];
        end
        fractionActivePerEmbryo = fractionActivePerEmbryo./mean(fractionActivePerEmbryo);

        if length(coveredEmbryos)>2
            R = corrcoef(nc12DurationPerEmbryo,fractionActivePerEmbryo); R = R(1,2);
            Correlations = [Correlations R];
            RvaluesForLegened - [RvaluesForLegened R];
            fractionDataForPlot = [fractionDataForPlot fractionActivePerEmbryo];
            durationDataForPlot = [durationDataForPlot nc12DurationPerEmbryo];
        end
        %plot(nc12DurationPerEmbryo,fractionActivePerEmbryo,'o','MarkerFaceColor',Color,'MarkerEdgeColor','none')
    end
     plot(durationDataForPlot./60,fractionDataForPlot,'ko','MarkerFaceColor','k','MarkerEdgeColor','none')
end
hold off
legend(forLegend)
ylabel({'fraction active','normalized by the mean of bin and enhancer'})
xlabel('embryo nc 12 duration (minutes)')


figure
histogram(Correlations)
xlabel({'Pearson correlation between','NC length and fraction active'})
ylabel('counts')


%             
%             
%             
%             durationDataForPlot = [durationDataForPlot nc12DurationPerEmbryo]
%             fractionDataForPlot = [fractionDataForPlot fractionActivePerEmbryo]
%             
%             if length(nc12DurationPerEmbryo) > 2
%                 R = corrcoef(nc12DurationPerEmbryo,fractionActivePerEmbryo); R = R(1,2);
%                 durationDataForPlot = [durationDataForPlot nc12DurationPerEmbryo];
%                 fractionDataForPlot = [fractionDataForPlot fractionActivePerEmbryo];
%                 Correlations = [Correlations R];
%                 plot(nc12DurationPerEmbryo./60,fractionActivePerEmbryo,...
%                     'o','MarkerFaceColor',Color,'MarkerEdgeColor','none')
%             end
%         end
%     end
% %     plot(durationDataForPlot,fractionDataForPlot,'o','MarkerFaceColor',Color,'MarkerEdgeColor','none')
% end
% hold off
% xlabel('NC12 duration (min)')
% ylabel({'fraction active','normalized to mean of this enhancer and Dl fluo bin'})


% figure
% histogram(Correlations)
% xlabel({'Pearson correlation between','NC length and fraction active'})
% ylabel('counts')









% 
% 
% %% load data
% AllDataPath = 'S:\Simon\Dropbox\DorsalSyntheticsDropbox\';
% load(AllDataPath + "dorsalResultsDatabase.mat");
% 
% 
% %%
% prefix = '2020-07-27-1DgVW_2xDl_5'
% MitosisFrames = LiveExperiment(prefix).anaphaseFrames
% FrameInfo = getFrameInfo(LiveExperiment(prefix))
