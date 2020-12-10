function averagesTake2(DataType,numBins)

% load everything
close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])

% we'll use the struct called 'combinedCompiledProjects_allEnhancers', wich
% contains one entry per nucleus.
% First, make a smaller struct called 'enhancerStruct'containing only the nc12 nuclei for the enhancer
% specified by the 'DataTye' argument.

enhancerStruct = combinedCompiledProjects_allEnhancers(strcmpi({combinedCompiledProjects_allEnhancers.dataSet},DataType)  &...
    [combinedCompiledProjects_allEnhancers.cycle]==12);

% define the bins and the bin the nuclei accordingly
dlfluobins = linspace(0,4500,numBins); 

for n = 1:length(enhancerStruct)
    nucleusFluo = enhnacerStruct(n).dorsalFluoFeature;
    

%% *** Maximum fluorescence
fluo95 = [];

fluo95_se = [];
for k = 1:length(dlfluobins)
    b = enhancerStruct([enhancerStruct.dorsalFluoBin]==k);
    fluo95(k) = nanmean([b.particleFluo95]);
    fluo95_se(k) = nanstd([b.particleFluo95])./sqrt(length([b.particleFluo95]));
    
end

plot(dlfluobins, fluo95);
errorbar(dlfluobins, fluo95, fluo95_se);


%% *** Fraction of active nuclei





%% *** Accumulated fluorescence



