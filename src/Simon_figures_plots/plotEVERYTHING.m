function plotEVERYTHING(metric)
% metric can be maxfluo, accfluo or fraction
%%
allenhancers = {'1Dg-8D','1Dg11','1DgW','1Dg-5','1DgVW','1DgS2',...
    '1DgAW3','1DgVVW3','1Dg-12','1DgSVW2','2Dgc','TwiPEv5'};
optoenhancers = {'1Dg11_noExport','1Dg11_Export_first4min','1Dg11_exportedAfter4min'};
affinity_enhancers = {'1Dg11_2xDl_FFF','1DgW_2x','1DgVW','1DgS2','1DgAW3','1DgVVW3','1DgSVW2'};

enhancers = affinity_enhancers;

%% the traditional way, binning
numBins = 18;
errorgroup = 'embryos'; %whether error is calculated across nuclei or across embryos
figure
tiledlayout('flow')
hold on
fiducialTime = 6;
for e = 1:length(enhancers)
    ax = nexttile;
    enhancerName = enhancers{e};
    averagesTake2(enhancerName,numBins,metric,fiducialTime,errorgroup,ax)
    legend(enhancers{e},'Location','northwest')
end
hold off


% %% moving average
% averagingWindow = 40; %number of nuclei
% figure
% tiledlayout('flow')
% hold on
% fiducialTime = 6;
% for e = 1:length(enhancers)
%     ax = nexttile;
%     enhancerName = enhancers{e};
%     movingAverage(enhancerName,metric,averagingWindow,fiducialTime,ax)
%     legend(enhancers{e},'Location','northwest')
% end
% hold off



% %% accumulated
% window = 20;
% figure
% tiledlayout('flow')
% hold on
% for e = 1:length(enhancers)
%     ax = nexttile;
%     enhancerName = enhancers{e};
%     integratedActivity2(enhancerName,metric,window,ax)
%     legend(enhancers{e},'Location','northwest')
% end
% hold off


end