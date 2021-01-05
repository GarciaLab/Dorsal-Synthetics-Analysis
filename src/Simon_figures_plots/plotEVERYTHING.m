function plotEVERYTHING(metric)
% metric can be maxfluo, accfluo or fraction
%%
enhancers = {'1Dg-8D','1Dg11','1DgW','1Dg-5','1DgVW','1DgS2',...
    '1DgAW3','1DgVVW3','1Dg-12','1DgSVW2','2Dgc','TwiPEv5'};

enhancers = {'1Dg11_noExport','1Dg11_Export_first4min'}


%% the traditional way, binning
numBins = 20;
figure
tiledlayout('flow')
hold on
fiducialTime = 6;
for e = 1:length(enhancers)
    ax = nexttile;
    enhancerName = enhancers{e};
    averagesTake2(enhancerName,numBins,metric,fiducialTime,ax)
    legend(enhancers{e},'Location','northwest')
end
hold off


% %% moving average
% averagingWindow = 50; %number of nuclei
% figure
% tiledlayout('flow')
% hold on
% for e = 1:length(enhancers)
%     ax = nexttile;
%     enhancerName = enhancers{e};
%     movingAverage(enhancerName,metric,averagingWindow,ax)
%     legend(enhancers{e},'Location','northwest')
% end
% hold off



% %% accumulated
% figure
% tiledlayout('flow')
% hold on
% for e = 1:length(enhancers)
%     ax = nexttile;
%     enhancerName = enhancers{e};
%     integratedActivity(enhancerName,metric,ax)
%     legend(enhancers{e},'Location','northwest')
% end
% hold off


end