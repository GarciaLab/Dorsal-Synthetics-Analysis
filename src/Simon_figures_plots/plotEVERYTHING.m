function plotEVERYTHING(metric,method)
% metric can be maxfluo, accfluo or fraction
%%
allenhancers = {'1Dg-8D','1Dg11','1DgW','1Dg-5','1DgVW','1DgS2',...
    '1DgAW3','1DgVVW3','1Dg-12','1DgSVW2','2Dgc','TwiPEv5'};
optoenhancers = {'1Dg11_noExport','1Dg11_Export_first4min','1Dg11_exportedAfter4min'};

affinity_enhancers = {'1Dg11_2xDl','1DgS2','1DgW_2x','1DgAW3','1DgSVW2','1DgVVW3','1DgVW'};
%Palette = brewermap(length(affinity_enhancers),'YlGnBu');
affinity_enhancers = {'1Dg11_2xDl_FFF'}%,'1DgS2','1DgW_2xDl_FFF','1DgAW3','1DgSVW2','1DgVVW3','1DG_VW_2xDl_FFF'};
%affinity_enhancers = {'1Dg11_2xDl','1Dg11_FFF'};
Palette = viridis(length(affinity_enhancers));
paperNames = {'6.23','5.81','5.39','5.13','4.8','4.73','4.29'};
scores = [6.23,5.81,5.39,5.13,4.8,4.73,4.29];
scores = scores(1:length(affinity_enhancers));
enhancers = affinity_enhancers;

%% the traditional way, binning
if contains(lower(method),'binning')

    numBins = 20;
    embryoRNAs = [];
    embryoRNAErrors= [];
    errorgroup = 'embryos'; %whether error is calculated across nuclei or across embryos
    fig = figure;
    ax = axes(fig);
    %tiledlayout(1,length(enhancers), 'TileSpacing', 'compact', 'Padding', 'compact')
    % tiledlayout('flow')
    hold on
    fiducialTime = 6;
    for e = 1:length(enhancers)
        Color = Palette(e,:);
        %ax = nexttile; %comment this out to plot everything in one graph
        enhancerName = enhancers{e};
        [embryoRNAs(e) embryoRNAErrors(e)] = averagesTake2(enhancerName,numBins,metric,fiducialTime,errorgroup,Color,ax);
        legend(paperNames{e},'Location','northwest')
    end
    hold off
    figure
    errorbar(scores,embryoRNAs,embryoRNAErrors,'ro-','CapSize',0,'MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','none')
    ylabel('mRNA produced per embryo')
    xlabel('binding site affinity (patser score)')


%% moving average
elseif contains(lower(method),'movingaverage')
    embryoMeanOverAll = [];
    embryoStdOverAll = [];
    
    averagingWindow = 30; %number of nuclei
    figure
    tiledlayout('flow')
    hold on
    fiducialTime = 6;
    for e = 1:length(enhancers)
        Color = Palette(e,:);
        ax = nexttile;
        enhancerName = enhancers{e};
        [embryoMeanOverAll(e) embryoStdOverAll(e)] = movingAverage(enhancerName,metric,averagingWindow,fiducialTime,Color,ax);
        legend(paperNames{e},'Location','northwest')
    end
    hold off
    
    figure
    errorbar(scores,embryoMeanOverAll,embryoStdOverAll,'r-o','CapSize',0,'MarkerFaceColor','r','MarkerEdgeColor','none','MarkerSize',8)
    title(metric)


%% accumulated
elseif contains(lower(method),'cumulative')
    %window = 20;
    numBins = 40;
    fiducialTime = 6; %minutes into nc12
    fig = figure;
    ax = axes(fig);
    tiledlayout('flow')
    hold on
    for e = 1:length(enhancers)
        Color = Palette(e,:);
        %ax = nexttile; %comment this out to plot everything in one graph
        enhancerName = enhancers{e};
        integratedActivity2(enhancerName,metric,fiducialTime,numBins,Color,ax)
        legend(paperNames{e},'Location','northwest')
    end
    hold off
end

end