function movingAverage(DataType,metric,averagingWindow,fiducialTime,ax)
% DataType should be an enhancer such as '1Dg11'
% metrics can be  maxfluo, accumulatedfluo or fraction
% averaging window is the number of nuclei
% ax is for plotting from the plotEVERYTHING funcion

% set up stuff
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat']);%, 'dorsalResultsDatabase')
% combinedCompiledProjects_allEnhancers is a struct with one nucleus per
% entry for all enhancers.
AllNC12Struct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);
AllNC12ThisEnhancer = AllNC12Struct(contains({AllNC12Struct.dataSet}, DataType));
AllNC12ThisEnhancer = DorsalFluoArbitraryTime(AllNC12ThisEnhancer,fiducialTime);


% sort the struct according to dorsal fluorescence
tempTable = struct2table(AllNC12ThisEnhancer); % convert the struct to a table
sortedT = sortrows(tempTable, 'dorsalFluoFeature'); % sort the table by dorsal fluo
sortedStruct = table2struct(sortedT); % change it back to struct array 

%get rid of nans
sortedStruct = sortedStruct(~isnan([sortedStruct.dorsalFluoFeature]));


% to store everything
maxFluo = [];
accFluo=[];
fraction=[];
timeon=[];
dorsalFluos = [sortedStruct.dorsalFluoFeature];

    
%loop over nuclei of this enhancer to populate the vectors
for i = 1:length(sortedStruct)    
    isOn(i) = ~isempty(sortedStruct(i).particleAccumulatedFluo);    
    if isOn(i)
        accFluo(i) = sortedStruct(i).particleAccumulatedFluo;
        maxFluo(i) = max(sortedStruct(i).particleFluo95);
        timeOn(i) = sortedStruct(i).particleTimeOn;
    else
        accFluo(i) = 0;
        maxFluo(i) = 0;
        timeOn(i) = nan;
    end

end

% convert 0s to NaNs in maxFluo and accFluo so that we only take active
% nucleiS into account 
accFluo(accFluo==0)=nan;
maxFluo(maxFluo==0)=nan;

% do a moving window with the same number of nuclei.
dorsalFluoWindow = [];
fractionOnWindow =[];
accFluoWindow = [];
maxFluoWindow = [];
timeOnWindow = [];

for i = 1:length(sortedStruct) 
    windowEnds = min(i+averagingWindow,length(sortedStruct));
    dorsalFluoWindow(i) = mean(dorsalFluos(i:windowEnds));
    fractionOnWindow(i) = mean(isOn(i:windowEnds));
    accFluoWindow(i) = nanmean(accFluo(i:windowEnds));
    maxFluoWindow(i) = nanmean(maxFluo(i:windowEnds));
    timeOnWindow(i) = nanmean(timeOn((i:windowEnds)));
end

% do a moving window with the same 'jump' in dorsal fluorescence every time
dorsalFluoWindow2 = [];
fractionOnWindow2 =[];
accFluoWindow2 = [];
maxFluoWindow2 = [];
timeOnWindow2 = [];
dorsalFluoBinWidth = 150;

for i = 1:max(dorsalFluos-(averagingWindow/3))
    windowStarts = i;
    windowEnds = i+dorsalFluoBinWidth;
    windowNuclei = logical((dorsalFluos>windowStarts) .* (dorsalFluos<windowEnds));
    dorsalFluoWindow2(i) = mean(dorsalFluos(windowNuclei));
    fractionOnWindow2(i) = mean(isOn(windowNuclei));
    accFluoWindow2(i) = nanmean(accFluo(windowNuclei));
    maxFluoWindow2(i) = nanmean(maxFluo(windowNuclei));
    timeOnWindow2(i) = nanmean(timeOn(windowNuclei));
end


fluoOnNuclei = dorsalFluos(isOn);
fluoOffNuclei = dorsalFluos(~isOn);

% plot results
if strcmpi(metric,'maxfluo') 
hold on
plot(ax,dorsalFluos,maxFluo,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
plot(ax,dorsalFluoWindow,maxFluoWindow,'b','LineWidth',2)
ylim([0 600])
    
    
elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
hold on
plot(ax,dorsalFluos,accFluo,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
plot(ax,dorsalFluoWindow,accFluoWindow,'b')
ylim([0 800])    
    
elseif contains(lower(metric),'fraction')
yyaxis left
hold on
histogram(fluoOffNuclei,'BinWidth',200,'FaceColor',[.4 .4 .4],'EdgeColor','none')
histogram(fluoOnNuclei,'BinWidth',200,'FaceColor',[.4 .4 1],'EdgeColor','none')
%plot(ax,dorsalFluos,isOn,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
yyaxis right
plot(ax,dorsalFluoWindow,fractionOnWindow,'b')
ylim([0 1.1])
%plot(dorsalFluoWindow2,fractionOnWindow2,'r')
%hold off

elseif contains(lower(metric),'timeon')
hold on
plot(ax,dorsalFluos,timeOn,'o','MarkerEdgeColor','none','MarkerFaceColor',[.7 .7 .7])
plot(ax,dorsalFluoWindow,timeOnWindow,'b','LineWidth',2) 
ylim([0 10])
    
end









end


