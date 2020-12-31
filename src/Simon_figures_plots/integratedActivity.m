function integratedActivity(DataType,metric,ax)
% DataType should be an enhancer such as '1Dg11'
% metrics can be  maxfluo, accumulatedfluo or fraction
% ax is for plotting from the plotEVERYTHING funcion

% set up stuff
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat']);%, 'dorsalResultsDatabase')
% combinedCompiledProjects_allEnhancers is a struct with one nucleus per
% entry for all enhancers.
AllNC12Struct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);
AllNC12ThisEnhancer = AllNC12Struct(contains({AllNC12Struct.dataSet}, DataType));

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
timeOn = [];
isOn=[];
dorsalFluosForOnNuclei=[];
dorsalFluos = [sortedStruct.dorsalFluoFeature];

    
%loop over nuclei of this enhancer to populate the vectors
counter=1;
for i = 1:length(sortedStruct)    
    isOn(i) = ~isempty(sortedStruct(i).particleAccumulatedFluo);    
    if isOn(i)
        accFluo(counter) = sortedStruct(i).particleAccumulatedFluo;
        maxFluo(counter) = max(sortedStruct(i).particleFluo95);
        timeOn(counter) = sortedStruct(i).particleTimeOn;
        dorsalFluosForOnNuclei(counter) = dorsalFluos(i);
        counter = counter+1;
    end
end

% convert 0s to NaNs in maxFluo and accFluo so that we only take active
% nucleiS into account 
accFluo(accFluo==0)=nan;
maxFluo(maxFluo==0)=nan;

int_maxFluo = cumsum(maxFluo,'omitnan');
int_accFluo = cumsum(accFluo,'omitnan');
int_fraction = cumsum(isOn);
int_timeon = cumsum(timeOn);

% plot results
if strcmpi(metric,'maxfluo') 
plot(ax,dorsalFluosForOnNuclei,int_maxFluo,'b')
    
    
elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
plot(ax,dorsalFluosForOnNuclei,int_accFluo,'b')
    
    
elseif contains(lower(metric),'fraction')
plot(ax,dorsalFluos,int_fraction,'b')

elseif contains(lower(metric),'timeon')
plot(ax,dorsalFluosForOnNuclei,int_timeon,'b')
    
    
end









end




