function integratedActivity2(DataType,metric,window,ax)
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


% bin nuclei according to Dorsal-Venus fluorescence
DorsalFluos = [AllNC12ThisEnhancer.dorsalFluoFeature];
binLimits = (0:100:nanmax(DorsalFluos));
BinnedDorsalFluos = BinData(DorsalFluos,binLimits);
BinnedDorsalFluos(isnan(DorsalFluos)) = nan;

for i = 1:length(AllNC12ThisEnhancer)
    AllNC12ThisEnhancer(i).dorsalFluoBin2 = BinnedDorsalFluos(i);
end

% now go over bins and calculate the fraction of active nuclei

FractionPerBin = [];
MaxFluoPerBin = [];
AccumulatedFluoPerBin = [];
TimeOnPerBin = [];
for bin = 1:max(BinnedDorsalFluos)
    BinOnNuclei = 0;
    BinOffNuclei = 0;
    BinMaxFluo = nan;
    BinAccumulatedFluo = nan;
    BinTimeOn = nan;
    
    for n = 1:length(AllNC12ThisEnhancer)
        if AllNC12ThisEnhancer(n).dorsalFluoBin2 == bin
            if ~isempty(AllNC12ThisEnhancer(n).particleFluo)
                BinOnNuclei = BinOnNuclei+1;
                BinMaxFluo = [BinMaxFluo AllNC12ThisEnhancer(n).particleFluo95];
                BinAccumulatedFluo = [BinAccumulatedFluo AllNC12ThisEnhancer(n).particleAccumulatedFluo];
                BinTimeOn = [BinTimeOn AllNC12ThisEnhancer(n).particleTimeOn];
            else
                BinOffNuclei = BinOffNuclei+1;
            end
        end
    end
    
    FractionPerBin = [FractionPerBin BinOnNuclei/(BinOnNuclei+BinOffNuclei)];
    MaxFluoPerBin = [MaxFluoPerBin nanmean(BinMaxFluo)];
    AccumulatedFluoPerBin = [AccumulatedFluoPerBin nanmean(BinAccumulatedFluo)];
    TimeOnPerBin = [TimeOnPerBin nanmean(BinTimeOn)];
end


integrated_fraction = cumsum(FractionPerBin,'omitnan');
integrated_maxfluo = cumsum(MaxFluoPerBin,'omitnan');
integrated_accumulatedfluo = cumsum(AccumulatedFluoPerBin,'omitnan');
integrated_timeon = cumsum(TimeOnPerBin,'omitnan');



if contains(lower(metric),'fraction')
plot(binLimits,integrated_fraction,'b','LineWidth',2)
ylabel('fraction active')

elseif strcmpi(metric,'maxfluo') 
plot(binLimits,integrated_maxfluo,'b','LineWidth',2)
ylabel('maximum fluo')

elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
plot(binLimits,integrated_accumulatedfluo,'b','LineWidth',2)
ylabel('accumulated fluo')

elseif contains(lower(metric),'timeon')
plot(binLimits,integrated_timeon,'b','LineWidth',2)
ylabel('time on')


% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % sort the struct according to dorsal fluorescence
% tempTable = struct2table(AllNC12ThisEnhancer); % convert the struct to a table
% sortedT = sortrows(tempTable, 'dorsalFluoFeature'); % sort the table by dorsal fluo
% sortedStruct = table2struct(sortedT); % change it back to struct array 
% 
% %get rid of nans
% sortedStruct = sortedStruct(~isnan([sortedStruct.dorsalFluoFeature]));
% 
% 
% % to store everything
% maxFluo = [];
% accFluo=[];
% fraction=[];
% timeOn = [];
% isOn=[];
% dorsalFluosForOnNuclei=[];
% dorsalFluos = [sortedStruct.dorsalFluoFeature];
% 
%     
% %loop over nuclei of this enhancer to populate the vectors
% counter=1;
% for i = 1:length(sortedStruct)    
%     isOn(i) = ~isempty(sortedStruct(i).particleAccumulatedFluo);    
%     if isOn(i)
%         accFluo(counter) = sortedStruct(i).particleAccumulatedFluo;
%         maxFluo(counter) = max(sortedStruct(i).particleFluo95);
%         timeOn(counter) = sortedStruct(i).particleTimeOn;
%         dorsalFluosForOnNuclei(counter) = dorsalFluos(i);
%         counter = counter+1;
%     end
% end
% 
% % convert 0s to NaNs in maxFluo and accFluo so that we only take active
% % nucleiS into account 
% accFluo(accFluo==0)=nan;
% maxFluo(maxFluo==0)=nan;
% 
% int_maxFluo = cumsum(maxFluo,'omitnan');
% int_accFluo = cumsum(accFluo,'omitnan');
% int_fraction = cumsum(isOn);
% int_timeon = cumsum(timeOn);
% 
% 
% % we want the integration window to be always the same number of nuclei
% 
% %deal with the active nuclei only first
% int_maxFluo2=maxFluo(1);
% int_accFluo2=accFluo(1);
% int_timeon2=timeOn(1);
% dorsalFluosForOnNuclei2=dorsalFluosForOnNuclei(1);
% for j = 0:window:length(maxFluo)-window
%     int_maxFluo2 = [int_maxFluo2 int_maxFluo2(end)+sum(maxFluo(j+1:j+window))];
%     int_accFluo2 = [int_accFluo2 int_accFluo2(end)+sum(accFluo(j+1:j+window))];
%     int_timeon2 = [int_timeon2 int_timeon2(end)+sum(timeOn(j+1:j+window))];
%     dorsalFluosForOnNuclei2 = [dorsalFluosForOnNuclei2 mean(dorsalFluosForOnNuclei(j+1:j+window))];
% end
% %deal with all nuclei now
% int_fraction2=isOn(1);
% dorsalFluosForAll2=dorsalFluos(1);
% for j = 0:window:length(dorsalFluos)-window
%     int_fraction2 = [int_fraction2 int_fraction2(end)+sum(isOn(j+1:j+window))];
%     dorsalFluosForAll2 = [dorsalFluosForAll2 mean(dorsalFluos(j+1:j+window))];
% end
% 
% 
% 
% % plot results
% if strcmpi(metric,'maxfluo') 
% plot(ax,dorsalFluosForOnNuclei2,int_maxFluo2,'b','LineWidth',2)
%     
% elseif strcmpi(metric,'accumulatedfluo') || strcmpi(metric,'accfluo')
% plot(ax,dorsalFluosForOnNuclei2,int_accFluo2,'b''LineWidth',2)
%     
% elseif contains(lower(metric),'fraction')
% plot(ax,dorsalFluosForAll2,int_fraction2,'b','LineWidth',2)
% 
% elseif contains(lower(metric),'timeon')
% plot(ax,dorsalFluosForOnNuclei2,int_timeon2,'b','LineWidth',2)
%     
%     
% end
% 
% 







end




