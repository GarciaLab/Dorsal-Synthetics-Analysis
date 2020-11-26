%% load stuff

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])%, 'dorsalResultsDatabase')
AllNC12Struct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);
% b = load('S:\Simon\Dropbox\DorsalSyntheticsDropbox\dorsalResultsDatabase.mat');

%Datasets = {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3','1DgVW'};
Datasets = {'1Dg11','TwiPE'}
PatserScores = [6.23,5.81,5.39,5.13,4.8,4.73,4.29];
% Nbins = 19; %number of Dorsal concentration bins
% NNuclei = 25; %maximum number of nuclei in nc12

% sort the struct according to nucleus fluorescence
AllNC12Table = struct2table(AllNC12Struct); % convert the struct to a table
sortedT = sortrows(AllNC12Table, 'dorsalFluoFeature'); % sort the table by dorsal fluo
sortedStruct = table2struct(sortedT); % change it back to struct array 
pseudoBinN = 15;
count=1;
figure
hold on
for e = 1:length(Datasets)
    
    EnhancerStruct = sortedStruct(contains({sortedStruct.dataSet}, Datasets{e})); %make a substruct with just this enhancer
    DorsalFluos = [];
    OutputmRNA = [];
    
    for i = 1:length(EnhancerStruct)
        Dorsal = EnhancerStruct(i).dorsalFluoFeature;
        mRNA = EnhancerStruct(i).particleFluo;
        if ~isempty(mRNA)
            %mRNA = mean(mRNA);
            mRNA = 1;
        else
            mRNA = 0; 
        end
        DorsalFluos(i) = Dorsal;
        OutputmRNA(i) = mRNA;
    end
    
    % correct for the number of nuclei
    DorsalFluos2=[];
    OutputmRNA2 = [];
    count=1;
    for i = 0:pseudoBinN:length(DorsalFluos)
        DorsalFluos2(count) = nanmean(DorsalFluos(i+1:min(i+pseudoBinN,length(DorsalFluos))));
        OutputmRNA2(count) = sum(OutputmRNA(i+1:min(i+pseudoBinN,length(DorsalFluos)))); 
        count=count+1;
    end
    
%     CumDist = cumsum(OutputmRNA);
    CumDist2 = cumsum(OutputmRNA2);
    %CumDist = CumDist./CumDist(end);
    CumDist2 = CumDist2./CumDist2(end);

    plot(DorsalFluos2,CumDist2,'-','LineWidth',2)
end
legend(Datasets)
hold off
ylabel('dNactive/d[Dl]')
xlabel('Dorsal fluorescence (AU)')


%% same but only consider nuclei with spots
AllNC12Table = struct2table(AllNC12Struct); % convert the struct to a table
sortedT = sortrows(AllNC12Table, 'dorsalFluoFeature'); % sort the table by dorsal fluo
sortedStruct = table2struct(sortedT); % change it back to struct array 
pseudoBinN = 15;
count=1;
figure
hold on
for e = 1:length(Datasets)   
    EnhancerStruct = sortedStruct(contains({sortedStruct.dataSet}, Datasets{e})); %make a substruct with just this enhancer
    DorsalFluos = [];
    OutputmRNA = [];
    count1=1;
    for i = 1:length(EnhancerStruct)
        Dorsal = EnhancerStruct(i).dorsalFluoFeature;
        mRNA = EnhancerStruct(i).particleFluo;
        if ~isempty(mRNA)
            OutputmRNA(count1) = mean(mRNA);
            DorsalFluos(count1) = Dorsal;
            count1=count1+1;
        end
    end
    
    
    % correct for the number of nuclei
    DorsalFluos2=[];
    OutputmRNA2 = [];
    count2=1;
    for i = 0:pseudoBinN:length(DorsalFluos)
        DorsalFluos2(count2) = nanmean(DorsalFluos(i+1:min(i+pseudoBinN,length(DorsalFluos))));
        OutputmRNA2(count2) = sum(OutputmRNA(i+1:min(i+pseudoBinN,length(DorsalFluos)))); 
        count2=count2+1;
    end
    
%     CumDist = cumsum(OutputmRNA);
    CumDist2 = cumsum(OutputmRNA2);
    %CumDist = CumDist./CumDist(end);
    %CumDist2 = CumDist2./CumDist2(end);

    plot(DorsalFluos2,CumDist2,'-','LineWidth',2)
end
legend(Datasets)
hold off
ylabel('dNactive/d[Dl]')
xlabel('Dorsal fluorescence (AU)')



%% Bin Dorsal such that there are linear increments in the accumulated mRNA
Nbins = 15;
FluoBinsValues = [0];
for c = linspace(1/Nbins,1,Nbins)
    % find position in the CDF array that is closest to this bin
    [~,idx] = min(abs(CumDist2-c));
    % find what Dl concentration this corresponds to
    FluoBinsValues = [FluoBinsValues DorsalFluos2(idx)];
end
FluoBinsValues

figure
plot(DorsalFluos2,CumDist2,'-ko','LineWidth',1)
hold on
for b = FluoBinsValues
    plot([b b],[0 1],'r-')
end
hold off

ylabel('integrated mRNA (normalized to max)')
xlabel('Dorsal fluorescence (AU)')
