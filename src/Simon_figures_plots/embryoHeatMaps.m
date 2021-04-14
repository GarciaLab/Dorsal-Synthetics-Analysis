function embryoHeatMaps(prefix,numBins)


%% load stuff
FrameInfo = getFrameInfo(LiveExperiment(prefix));

CompiledParticles = getCompiledParticles(LiveExperiment(prefix));
CompiledParticles = CompiledParticles.CompiledParticlesStruct;

Schnitzcells = getSchnitzcells(LiveExperiment(prefix));

for s = 1:length(Schnitzcells)
    Schnitzcells(s).originalSchnitzIdx = s;
end

NC12Schnitzcells = Schnitzcells([Schnitzcells.cycle]==12);
for i = 1:length(NC12Schnitzcells)
    if isempty(NC12Schnitzcells(i).FluoFeature)
        NC12Schnitzcells(i).FluoFeature = nan;
    end
end
NC12Schnitzcells2 = NC12Schnitzcells(~isnan([NC12Schnitzcells.FluoFeature]));

%% get the NC12 frames
NCs = [FrameInfo.nc];
AbsTime = [FrameInfo.Time];
NC12Start = find(NCs==12,1,'first');
NC12End = find(NCs==12,1,'last');
numFrames = NC12End-NC12Start;
AbsNC12Time = AbsTime(NC12Start:NC12End);
AbsNC12Time = AbsNC12Time(2:end)-NC12Start;
%% bin nuclei from scratch
nucleiFluorescence = [NC12Schnitzcells2.FluoFeature];
binValues = linspace(0,4500,numBins);
binnedNuclearFluo = BinData(nucleiFluorescence,binValues);
for n = 1:length(NC12Schnitzcells2)
    NC12Schnitzcells2(n).dorsalFluoBin2 = binnedNuclearFluo(n);
end
coveredBins = unique([NC12Schnitzcells2.dorsalFluoBin2]);
binValues = binValues(coveredBins);

%% go bin by bin and make a heatmap of the nuclei in that bin

for bi = coveredBins
    figure
    nucleiThisBin = find([NC12Schnitzcells2.dorsalFluoBin2]==bi);
    nucleiPerBin = sum([NC12Schnitzcells2.dorsalFluoBin2]==bi);
    binNucleiIdx = [NC12Schnitzcells2.originalSchnitzIdx];
    binNucleiIdx = binNucleiIdx(nucleiThisBin);
    AllNucleiArray = zeros(numFrames,nucleiPerBin);
    counter=1;
    for n = binNucleiIdx
        particleIdx = find([CompiledParticles.schnitz]==n);
        if particleIdx
            particleFluo = CompiledParticles(particleIdx).Fluo;
            particleFrames = CompiledParticles(particleIdx).Frame;
            particleFrames = particleFrames-NC12Start;
            AllNucleiArray(particleFrames,counter) = particleFluo;
        else
            AllNucleiArray(:,counter) = 0;
        end
        counter = counter+1;
    end
    imagesc(AllNucleiArray')
    title(['bin = ' num2str(bi)])
end


%% pool together all bins (all nucleaise)

figure
%nucleiThisBin = find([NC12Schnitzcells2.dorsalFluoBin2]==bi);
%nucleiPerBin = sum([NC12Schnitzcells2.dorsalFluoBin2]==bi);
%binNucleiIdx = [NC12Schnitzcells2.originalSchnitzIdx];
%binNucleiIdx = binNucleiIdx(nucleiThisBin);
AllNucleiArray = zeros(length(NC12Schnitzcells2),numFrames);
for n = 1:length(NC12Schnitzcells2)
    SchnitzID = NC12Schnitzcells2(n).originalSchnitzIdx;
    particleIdx = find([CompiledParticles.schnitz]==SchnitzID);
    if particleIdx
        particleFluo = CompiledParticles(particleIdx).Fluo;
        particleFrames = CompiledParticles(particleIdx).Frame;
        particleFrames = particleFrames-NC12Start;
        AllNucleiArray(n,particleFrames) = particleFluo;
    else
        AllNucleiArray(n,:) = 0;
    end
end
AllNucleiArray(AllNucleiArray<0)=0;

%% sort nuclei by their particle fluorescence
sumFluo = sum(AllNucleiArray,2);
AllNucleiArray(:,end) = sumFluo;
SortedAllNucleiArray = sortrows(AllNucleiArray,size(AllNucleiArray,2));
SortedAllNucleiArray(:,end)=0;
imagesc(flip(SortedAllNucleiArray))
title(['All nc12 nuclei from ' prefix])
xlabel('time into nc12 (min)')
ylabel('nuclei')
%% deal with labels
% timeLabels = 0:5:15;
% find(AbsNC12Time./60==
% labels = AbsNC12Time(xticks)./60;
% xticklabels(num2str(labels))
%% 








end
