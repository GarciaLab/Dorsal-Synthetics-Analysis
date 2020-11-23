function checkSchnitzAssignmentToParticles(Prefix)

% [~,~,DropboxFolder,~, PreProcPath,...
%     ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);
% resultsFolder = [DropboxFolder, filesep, Prefix];

liveExperiment = LiveExperiment(Prefix);
resultsFolder = liveExperiment.resultsFolder;

load([resultsFolder,'CompiledParticles.mat'],'CompiledParticles');
load([resultsFolder, 'FrameInfo.mat']);
load([resultsFolder, Prefix, '_lin.mat'])

ch = 1;
CompiledParticles = CompiledParticles{ch};

%%

if ~isempty(CompiledParticles)
    
    particleNCs = [CompiledParticles.nc];
    nc12Idx = [1:length(CompiledParticles)] .* (particleNCs==12);
    
    minDistance = 50; %pixels of distance between a spot and a nucleus centroid to be matched
    for p = 1:nc12Idx
        
        particleXPosPerFrame = CompiledParticles(p).xPos;
        particleYPosPerFrame = CompiledParticles(p).yPos;
        particleFrames = CompiledParticles(p).Frame;
        closestNucleusPerFrame = nan(1,length(particleFrames));
        
        for f = 1:length(particleFrames)
            
            frame = particleFrames(f);
            particleXPosThisFrame = particleXPosPerFrame(f);
            particleYPosThisFrame = particleYPosPerFrame(f);
            distanceToNucleiThisFrame = nan(1,length(schnitzcells));
            %         nucleiThisFrameIDs = [];
            
            for n = 1:length(schnitzcells)
                schnitzFrames = schnitzcells(n).frames;
                if ismember(frame,schnitzFrames)
                    schnitzXPos = schnitzcells(n).cenx;
                    schnitzXPos = schnitzXPos(schnitzFrames==frame);
                    schnitzYPos = schnitzcells(n).ceny;
                    schnitzYPos = schnitzYPos(schnitzFrames==frame);
                    distanceToNucleiThisFrame(n) = sqrt((particleXPosThisFrame-schnitzXPos)^2+(particleYPosThisFrame-schnitzYPos)^2);
                end
            end
            [dist,idx] = min(distanceToNucleiThisFrame);
            if dist < minDistance
                closestNucleusPerFrame(f) = idx;
            end
        end
        closestNucleusEver = mode(closestNucleusPerFrame);
        if ~isnan(closestNucleusEver)
            assert(CompiledParticles(p).schnitz == closestNucleusEver,'assigned schnitz is not the right one')
        end
        
    end
    
    display('congrats! schnitzes were correctly assigned to particles')
    %%
    
end

end

