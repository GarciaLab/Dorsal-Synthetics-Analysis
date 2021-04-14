
% DESCRIPTION
% Sorts particle traces into transcriptional and non-transcriptional
% groups, then calculates the mean and 95th percentile intensity for
% transcriptional and non-transcriptional events. Intended for datasets
% with a DNA label channel and a transcription channel, where
% transcriptional measurements are made in active and inactive nuclei.
% 
% ARGUMENTS
% prefix1: dataset segmented using DNA label channel
% prefix2: dataset segmented using transcription channel
%
% OUTPUT
% 'max_intensities_non': single array containing 95th percentile intensities
% for non-transcriptional events
%
% 'max_intensities_transcription': single array containing 95th percentile intensities
% for transcriptional events
%
% 'mean_intensities_non': single array containing mean intensities
% for non-transcriptional events
%
% 'mean_intensities_transcription': single array containing mean intensities
% for transcriptional events

function [max_intensities_non,max_intensities_transcription,mean_intensities_non, ...
    mean_intensities_transcription] = MeanMax2Channel(prefix1,prefix2)

minDist = 6;

% get experimental info for each dataset
    liveExperiment1 = LiveExperiment(prefix1);
    liveExperiment2 = LiveExperiment(prefix2);

    Spots1 = getSpots(liveExperiment1);
    Spots2 = getSpots(liveExperiment2);

    Particles1 = getParticles(liveExperiment1);
    Particles2 = getParticles(liveExperiment2);

 % find approved particles for both datasets

    %these are the not disapproved MS2 particles
    if iscell(Particles2)
        approved_parts2 = find([Particles2{2}.Approved] ~= -1);
    else 
        approved_parts2 = find([Particles2.Approved] ~= -1);
    end
    
    %these are the not disapproved parB particles
    if isfield(Particles1{2}, 'Approved')
        approved_parts1 = find([Particles1{2}.Approved] ~= -1);
    else
        approved_parts1 = 1:length(Particles1{2}); 
    end

% create array of all spot positions in transcription channel of the
% transcriptional-channel-segmented dataset
    ch2_spot_positions = [];
    %loop over MS2 particles
    for i = 1:length(approved_parts2)       
        par_num2 = approved_parts2(i);        
        %loop over frames of this MS2 particle
        if iscell(Particles2) %for some reason sometimes it is not
            for j = 1:length(Particles2{1,2}(par_num2).Frame)
                pos = [ Particles2{1,2}(par_num2).xPos(j),Particles2{1,2}(par_num2).yPos(j), Particles2{1,2}(par_num2).zPos(j)];
                ch2_spot_positions = [ch2_spot_positions; pos];
            end
        else
            for j = 1:length(Particles2(par_num2).Frame)
                pos = [ Particles2(par_num2).xPos(j),Particles2(par_num2).yPos(j), Particles2(par_num2).zPos(j)];
                ch2_spot_positions = [ch2_spot_positions; pos];
            end
        end
    end

    max_intensities_transcription = [];

    max_intensities_non = [];

    mean_intensities_transcription = [];

    mean_intensities_non = [];

    spot_particle_nums = [];

    non_spot_particle_nums = [];

% iterate through all approved particles in the transcriptional channel of
% the DNA-label-segmented dataset

     for k = 1:length(approved_parts1)

        par1_num = approved_parts1(k);

        frames = Particles1{1,2}(par1_num).Frame;

        indices = Particles1{1,2}(par1_num).Index;

        transcription_or_not = [];

        particle_intensities = [];

% iterate through all frames of each approved particle in the
% DNA-label-segmented dataset

        for h = 1:length(frames)

            ch1_filt_pos = [Particles1{1,2}(par1_num).xPos(h),Particles1{1,2}(par1_num).yPos(h), Particles1{1,2}(par1_num).zPos(h)];

            particle_intensities = [particle_intensities, Spots1{1,2}(frames(h)).Fits(indices(h)).FixedAreaIntensity3];

           
 % iterate through all spot positions in the
 % transcription-channel-segmented dataset. If a spot in the transcription channel of the
 % transcription-channel-segemnted dataset is within 2 pixels of a spot in the transcription
 % channel of the DNA-label-segmented data set, the variable
 % 'transcription_or_not' will get a 1 appended to it, if not, a 0 will get
 % appended. For each particle, if 'transcription_or_not' is all 0s, the
 % particle will be classified as non-transcriptional. If not, the particle
 % will be classified as transcriptional. 
 
            for p = 1:length(ch2_spot_positions)

               spot2_pos = ch2_spot_positions(p,:,:);

               dist  = sqrt(sum((ch1_filt_pos - spot2_pos).^2));

               if dist < minDist

                   transcription_or_not = [transcription_or_not,1];

               else 

                   transcription_or_not = [transcription_or_not,0];

               end

            end
        end    

       particle_intensities = particle_intensities(~isnan(particle_intensities)); 

 % calculate mean and 95th percentile for each particle and append it to
 % the appropriate output variable based on whether or not the particle is
 % transcriptional
       
       if sum(transcription_or_not) == 0

           non_spot_particle_nums = [non_spot_particle_nums, par1_num];

           max_intensities_non = [max_intensities_non, particle_intensities(particle_intensities >= prctile(particle_intensities,95))];

           mean_intensities_non = [mean_intensities_non, mean(particle_intensities)];

       else

           spot_particle_nums = [spot_particle_nums, par1_num];

           max_intensities_transcription = [max_intensities_transcription, particle_intensities(particle_intensities >= prctile(particle_intensities,95))];

           mean_intensities_transcription = [mean_intensities_transcription, mean(particle_intensities)]; 

       end
     end
end