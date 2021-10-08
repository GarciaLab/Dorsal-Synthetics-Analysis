function fraction_onset = ParallelModel(dorsalVals,theta,modelOpts)

% This is for a model where the promoter has n switches in parallel, all of
% which can transition with some rate to from OFF to ON. The promoter turns
% on only if all the switches are on.

% the function solves the probability that each switch turns on and then
% multiplies them to obtain the probability that the promoter is on.

% theta specifies model parameters:
% theta = [c,kd,nSwitches,DlIndependentK,tcycle]

% To compare with data, we generated per embryo fraction active and mean onset
% times using generatePerEmbryoDataForFits(numBins) with numBins = 18. The
% output is stored in Dropbox:
% '/Users/simon_alamos/Dropbox/DorsalSyntheticsDropbox/manuscript/window/basic/dataForFitting/archive';

% use modelOpts.nSims = nSims;

% this function uses constant Dl values.

% if isempty(modelOpts.TimeVariantDorsalValues)
%     ARPath = 'C:\Users\owner\Dropbox\DorsalSyntheticsDropbox\manuscript\window/basic/dataForFitting/archive';
%     SAPath = '/Users/simon_alamos/Dropbox/DorsalSyntheticsDropbox/manuscript/window/basic/dataForFitting/archive';
%     fullMatFileName = [SAPath '/DorsalFluoTraces.mat'];
%     load(fullMatFileName);
%     modelOpts.TimeVariantDorsalValues = [DorsalFluoTraces.meanDorsalFluo];
%     modelOpts.TimeVariantAbsoluteTimes = DorsalFluoTraces(1).absoluteTime; %in seconds
%     modelOpts.middleBinValues = [DorsalFluoTraces.binValue];
% end

%% set up problem
if isempty(modelOpts)
    modelOpts.exitOnlyDuringOffStates = true;
    modelOpts.nSims = 1E3;
    modelOpts.modelType = 'entryexit';
    modelOpts.piForm = "cOcc"; %pi depends on Dl occupancy. The other option is "cdl"
end

%mcmcpred feeds in x as a struct. otherwise it's an array
if isstruct(dorsalVals)
    dorsalVals = dorsalVals.ydata(:, 1);
end

%% Simulation paramaters
numCells = modelOpts.nSims;   % the total number of nuclei in the simulation
TotalTime =  theta(end);% minutes, end of the simulation
transcriptionStart = 0.01; %minutes, delayed start of the transcriptional window
dt = TotalTime/80; %this 80 seems sufficient for all purposes. sorry for hardcoding.

NSwitches = round(theta(end-2));   % number of parallel switches 
c = theta(1); % scales the Dorsal occupancy
kd = theta(2); % dorsal binding affinity
basalK = theta(end-1); %the basal switching rate in the absence of Dl

%Create the matrix to store the number of nuclei in the ON state 
M(1:TotalTime/dt,1) = 0; %initialize it to zero everywhere
transcriptionStartRow = ceil(transcriptionStart/dt); %this is the row of M when we want transcription to stop

%Initial conditions: everyone is at the first state at the first dt.
%note that the system enters state one with a delay specified by "transcriptionStart"
M(transcriptionStartRow,1)=numCells; % 

time_vec = 2:TotalTime/dt; % just a vector to loop over

%this is to find the dorsal time trace that corresponds to each dorsal bin
%[~, idx] = min(abs(repmat(modelOpts.TimeVariantAbsoluteTimes, length(time_vec), 1)'- time_vec*dt*60 )); %t*dt*60 is actual time in sec

time_vec_2 = 0:dt:TotalTime-dt; %time vector 

n_dls = length(dorsalVals);

fraction_onset = nan(length(dorsalVals), 2); %to store the output


%% Do the calculation now
% figure
% hold on
% Palette = cbrewer('seq', 'Reds', length(1:n_dls));
% actualAbsTime = linspace(0,TotalTime,size(M,1)); %just for plotting

for d = 1:n_dls
    DlConc = dorsalVals(d);
    SpedUpK = basalK*(c*(DlConc./kd) ./ (1 + DlConc./kd)); % rate of the switch controlled by Dorsal
    
    SwitchesStates = zeros(NSwitches,length(time_vec_2));
    
    % calculate the Dorsal-regulated switch first
    SwitchesStates(1,:) = 1-exp(-SpedUpK.*time_vec_2);
    
    % calculate the rest of the unregulated switches
    for n = 2:NSwitches
        SwitchesStates(n,:) = 1-exp(-basalK.*time_vec_2);
    end
    
    % the fraction of promoters in the on state is the multiplication of
    % the probabilities of all switches combined.
    
    All3On = prod(SwitchesStates,1);
    
    % calculate the fraction active = the fraction of nuclei in the on
    % state by the end of the time window.
    fraction_onset(d,1) = All3On(end);  
    fraction_onset(d,2) = sum(diff(All3On).*time_vec_2(1:end-1)/sum(diff(All3On))); %expected value
    
end
% 
% figure
% subplot(1,2,1)
% plot(dorsalVals,fraction_onset(:,1),'r')
% xlabel('[Dl]')
% ylabel('fraction active')
% ylim([0 1.1])
% 
% subplot(1,2,2)
% plot(dorsalVals,fraction_onset(:,2),'b')
% xlabel('[Dl]')
% ylabel('mean turn on time')
% ylim([0 8.1])





% 
% 
% 
% 
% for d = 1:n_dls %loop over dorsal concentration bins
%     
%     [~,nearestBin] = min(abs(modelOpts.middleBinValues - dorsalVals(d)));
%     dorsalTraceFluo = modelOpts.TimeVariantDorsalValues(:,nearestBin);
%     
%     for t = time_vec(transcriptionStartRow:end) % loop over time steps        
%         %assert(abs(sum(M(t-1,:))-M(1,1))<errorTolerance,'the total probability across states should always add up to the initial one')        
%         dls = dorsalTraceFluo(idx(t-1)); % each dorsal bin has a corresponding concentration time trace
%         %dls = dorsalVals(d) + diff(dorsalVals(1:2)); % this is in case we want constant Dorsal       
%         if modelOpts.piForm == "cdl"
%             kdt_off = c*dls; % transition rate between off states
%         elseif modelOpts.piForm == "cOcc"
%             kdt_off = (c*(dls./kd) ./ (1 + dls./kd))*dt; % transition rate between off states
%         end
%         
%         %Calculate the first state
%         if NInactiveStates ~= 0 % if the first state is an inactive one          
%             M(t,1) = (1-kdt_inac)*M(t-1,1); % it transitions with a rate of kdt_inac
%             
%             %Calculate the evolution of all states minus the ones at the edges
%             % loop over inactive states
%             for s = 2:NInactiveStates
%                 M(t,s) = (1-kdt_inac)*M(t-1,s) + kdt_inac*M(t-1,s-1); % (stay-leave) + enter
%             end
%             
%             % do the off state that comes immediatly after the last inactive state
%             M(t,NInactiveStates+1) = (1-kdt_off)*M(t-1,NInactiveStates+1) + kdt_inac*M(t-1,NInactiveStates);       
%         
%         else % if the first state is an off one
%             M(t,1) = (1-kdt_off)*M(t-1,1); % it transitions with a rate of kdt_off
%         end
%         
%         % loop over the rest of the off states
%         for s = NInactiveStatesPlus2:NLinStates
%             M(t,s) = (1-kdt_off)*M(t-1,s) + kdt_off*M(t-1,s-1);
%         end
%         
%         %Calculate the last state
%         M(t,end) = M(t-1,end) + kdt_off*M(t-1,end-1);
%         
%     end
%     
%     fraction_onset(d,1) = M(end,end)/numCells;
%     yEnd = M(:,end); %number of nuclei in the last state as a function of time
%     fraction_onset(d,2) = sum(diff(yEnd).*time_vec_2(1:end-1)')/sum(diff(yEnd)); %expected value
%     
% %     plot(actualAbsTime,yEnd,'Color',Palette(d,:),'LineWidth',1)
% %     plot([fraction_onset(d,2) fraction_onset(d,2)],[0 0.5],'-','Color',Palette(d,:))
% %     %title(['Dorsal bin ' num2str(d)])
% end
% % hold off
% 
% 
% 
% 
% % %% Make a movie
% % MVector=0:NumStates;     %This is the vector of bins for the histogram
% % figure(1)
% % for t=1:TotalTime/dt
% %    bar(MVector,M(t,:))
% %    ylim([0,100])
% %    drawnow          %Force Matlab to draw the plot
% % end