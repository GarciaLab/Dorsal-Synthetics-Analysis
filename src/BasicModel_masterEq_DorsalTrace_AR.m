function fraction_onset = BasicModel_masterEq_DorsalTrace_AR(dorsalVals,theta, modelOpts)

% Solve the master equation for state of promoter before transcription onset

% compared to BasicModel_masterEq, in this function dorsalVals is still an array of one value per bin,
% but we'll use these values to retrieve their corresponding Dorsal time traces.

% now, at each dt, we use a time-variant Dorsal concentration

% these data was generated in advance and stored in a .mat file that we
% retrieve here.
if isempty(modelOpts.TimeVariantDorsalValues)
    ARPath = 'C:\Users\owner\Dropbox\DorsalSyntheticsDropbox\manuscript\window/basic/dataForFitting/archive';
    %SAPath = '/Users/simon_alamos/Dropbox/DorsalSyntheticsDropbox/manuscript/window/basic/dataForFitting/archive';
    fullMatFileName = [ARPath '/DorsalFluoTraces.mat'];
    load(fullMatFileName);
    TimeVariantDorsalValues = [DorsalFluoTraces.meanDorsalFluo];
    TimeVariantAbsoluteTimes = DorsalFluoTraces(1).absoluteTime; %in seconds
    middleBinValues = [DorsalFluoTraces.binValue];
end

%% set up problem
% Model parameters:
% c=1;    % transition rate 1/min
% kd = 1E3;
% dls = 1000;

if isempty(modelOpts)
    modelOpts.exitOnlyDuringOffStates = true;
    modelOpts.nSims = 1E3;
    modelOpts.modelType = 'entryexit';
end

%mcmcpred feeds in x as a struct. otherwise it's an array
if isstruct(dorsalVals)
    dorsalVals = dorsalVals.ydata(:, 1);
end

%Simulation paramaters
numCells = modelOpts.nSims;   % the total number of nuclei in the simulation
TotalTime =  theta(7);   % end of the simulation
% dt = 0.1;        %Smaller than all time scales in the system
dt = TotalTime/80; %this 80 seems sufficient for all purposes. sorry for hardcoding.
NOffStates = round(theta(4));   %number of states


c = theta(1);
kd = theta(2);
nEntryStates = round(theta(3));
pi_entry = theta(5);
pi_exit = theta(6);

%Create the matrix to store the results
M(1:TotalTime/dt,1:NOffStates+1) = 0; %initialize to zero everywhenre

%Initial conditions:
M(1,1)=numCells;    %everyone is at state 1 initially

fraction_onset = nan(length(dorsalVals), 2);

%% Do the calculation
for d = 1:length(dorsalVals)
    dorsalFluo = dorsalVals(d);
    [~,nearestBin] = min(abs(middleBinValues - dorsalFluo));
    dorsalTraceFluo = modelOpts.TimeVariantDorsalValues(:,nearestBin);
    
    for t=2:TotalTime/dt % loop over time steps
        
        [~, idx] = min(abs(modelOpts.TimeVariantAbsoluteTimes- (t*dt*60) )); %t*dt*60 is actualtimeinsecs
        dls = dorsalTraceFluo(idx);
        k = (c*(dls./kd) ./ (1 + dls./kd));
        
        %Calculate the evolution of all boxes minus the ones at the edges
        for s=2:NOffStates % loop over states
            stay = M(t-1,s);
            leave = k*dt*M(t-1,s);
            enter = k*dt*M(t-1,s-1);
            
            M(t,s) = stay + enter - leave;
        end
        
        %Calculate the first box
        M(t,1) = M(t-1,1) - k*dt*M(t-1,1);
        
        %Calculate the last box
        M(t,NOffStates+1) = M(t-1,NOffStates+1) + k*dt*M(t-1,NOffStates);
    end
    
    fraction_onset(d,1) = M(end,end)/numCells;
    y6 = M(:,end); %number of nuclei in the last state as a function of time
    t = 0:dt:TotalTime-dt;
    fraction_onset(d,2) = sum(diff(y6).*t(1:end-1)')/sum(diff(y6)); %expected value
end

% %% Make a movie
% MVector=0:NumStates;     %This is the vector of bins for the histogram
% figure(1)
% for t=1:TotalTime/dt
%    bar(MVector,M(t,:))
%    ylim([0,100])
%    drawnow          %Force Matlab to draw the plot
% end