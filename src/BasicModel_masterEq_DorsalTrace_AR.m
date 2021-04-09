function fraction_onset = BasicModel_masterEq_DorsalTrace_AR(dorsalVals,theta, modelOpts)

% Solves the master equation for state of promoter before transcription
% onset. It can handle active and inactive states

% To compare with data, we generated per embryo fraction active and onset
% times using generatePerEmbryoDataForFits(numBins) with numBins = 18. The
% output is stored in Dropbox:
% '/Users/simon_alamos/Dropbox/DorsalSyntheticsDropbox/manuscript/window/basic/dataForFitting/archive';

% An analytical solution to a subset of this model can be generated using the function 
% fiveOffSteps

% Theta contains parameters in this order: c, kd, Ninactive, Noff, piEntry, piExit, tCycle
% use modelOpts.nSims = nSims;

% Compared to BasicModel_masterEq, in this function dorsalVals is still an array of one value per bin,
% but we'll use these values to retrieve their corresponding Dorsal time traces.

% Now, at each dt, we'll use a time-variant Dorsal concentration
% these data was generated in advance using:
% generateDorsalTimeTraces(linspace(0,3800,18))
% it is stored in a .mat file that we retrieve here.

if isempty(modelOpts.TimeVariantDorsalValues)
    ARPath = 'C:\Users\owner\Dropbox\DorsalSyntheticsDropbox\manuscript\window/basic/dataForFitting/archive';
    SAPath = '/Users/simon_alamos/Dropbox/DorsalSyntheticsDropbox/manuscript/window/basic/dataForFitting/archive';
    fullMatFileName = [SAPath '/DorsalFluoTraces.mat'];
    load(fullMatFileName);
    modelOpts.TimeVariantDorsalValues = [DorsalFluoTraces.meanDorsalFluo];
    modelOpts.TimeVariantAbsoluteTimes = DorsalFluoTraces(1).absoluteTime; %in seconds
    modelOpts.middleBinValues = [DorsalFluoTraces.binValue];
end

%% set up problem
if isempty(modelOpts)
    modelOpts.exitOnlyDuringOffStates = true;
    modelOpts.nSims = 1E3;
    modelOpts.modelType = 'entryexit';
end

%mcmcpred feeds in x as a struct. otherwise it's an array
if isstruct(dorsalVals)
    dorsalVals = dorsalVals.ydata(:, 1);
end

%% Simulation paramaters
numCells = modelOpts.nSims;   % the total number of nuclei in the simulation
TotalTime =  theta(7);   % end of the simulation
dt = TotalTime/80; %this 80 seems sufficient for all purposes. sorry for hardcoding.
NOffStates = round(theta(4));   %number of off states
NInactiveStates = round(theta(3)); %number of inactive states
NLinStates = NOffStates+NInactiveStates;

% min(abs(modelOpts.TimeVariantAbsoluteTimes- (t*dt*60) ))


c = theta(1);
kd = theta(2);
nEntryStates = round(theta(3));
pi_entry = theta(5);
pi_exit = theta(6);

%Create the matrix to store the results
M(1:TotalTime/dt,1:(NInactiveStates+NOffStates+1)) = 0; %initialize to zero everywhenre

%Initial conditions:
M(1,1)=numCells;    %everyone is at the first state initially

time_vec = 2:TotalTime/dt; 
[~, idx] = min(abs(repmat(modelOpts.TimeVariantAbsoluteTimes, length(time_vec), 1)'- time_vec*dt*60 )); %t*dt*60 is actualtimeinsecs


time_vec_2 = 0:dt:TotalTime-dt;

n_dls = length(dorsalVals)-1;

fraction_onset = nan(length(dorsalVals)-1, 2);


%% Do the calculation now
for d = 1:n_dls %loop over dorsal bins
    
    dorsalFluo = dorsalVals(d);
    [~,nearestBin] = min(abs(modelOpts.middleBinValues - dorsalFluo));
    dorsalTraceFluo = modelOpts.TimeVariantDorsalValues(:,nearestBin);
    
    for t = time_vec % loop over time steps
        
        dls = dorsalTraceFluo(idx(t-1)); % each dorsal bin has a corresponding concentration time trace
        
        %dls = dorsalTraceFluo(400); % this is in case we want constant Dorsal
        
        k = (c*(dls./kd) ./ (1 + dls./kd));
            kdt_off = k*dt;
        kdt_inac = pi_entry*dt;

        %Calculate the evolution of all boxes minus the ones at the edges
        
        % loop over inactive states         
        for s = 2:NInactiveStates           
            M(t,s) = (1-kdt_inac)*M(t-1,s) + kdt_inac*M(t-1,s-1); %stay + enter - leave
        end
        
        % loop over off states
        for s = (NInactiveStates+2):(NInactiveStates+NOffStates)            
            M(t,s) = (1-kdt_off)*M(t-1,s) + kdt_off*M(t-1,s-1); %stay + enter - leave
        end

        
        %Calculate the first box
        if NInactiveStates
            M(t,1) = (1-kdt_inac)*M(t-1,1);
        else
            M(t,1) = (1-kdt_off)*M(t-1,1);
        end
                
        %Calculate the last box
        M(t,end) = M(t-1,end) + kdt_off*M(t-1,end-1);
        
    end
    
    fraction_onset(d,1) = M(end,end)/numCells;
    y6 = M(:,end); %number of nuclei in the last state as a function of time
    fraction_onset(d,2) = sum(diff(y6).*time_vec_2(1:end-1)')/sum(diff(y6)); %expected value
end


%plot(linspace(1,TotalTime,size(M,1)),M,'LineWidth',2)

% %% Make a movie
% MVector=0:NumStates;     %This is the vector of bins for the histogram
% figure(1)
% for t=1:TotalTime/dt
%    bar(MVector,M(t,:))
%    ylim([0,100])
%    drawnow          %Force Matlab to draw the plot
% end