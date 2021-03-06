function fraction_onset = BasicModel_masterEq_DorsalTrace(dorsalVals,c,kd,n_entry,n_off,piEntry,piExit,tCycle,options)


% Solve the master equation for state of promoter before transcription onset

% compared to BasicModel_masterEq, in this function dorsalVals is still an array of one value per bin, 
% but we'll use these values to retrieve their corresponding Dorsal time traces.

% now, at each dt, we use a time-variant Dorsal concentration

% these data was generated in advance and stored in a .mat file that we
% retrieve here.
Path = '/Users/simon_alamos/Dropbox/DorsalSyntheticsDropbox/manuscript/window/basic/dataForFitting/archive';
fullMatFileName = [Path '/DorsalFluoTraces.mat'];
load(fullMatFileName);
TimeVariantDorsalValues = [DorsalFluoTraces.meanDorsalFluo];
TimeVariantAbsoluteTimes = DorsalFluoTraces(1).absoluteTime; %in seconds

%% set up problem
% Model parameters:
% c=1;    % transition rate 1/min
% kd = 1E3;
% dls = 1000;

%Simulation paramaters
numCells = 100;   % the total number of nuclei in the simulation
dt = 0.1;        % Smaller than all time scales in the system
TotalTime = tCycle;   % end of the simulation in minutes
NOffStates = n_off;    %number of states

%Create the matrix to store the results
M(1:TotalTime/dt,1:NOffStates+1) = 0; %initialize to zero everywhenre

%Initial conditions:
M(1,1)=numCells;    %everyone is at state 1 initially

%% Do the calculation
for d = 1:(length(dorsalVals)-1)
    
    dorsalTraceFluo = TimeVariantDorsalValues(:,d);

    for t=2:TotalTime/dt % loop over time steps
        
        actualTimeInSec = t*dt*60; 
        [delta, idx] = min(abs(TimeVariantAbsoluteTimes-actualTimeInSec));
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

    fraction_onset(d,1) = M(end,end)/100;
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















