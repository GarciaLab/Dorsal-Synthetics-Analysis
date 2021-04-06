function [fraction,meanonset] = BasicModel_masterEq(c,kd,dorsalVals)

%Solve the master equation for state of promoter before transcription onset

%% set up problem
% Model parameters:
% c=1;    % transition rate 1/min
% kd = 1E3;
% dls = 1000;

%Simulation paramaters
numCells = 100;   % the total number of nuclei in the simulation
dt = 0.0005;        %Smaller than all time scales in the system
TotalTime = 8;   % end of the simulation
NumStates = 5;    %number of states

%Create the matrix to store the results
M(1:TotalTime/dt,1:NumStates+1) = 0; %initialize to zero everywhenre

%Initial conditions:
M(1,1)=numCells;    %everyone is at state 1 initially

%% Do the calculation
for d = 1:length(dorsalVals)
    dls = dorsalVals(d);
    k = (c*(dls./kd) ./ (1 + dls./kd));
    for t=2:TotalTime/dt % loop over time steps

        %Calculate the evolution of all boxes minus the ones at the edges
        for s=2:NumStates % loop over states
            stay = M(t-1,s);
            leave = k*dt*M(t-1,s);
            enter = k*dt*M(t-1,s-1);

            M(t,s) = stay + enter - leave;        
        end

        %Calculate the first box
        M(t,1) = M(t-1,1) - k*dt*M(t-1,1);

        %Calculate the last box
        M(t,NumStates+1) = M(t-1,NumStates+1) + k*dt*M(t-1,NumStates);
    end

    fraction(d) = M(end,end)/100;
    y6 = M(:,end); %number of nuclei in the last state as a function of time
    t = 0:dt:TotalTime-dt;
    meanonset(d) = sum(diff(y6).*t(1:end-1)')/sum(diff(y6)); %expected value
end

% %% Make a movie
% MVector=0:NumStates;     %This is the vector of bins for the histogram
% figure(1)
% for t=1:TotalTime/dt
%    bar(MVector,M(t,:))
%    ylim([0,100])
%    drawnow          %Force Matlab to draw the plot 
% end















