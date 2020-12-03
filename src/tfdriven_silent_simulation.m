
close all force;

cmap = colormap(viridis(nPlots));

dmax = 5000;
nPlots = 10;
dl = 500;
dls = logspace(1, log10(dmax)); %aus
t = linspace(0, 10)'; %mins
kd = 500;
kds = logspace(2, 4, nPlots);
cs = logspace(-1, 2, nPlots);
pi1s = logspace(-2, 1, nPlots);

R = 500;
% c = .002;
c = 20;

t_cycle = 10; %min

nSteps = 10;
nSims = 1E4;
nOffStates = 5;
onstate = nOffStates+1;
silentstate = onstate+1;
occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );
pi0 = 1; %min
pi1 = .1; %min

dls = 700; %au
mfpts = nan(length(dls), length(kds), length(pi1s));

for k = 1:length(pi1s)
    for i = 1:length(dls)
        for j = 1:length(kds)
            pi0 = c.*occupancy(dls(i), kds(j)); %min-1
            pi1 = pi1s(k); %min-1

            fpt_on_observed = averagePaths(nSims, nSteps, pi0, pi1, onstate, silentstate, t_cycle);

            mfpts(i, j, k) = mean(fpt_on_observed);
        end
    end
end



figure;
for k = 1:length(pi1s)
%     plot(kds, mfpts(40, :, 1));
    plot(kds, mfpts(1, :, k), 'LineWidth', 2, 'Color', cmap(k, :));
    hold on
end
xlabel('K_D (au)')
ylabel('mean time to turn on (min)')
leg =legend(num2str(round2(pi1s')));
title(leg, '\pi_1');

figure;
plot(pi1s, mfpts(:, 5));
xlabel('super off transition rate (min-1)')
ylabel('mean time to turn on (min)')

% 
% figure;
% plot(dls, mfpts(:, 5));
% xlabel('[Dl] (au)')
% ylabel('mean time to turn on (min)')

function [states, times] = makePath(nSteps, pi0, pi1, onstate, silentstate)

% times = 0;
% states = 1;
state = 1;

states = nan(1, nSteps);
times = nan(1, nSteps);
r0 = exprnd(pi0^-1, [1, nSteps]); %min
r1 = exprnd(pi1^-1, [1, nSteps]); %min

for step = 1:nSteps
    
    [~, ind] = min([r0(step), r1(step)]);
    
    if ind == 1 && state < onstate
        state =  state + 1;
        tau = r0(step);
    elseif ind == 2 || state == onstate
        state = silentstate;
        tau = r1(step);
    end
    
    states(step) = state;
    
    if step ~= 1
         times(step) = times(step-1) + tau;
    else
        times(step) = 1;
    end
    
end

end

function fpt_on_observed = averagePaths(nSims, nSteps, pi0, pi1, onstate, silentstate, t_cycle)

    fpt_on = [];
%     fpt_on_observed = [];
%     duration = [];
    
    for k = 1:nSims
        
        [states, times] = makePath(nSteps, pi0, pi1, onstate, silentstate);
        
        % plot(time, states);
        % xlim([0, 10]);
        ton =  times(find(states==onstate, 1 ));
        fpt_on = [fpt_on, ton];
        
    end
    
    fpt_on_observed = fpt_on(fpt_on < t_cycle);
%     duration = t_cycle-fpt_on_observed;
    
        %     histogram(fpt_on_observed, 'Normalization', 'pdf');
    %     legend(['<\tau>=', num2str(round2(mean(fpt_on_observed))), ' min']);
    %     xlabel('\tau (min)')
    %     ylabel('probability density function')
    
    
end

