function [fpt_on_observed,factive] = averagePaths_entryexit(nSims,...
    nSteps, pi0, pi1,pi2, onstate, silentstate,...
    t_cycle, firstoffstate, tau_entry, tau_exit, tau_on)
%subfunction for tfdrivenentryexit


fpt_on = nan(1, nSims);

if isempty(tau_entry)
    tau_entry = exprnd(pi2^-1, [5, nSims]);
end

if isempty(tau_exit)
    tau_exit = exprnd(pi1^-1, [nSteps, nSims]); 
end

if isempty(tau_on)
    tau_on = exprnd(pi0^-1, [nSteps, nSims]); 
end


[~, whichTransition] = min(cat(3,tau_on,tau_exit), [], 3);
whichTransition = squeeze(whichTransition);

initStates = [repmat([1, 2, 3, 4, 5, 6], 1, 1), zeros(1,6)];
initTimes = [zeros(nSims, 1), cumsum(tau_entry)', zeros(nSims,6)];

for k = 1:nSims
%parfor k = 1:nSims
    
    
    if initTimes(k, 6) < t_cycle
        
        [states, times] = makePath(nSteps,onstate, silentstate, firstoffstate,...
            tau_on(:, k), tau_exit(:, k), whichTransition(:, k), initStates, initTimes(k,:));
        
        % plot(time, states);
        % xlim([0, 10]);
        %     ton =  times(find(states==onstate, 1 ));
        ton = times(states==onstate);
        %     fpt_on = [fpt_on, ton];
        if ~isempty(ton)
            fpt_on(k) = ton;
        end
        
    end
    
end



fpt_on_observed = fpt_on(fpt_on < t_cycle);


factive = length(fpt_on_observed) / nSims;

% if ~isempty(fpt_on_observed)
%     
%     1
% end

%     duration = t_cycle-fpt_on_observed;

%     histogram(fpt_on_observed, 'Normalization', 'pdf');
%     legend(['<\tau>=', num2str(round2(mean(fpt_on_observed))), ' min']);
%     xlabel('\tau (min)')
%     ylabel('probability density function')


end

function [states, times] = makePath(nSteps, onstate, silentstate,firstoffstate,...
    tau_on, tau_exit, whichTransition, states, times)


n = 6;
for step = 1:nSteps
    n = n + 1;
    if whichTransition(step) == 1 && states(n-1) < onstate
        states(n) = states(n-1) + 1;
        times(n) = times(n-1) + tau_on(step);
    elseif whichTransition(step) == 2 || states(n-1) == onstate
        states(n) = silentstate;
        times(n) = times(n-1) + tau_exit(step);
        return;
    end
end

end

