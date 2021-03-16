%leaving this here for to remember for testing
%setpref('profiler','showJitLines',1);
%profile -memory on;

tic

% model = "basic";
% model = "entry";
% model = "exit";
model = "entryexit";

rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period
dmax = 4000;


t_cycle = 8; %min

nOffStates = 5;
nEntryStates = 5;
nOffEntryStates = nOffStates + nEntryStates;
firstoffstate = nEntryStates+1;
onstate = nEntryStates + nOffStates+1;
silentstate = onstate+1;
nStates = nEntryStates + nOffStates + 1 + 1;
occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );


exitOnlyDuringOffStates = true;
nSims = 1E5;
nPlots = 20;

dls = logspace(0, log10(dmax), 80);
kds = logspace(2, 7, nPlots);
cs = logspace(-5, 2, nPlots);
pi1s = logspace(-3, 2, nPlots);
pi2s = logspace(-1, 2, nPlots);

switch model
    case "entry"
        pi1s = 0;
    case "basic"
        pi1s = 0;
        pi2s = 1E10;
    case "exit"
        pi2s = 1E10;
end

nParams = numel(dls)*numel(kds)*numel(pi1s)*numel(pi2s)*numel(cs);

tau_exit = nan(nStates-1, nSims, numel(pi1s), 'double');
for k = 1:length(pi1s)
    tau_exit(:, :, k) = exprnd(pi1s(k)^-1, [nStates-1, nSims]);
end

tau_entry = nan(nEntryStates, nSims, numel(pi2s), 'double');
for k = 1:length(pi2s)
    tau_entry(:, :, k) = exprnd(pi2s(k)^-1, [nEntryStates, nSims]);
end

clear params;
params.dls = dls;
params.kds = kds;
params.cs = cs;
params.pi1s = pi1s;
params.pi2s = pi2s;
params.model = model;
params.nEntryStates = nEntryStates;
params.nOffStates = nOffStates;
params.nStates = nStates;
params.exitOnlyDuringOffStates = exitOnlyDuringOffStates;


%%
mfpts = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));
factive = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));
fpts_std = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));

N_cs = length(params.cs);
N_dls = length(params.dls);
N_kds = length(params.kds);
N_pi1s = length(params.pi1s);
N_pi2s = length(params.pi2s);

%dls, kds, pi1s, cs, pi2s
for m = 1:N_cs
    
    display("tfdrivenentryexit progress: "+ num2str(( (m-1) / N_cs)*100)+"%" )
    
    for i = 1:N_dls
        for j = 1:N_kds
            
            pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1
            tau_on = exprnd(pi0^-1, [nOffStates+1, nSims]);
            
            for n = 1:N_pi2s
                
                tau_entry_off = [squeeze(tau_entry(:, :, n)); tau_on];
                
                for k = 1:N_pi1s
                    
                    %let's determine if the transition is to the silent
                    %state or other.
                    [~, whichTransition] = min(cat(3,tau_entry_off,squeeze(tau_exit(:, :, k))), [], 3);
                    
                    if exitOnlyDuringOffStates
                        whichTransition(1:nEntryStates, :) = 1;
                    end
                    %the simulations that reached "on" are the ones that
                    %never reached silent.
                    reachedOn = sum(whichTransition(1:nOffEntryStates, :), 1) ==...
                        nOffEntryStates;

                    %let's get the total duration of the successful trajectories up to
                    %the on state
                    onsets_sim = sum(tau_entry_off(1:nOffEntryStates, reachedOn), 1);
                    onsets_sim_truncated = onsets_sim(onsets_sim < t_cycle);
                    
                    reachedOn_truncated = reachedOn(onsets_sim < t_cycle);
                    
                    factive(i,j,k,m,n) = sum(reachedOn_truncated)/nSims;
                    mfpts(i, j, k, m, n) = mean(onsets_sim_truncated);
                    fpts_std(i, j, k, m, n) = std(onsets_sim_truncated);
                    
                    %                     [fpt_on_observed, factive_temp] = averagePaths_entryexit(...
                    %                         nSims, nSteps, pi0, pi1s(k), pi2s(n),onstate, silentstate, t_cycle,...
                    %                         firstoffstate,...
                    %                         squeeze(tau_entry(:, :, n)), squeeze(tau_exit(:, :, k)), tau_on, params );
                    %
                end %for pi2
            end %for pi1
        end% for kd
    end %for dl
end %for c

dt = mfpts(:, nearestIndex(kds, 1E4), :, :, :) -...
    mfpts(:, nearestIndex(kds, 400), :, :, :);

dt = repmat(dt, [1 length(kds) 1 1 1]);


[~, dropboxfolder] = getDorsalFolders;

saveStr = model;
if exitOnlyDuringOffStates
    saveStr = saveStr + "_exitOnlyOnOffStates";
end
save(dropboxfolder + "\" + "tf_paramsearch_"+saveStr+"_.mat")

% %Also save a timestamped copy.
% dttm = strrep(strrep(string(datetime), " ", "_"), ":", "_");
% save(dropboxfolder + "\" + "tf_paramsearch_"+model+"_"+dttm+"_.mat")


toc
% load(dropboxfolder + "\" + "tf_paramsearch_"+saveStr+"_.mat")

figure;
try
plotTFDrivenParams(factive, dt, mfpts, 'nPoints', nPoints, 'dim', 2, 'params', params)
catch
    
plotTFDrivenParams(factive, dt, mfpts, 'dim', 2, 'params', params)
end
goodMatrixIndices = plotTFDrivenParams(factive, dt, mfpts, 'dim', 2, 'params', params);
plotGoodCurves(factive, dt, mfpts, params, goodMatrixIndices)

%%
% try
%     figure;
%     t = tiledlayout('flow');
%     for k = 1:2:length(dls)
%         nexttile;
%         plotTFDrivenParams(factive(k, :, :, :, :) , dt(k, :, :, :, :),...
%             mfpts(k, :, :, :, :), 'params', params, 'fig', gcf);
%     end
%     title(t, 'Effect of [Dorsal] on parameter space');
%
%
%     figure;
%     t = tiledlayout('flow');
%     for k = 1:1:length(cs)
%         nexttile;
%         plotTFDrivenParams(factive(:, :, :, k, :) , dt(:, :, :, k, :),...
%             mfpts(:, :, :, k, :),'params', params, 'fig', gcf);
%     end
%     title(t, 'Effect of c on parameter space');
%
%     figure;
%     t = tiledlayout('flow');
%     for k = 1:1:length(kds)
%         nexttile;
%         plotTFDrivenParams(factive(:, k, :, :, :) ,...
%             dt(:, k, :, :, :), mfpts(:, k, :, :, :),'params', params,'nPoints', 1E3, 'fig', gcf);
%     end
%     title(t, 'Effect of KD on parameter space');
%
%     figure;
%     t = tiledlayout('flow');
%     for k = 1:1:length(pi1s)
%         nexttile;
%         try
%             plotTFDrivenParams(factive(:, :, k, :, :) ,...
%                 dt(:, :, k, :, :), mfpts(:, :, k, :, :),'params', params, 'nPoints', 1E3, 'fig', gcf);
%             %          title(['\pi_{exit} = ', num2str(round2(pi1s(k))), ' min^{-1}'])
%             title(num2str(round2(pi1s(k))))
%         end
%     end
%     title(t, 'Effect of pi_exit on parameter space');
%
%
%
%     figure;
%     t = tiledlayout('flow');
%     for k = 1:1:length(pi2s)
%         nexttile;
%         %     try
%         plotTFDrivenParams(factive(:, :, :, :, k),...
%             dt(:, :, :, :, k), mfpts(:, :, :, :, k), 'params', params, 'fig', gcf, 'shouldRound', true);
%         %     end
%         %         title(['\pi_{entry} = ', num2str(round2(pi2s(k))), ' min^{-1}'])
%         title(num2str(round2(pi2s(k))))
%     end
%     title(t, 'Effect of pi_entry on parameter space');
% end