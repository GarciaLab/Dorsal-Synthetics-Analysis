
model = "basic";
% model = "entry";
% model = "exit";
% model = "entryexit";

rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period
dmax = 5000;
nPlots = 8;


t_cycle = 10; %min

nSteps = 7;
nSims = 1E3;
nOffStates = 5;
nEntryStates = 5;
firstoffstate = nEntryStates+1;
onstate = nEntryStates + nOffStates+1;
silentstate = onstate+1;
nStates = nEntryStates + nOffStates + 1 + 1;
occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );

if model == "entry"
    dls = logspace(log10(1), log10(dmax), 20);
    kds = logspace(2, 5, nPlots);
    cs = logspace(1, 3, nPlots);
    %     cs = 1;
    %     pi1s = logspace(-2, 1, nPlots);
    pi1s = 0;
    pi2s = logspace(-2, 1, nPlots);
    %     pi2s = 100;
elseif model == "entryexit"
    nSims = 1E3;
    %     dls = linspace(1, dmax, 20);
    dls = logspace(log10(1), log10(dmax), 40);
    kds = logspace(2, 6, nPlots);
    cs = logspace(0, 4, nPlots);
    %     cs = 1;
    pi1s = logspace(-2, 1, nPlots);
    pi2s = logspace(-2, 1, nPlots);
    %     pi2s = 100;
elseif model == "basic"
    nSims = 1E3;
    %     dls = linspace(1, dmax, 20);
    dls = logspace(log10(1), log10(dmax), 20);
    kds = logspace(2, 6, nPlots*3);
    cs = logspace(0, 4, nPlots*3);
    pi1s = 0;
    pi2s = 1E10;
elseif model == "exit"
    nSims = 1E4;
    dls = logspace(log10(1), log10(dmax), 100);
    kds = logspace(2, 6, nPlots);
    cs = logspace(0, 4, nPlots);
    pi1s = logspace(-2, 1, nPlots);
    pi2s = 1E10;
end

nParams = numel(dls)*numel(kds)*numel(pi1s)*numel(pi2s)*numel(cs);

[userview,~] = memory;
userview.MaxPossibleArrayBytes;

tau_exit = nan(nSteps, nSims, numel(pi1s), 'single');
for k = 1:length(pi1s)
    tau_exit(:, :, k) = exprnd(pi1s(k)^-1, [nSteps, nSims]);
end

tau_entry = nan(nEntryStates, nSims, numel(pi2s), 'single');
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


%%
mfpts = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));
factive = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));

%dls, kds, pi1s, cs, pi2s
for m = 1:length(cs)
    for i = 1:length(dls)
        for j = 1:length(kds)
                                
            pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1    
            tau_on = exprnd(pi0^-1, [nSteps, nSims]);

            for n = 1:length(pi2s)
                for k = 1:length(pi1s)

                    
                    [fpt_on_observed, factive_temp] = averagePaths_entryexit(...
                        nSims, nSteps, pi0, pi1s(k), pi2s(n),onstate, silentstate, t_cycle,...
                        firstoffstate,...
                        squeeze(tau_entry(:, :, n)), squeeze(tau_exit(:, :, k)), tau_on );
                    
                    
                    mfpts(i, j, k, m, n) = nanmean(fpt_on_observed);
                    
                    factive(i,j,k,m,n) = factive_temp;
                    
                end %for pi1
            end %for pi2
        end% for kd
    end %for dl
end %for c

dt = mfpts(:, nearestIndex(kds, 1E4), :, :, :) -...
    mfpts(:, nearestIndex(kds, 400), :, :, :); %kd(10)=10k, kd(4)=400

dt = repmat(dt, [1 length(kds) 1 1 1]);


[~, dropboxfolder] = getDorsalFolders;

save(dropboxfolder + "\" + "tf_paramsearch_"+model+"_.mat")

load(dropboxfolder + "\" + "tf_paramsearch_"+model+"_.mat")

figure;
try
    %     plotTFDrivenParams(factive, dt, mfpts, 'nPoints', 2E4)
    plotTFDrivenParams(factive, dt, mfpts, 'nPoints', 2E4, 'dim', 2, 'params', params)
    
catch
    plotTFDrivenParams(factive, dt, mfpts, 'dim', 2, 'params', params);
    %         plotTFDrivenParams(factive, dt, mfpts, 'dim', 2, 'params', params, 'nPoints', 100);
    
end

plotGoodCurves(factive, dt, mfpts, params)

%%
figure;
t = tiledlayout('flow');
for k = 1:2:length(dls)
    nexttile;
    plotTFDrivenParams(factive(k, :, :, :, :) , dt(k, :, :, :, :), mfpts(k, :, :, :, :), 'fig', gcf);
end
title(t, 'Effect of [Dorsal] on parameter space');


figure;
t = tiledlayout('flow');
for k = 1:1:length(cs)
    nexttile;
    plotTFDrivenParams(factive(:, :, :, k, :) , dt(:, :, :, k, :),...
        mfpts(:, :, :, k, :), 'fig', gcf);
end
title(t, 'Effect of c on parameter space');

figure;
t = tiledlayout('flow');
for k = 1:1:length(kds)
    nexttile;
    plotTFDrivenParams(factive(:, k, :, :, :) ,...
        dt(:, k, :, :, :), mfpts(:, k, :, :, :),'nPoints', 1E3, 'fig', gcf);
end
title(t, 'Effect of KD on parameter space');

figure;
t = tiledlayout('flow');
for k = 1:1:length(pi1s)
    nexttile;
    try
        plotTFDrivenParams(factive(:, :, k, :, :) ,...
            dt(:, :, k, :, :), mfpts(:, :, k, :, :), 'nPoints', 1E3, 'fig', gcf);
        %          title(['\pi_{exit} = ', num2str(round2(pi1s(k))), ' min^{-1}'])
        title(num2str(round2(pi1s(k))))
    end
end
title(t, 'Effect of pi_exit on parameter space');



figure;
t = tiledlayout('flow');
for k = 1:1:length(pi2s)
    nexttile;
    %     try
    plotTFDrivenParams(factive(:, :, :, :, k),...
        dt(:, :, :, :, k), mfpts(:, :, :, :, k), 'fig', gcf, 'shouldRound', true);
    %     end
    %         title(['\pi_{entry} = ', num2str(round2(pi2s(k))), ' min^{-1}'])
    title(num2str(round2(pi2s(k))))
end
title(t, 'Effect of pi_entry on parameter space');
%%
goodMatrixIndices = plotTFDrivenParams(factive, dt, mfpts);

%extract parameter set good at low Dorsal
figure;
tiledlayout('flow')

for m = 1:5
    temp1 = sub2ind(goodMatrixIndices, find(goodMatrixIndices(:, 1) == m));
    nexttile;
    for j = 1:size(temp1, 1)
        goodIndex = goodMatrixIndices(temp1(j), :);
        temp2 = num2cell(goodMatrixIndices(temp1(j), :));
        % factive(temp2{:})
        for k = 1:length(params.dls)
            
            factive_theory(k) = factive(k, goodIndex(2), goodIndex(3),...
                goodIndex(4), goodIndex(5));
            
        end
        plot(params.dls, factive_theory, 'LineWidth', 2)
        hold on
        % title({"K_D = "+kds(g(2)),...
        %     "\pi_{exit} = " + pi1s(g(3)),...
        %     "c = " + cs(g(4)),...
        %     "\pi_{entry} = " + pi2s(g(5))})
    end
end




end