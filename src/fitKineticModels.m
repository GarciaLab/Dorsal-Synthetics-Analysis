function [results, chain, s2chain, data, modelOpts] = fitKineticModels(varargin)

% To do:
% 2. Get variable state number walkers to walk properly
% 3. Fit multiple KDs simultaneously using batched fits.


wb = true;
nSteps = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats
nSims = 1E3; %number of simulations for the kinetic barrier model (not the number of mcmc walker steps).
exitOnlyDuringOffStates = true; %determines connectivity of the markov graph
modelType = "entryexit"; %choices- entryexit, entry, exit, basic
fun= "table"; %also 'sim'
t_cycle = 8;
variableStateNumber = false;
fixKD = false;
batchedAffinities = false;


%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

[~, dropboxfolder] = getDorsalFolders;

%1Dg data (highest affinity)
datPath = dropboxfolder + "\manuscript\window\basic\dataForFitting\archive\";
load(datPath + "DorsalFluoValues.mat", "DorsalFluoValues");
load(datPath + "FractionsPerEmbryo.mat", "FractionsPerEmbryo");
load(datPath + "TimeOnsPerEmbryo.mat", "TimeOnsPerEmbryo");


if batchedAffinities
    
    %this setting being true wouldn't make sense, so let's fix it.
    fixKD = false;
    
    DorsalFluoValues0 = DorsalFluoValues;
    FractionsPerEmbryo0 = FractionsPerEmbryo;
    TimeOnsPerEmbryo0 = TimeOnsPerEmbryo;
    
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
    scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
    data_batch = {};
    
    
    for k = 1:numel(enhancers)
        %         datPath = dropboxfolder + "\manuscript\window\basic\dataForFitting\archive\";
        %         load(datPath + "DorsalFluoValues.mat", "DorsalFluoValues");
        %         load(datPath + "FractionsPerEmbryo.mat", "FractionsPerEmbryo");
        %         load(datPath + "TimeOnsPerEmbryo.mat", "TimeOnsPerEmbryo");
        
        %needs to be nobs x ny
        
        %         X = repmat(DorsalFluoValues, 1, max(size(FractionsPerEmbryo)));
        %         data{k}.ydata = [X; FractionsPerEmbryo(:)'; TimeOnsPerEmbryo(:)']';
        %       X = repmat(DorsalFluoValues, 1, max(size(FractionsPerEmbryo)));
        
        %     data{k}.ydata = [X; FractionsPerEmbryo(:)'; TimeOnsPerEmbryo(:)']';
    end
    %1Dg data (highest affinity)
    %     datPath = dropboxfolder + "\manuscript\window\basic\dataForFitting\archive\";
    %     load(datPath + "DorsalFluoValues.mat", "DorsalFluoValues");
    datf1 = load(datPath + "FractionsPerEmbryo.mat", "FractionsPerEmbryo");
    datt1 = load(datPath + "TimeOnsPerEmbryo.mat", "TimeOnsPerEmbryo");
    X = repmat(DorsalFluoValues, 1, max(size(FractionsPerEmbryo)));
    data{1}.ydata = [X; FractionsPerEmbryo(:)'; TimeOnsPerEmbryo(:)']';
    
    clear FractionsPerEmbryo
    clear TimeOnsPerEmbryo
    FractionsPerEmbryo{1} = datf1.FractionsPerEmbryo;
    TimeOnsPerEmbryo{1} = datt1.TimeOnsPerEmbryo;
    
    datf2 = load(datPath + "FractionsPerEmbryo_1DgVW.mat", "FractionsPerEmbryo");
    datt2 = load(datPath + "TimeOnsPerEmbryo_1DgVW.mat", "TimeOnsPerEmbryo");
    FractionsPerEmbryo{2} = datf2.FractionsPerEmbryo;
    TimeOnsPerEmbryo{2} = datt2.TimeOnsPerEmbryo;
    X2 = repmat(DorsalFluoValues, 1, max(size(datf2.FractionsPerEmbryo)));
    data{2}.ydata = [X2; datf2.FractionsPerEmbryo(:)'; datt2.TimeOnsPerEmbryo(:)']';
    
    
else
    
    %needs to be nobs x ny
    X = repmat(DorsalFluoValues, 1, max(size(FractionsPerEmbryo)));
    data.ydata = [X; FractionsPerEmbryo(:)'; TimeOnsPerEmbryo(:)']';
end



%%
saveStr = modelType;
if exitOnlyDuringOffStates
    saveStr = saveStr + "_exitOnlyOnOffStates";
end
if variableStateNumber
    saveStr = saveStr + "variableStateNumber";
end
% sims = load(dropboxfolder +  "\simulations\archive\" + "tf_paramsearch_"+saveStr+"_.mat", 'params', 'factive', 'mfpts');
if fun == "table"
    sims = load(dropboxfolder +  "\simulations\" + "tf_paramsearch_"+saveStr+"_.mat", 'params', 'factive', 'mfpts');
end
%%
rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period
% options.drscale = 2; % a high value (5) is important for multimodal parameter spaces.
options.drscale = 5;
options.waitbar = wb; %the waitbar is rate limiting sometimes
options.nsimu = nSteps; %should be between 1E3 and 1E6
options.updatesigma = 1; %honestly don't know what this does
% options.method = 'mh';
%
names = ["c", "kd" , "nentrystates", "moffstates", "pentry", "pexit", "tcycle"];
if modelType == "entryexit"
    if ~batchedAffinities
        p0 = [10, 1E3, 5, 5, 1, 1, 8];
        lb = [1E-1, 1E2, 0, 1, 1E-1, 1E-2, 5];
        ub = [1E2, 1E6, 12, 12, 1E3, 1E1, 10];
    else
        p0 = [1, 1E3, 1, 1, 1, .01, 8];
        lb = [1E-1, 1E1, 0, 1, 1E0, 1E-3, 7];
        ub = [1E2, 1E4, 12, 12, 1E2, 1E-2, 9];
    end
elseif modelType == "entry"
    p0 = [10, 1E3, 5, 5, 1, 0, 8];
    lb = [1E-2, 1E0, 0, 1, 1E-1, 0, 5];%pentry lower than .1 causes crash
    ub = [1E2, 1E6, 12, 12, 1E1, 0, 10];
elseif modelType == "basic"
    p0 = [100, 1E3, 0, 1, 1E10, 0, 8];
    lb = [1E-2, 1E0, 0, 1, 1E10, 0, 5];
    ub = [1E3, 1E6, 0, 12, 1E10, 0, 10];
end

if variableStateNumber
    lb(1) = 10;
end

params = cell(1, length(p0));
pri_mu = NaN; %default prior gaussian mean
pri_sig = Inf; %default prior gaussian variance

for k = 1:length(names)
    
    targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.
    localflag = 0; %is this local to this dataset or shared amongst batches?
    
    if  contains(names(k), "states") && ~variableStateNumber
        targetflag = 0;
    end
    
    if ~contains(modelType, "exit")
        if names(k) == "pexit"
            targetflag = 0;
            p0(k) = 0;
        end
    end
    
    if ~contains(modelType, "entry")
        if names(k) == "pentry" || names(k) == "nentrystates"
            targetflag = 0;
            p0(k) = 0;
        end
    end
    
    if fixKD && names(k) == "kd"
        targetflag = 0;
    end
    
    %we allow each enhancer to have a different kd but share every other
    %param
    if batchedAffinities && names(k) == "kd"
        localflag = 1;
    end
    
    params{1, k}= {names(k), p0(k), lb(k), ub(k), pri_mu, pri_sig, targetflag, localflag};
end


modelOpts = struct;
modelOpts.modelType = modelType;
modelOpts.t_cycle = t_cycle;

if fun == "table"
    modelOpts.sims = sims;
else
    modelOpts.exitOnlyDuringOffStates = true;
    modelOpts.nSims = nSims;
    gpurng(1, "ThreeFry"); %fastest gpu rng
    
    nSilentStates = contains(modelType, 'exit');
    if variableStateNumber
        nentries = ub(3);
        moffs = ub(4);
    else
        nentries = p0(3);
        moffs = p0(4);
    end
    nStates = nentries + moffs + 1 + nSilentStates;
    n_dls = length(DorsalFluoValues);
    modelOpts.r_vec = gpuArray(rand(1, nentries*nSims +...
        ( (moffs+1) * nSims * n_dls ) +...
        nSilentStates*((nStates-1) * nSims), 'single'));
end

model = struct;

if fun == "table"
    mdl = @(x, p) kineticFunForFits_table(x, p, modelOpts);
elseif fun== "sim"
    mdl = @(x, p) kineticFunForFits_sim_vec_gpu_customrnd(x, p, modelOpts);
end
% mdl = @(x, p) kineticFunForFits_sim(x, p, modelOpts);
model.modelfun   = mdl;  %use mcmcrun generated ssfun

if ~batchedAffinities
    model.ssfun = @(theta, data) sum( (data.ydata(:, 2:end) -...
        mdl(data.ydata(:, 1), theta) ).^2, 'omitnan');
end

results = [];
% [results,~,~,~]=mcmcrun(model,data,params,options,results);
% [results,~,~,~]=mcmcrun(model,data,params,options,results);
[results,chain,s2chain,~]=mcmcrun(model,data,params,options,results);

burnInTime = .25; %let's burn the first 25% of the chain just in case
chain = chain(round(burnInTime*nSteps):nSteps, :);
if ~isempty(s2chain)
    s2chain = s2chain(round(.25*nSteps):nSteps, :);
end



%%
close all force;

chainfig = figure(); clf
mcmcplot(chain,[],results,'chainpanel')

% Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.
%geweke is a measure of whether things converged between 0 and 1.
chainstats(chain,results)


try
    %ideally, these guys look like ellipses. if certain parameters give weird
    %shapes, it might mean those parameters should be removed from the model if
    %possible
    pairFig = figure; clf
    % mcmcplot(chain,[],results,'pairs', .5);
    mcmcplot(chain,[],results, 'pairs', 4);
    
    %Make the scatter plot more visible
    ax = gca;
    ax.Children(4).MarkerSize = 4;
    ax.Children(4).Color = [.4 .4 1];
end
%
% figure;
% out = mcmcpred(results,chain,[],data, mdl);
% mcmcpredplot(out);

%%
if isempty(results.mean)
    results.mean = mean(chain);
    warning('Results mean was empty. Using chain mean.')
end

theta_mean = getAllThetas(names, results);

% if modelType == "entryexit"
%     if ~fixKD
%         if ~batchedAffinities
%             
%             theta_mean = results.mean;
%         else
%             for k = 1:length(data)
%                 if ~variableStateNumber
%                     theta_mean{k} = [results.mean(1), results.mean(k+1),p0(4), p0(5),...
%                         results.mean(end-2:end)];
%                 else
%                     %                     theta_mean{k} = [results.mean(1), results.mean(k+1),...
%                     %                         results.
%                 end
%             end
%             
%         end
%     else
%         theta_mean = [results.mean(1), p0(2), results.mean(2:end)];
%     end
% elseif modelType == "entry"
%     if ~fixKD
%         if ~batchedAffinities
%             theta_mean = [results.mean(1:5) 0, results.mean(6)];
%         else
%             for k = 1:length(data)
%                 if ~variableStateNumber
%                     theta_mean{k} = [results.mean(1), results.mean(k+1),p0(4), p0(5),...
%                         results.mean(end-1), 0, results.mean(end)];
%                 else
%                     theta_mean{k} = [results.mean(1), results.mean(k+1),...
%                         results.mean(end-3), results.mean(end-2), results.mean(end-1), 0, results.mean(end)];
%                 end
%             end
%         end
%     else
%         theta_mean = [results.mean(1), p0(2), results.mean(2:4) 0, results.mean(5)];
%     end
% elseif modelType == "basic"
%     %c kd nentry moff pientry piexit tcycle
%     if ~fixKD
%         theta_mean = [results.mean(1:2) 0, results.mean(3) results.mean(5)];
%     else
%         theta_mean = [results.mean(1), p0(2), 0, results.mean(2), 1E10, 0, results.mean(3)];
%     end
% end

if ~iscell(theta_mean)
    yy = results.modelfun(DorsalFluoValues, theta_mean);
else
    for k = 1:length(theta_mean)
        yy{k} = results.modelfun(DorsalFluoValues, theta_mean{k});
    end
end

if ~iscell(yy)
    figure;
    tiledlayout('flow')
    nexttile;
    errorbar(DorsalFluoValues, nanmean(FractionsPerEmbryo, 1), nanstd(FractionsPerEmbryo, 1))
    hold on
    plot(DorsalFluoValues, yy(:, 1))
    
    ylabel('factive')
    xlabel('dl')
    legend('data', 'sim')
    nexttile;
    errorbar(DorsalFluoValues, nanmean(TimeOnsPerEmbryo, 1), nanstd(TimeOnsPerEmbryo, 1))
    
    hold on
    plot(DorsalFluoValues, yy(:, 2))
    
    ylabel('onset')
    xlabel('dl')
    legend('data', 'sim')
    
    % nexttile
    % scatter(sims.factive(:), sims.mfpts(:))
    % scatter(yy(:, 1), yy(:, 2))
    % xlim([0, 1.1])
    % ylim([0, 8.1])
    
    % nexttile;
    % plot(data.ydata(1, :), data.ydata(2, :), data.ydata(2, :))
    % ylabel('fraction')
    % xlabel('dl')
    % zlabel('onset')
else
    for k = 1:length(yy)
        figure;
        tiledlayout('flow')
        nexttile;
        errorbar(DorsalFluoValues, nanmean(FractionsPerEmbryo{k}, 1), nanstd(FractionsPerEmbryo{k}, 1))
        hold on
        plot(DorsalFluoValues, yy{k}(:, 1))
        
        ylabel('factive')
        xlabel('dl')
        legend('data', 'sim')
        nexttile;
        errorbar(DorsalFluoValues, nanmean(TimeOnsPerEmbryo{k}, 1), nanstd(TimeOnsPerEmbryo{k}, 1))
        
        hold on
        plot(DorsalFluoValues, yy{k}(:, 2))
        
        ylabel('onset')
        xlabel('dl')
        legend('data', 'sim')
        
    end
end



%%
try
    n = size(chain, 1);
    kd = results.theta(2)*ones(n, 1);
    nentries = results.theta(3)*ones(n, 1);
    moffs = results.theta(4)*ones(n, 1);
    piexits = results.theta(6)*ones(n, 1);
    theta = horzcat(chain(:, 1), kd, nentries, moffs, chain(:, 2), piexits);
    
    yyy = nan(length(DorsalFluoValues),length(results.mean),n);
    for k = 1:n
        yyy(:, :, k) = results.modelfun(DorsalFluoValues, theta(k, :));
    end
end

nexttile;
yyy_f = squeeze(yyy(:, 1, :));
yyy_o = squeeze(yyy(:, 2, :));

scatter(yyy_f(:), yyy_o(:));
xlim([0, 1.1])
ylim([0, 8.2])
%
% nexttile;
% % scatter(sims.factive(:), sims.mfpts(:))
% xlim([0, 1.1])
% ylim([0, 8.2])



disp('debug stop')
end