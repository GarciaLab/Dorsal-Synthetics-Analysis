function results = fitKineticModels(varargin)

% To do:
% 1. DONE. Fit all points, not just the averages. 
% 2. Fit multiple KDs simultaneously using batched fits. 
% 3. Port to R Stan
% 4. Re-do sims letting n states float and adjust parameter ranges,
% especially Dl. 

close all force;

wb = true;
nSteps = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats
exitOnlyDuringOffStates = true; %determines connectivity of the markov graph
modelType = "entry"; %choices- entryexit, entry, exit, basic
fun= "table"; %also 'sim'
t_cycle = 8;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

[~, dropboxfolder] = getDorsalFolders;

datPath = dropboxfolder + "\manuscript\window\basic\dataForFitting\archive\";
load(datPath + "DorsalFluoValues.mat", "DorsalFluoValues");
load(datPath + "FractionsPerEmbryo.mat", "FractionsPerEmbryo");
load(datPath + "TimeOnsPerEmbryo.mat", "TimeOnsPerEmbryo");


%needs to be nobs x ny
X = repmat(DorsalFluoValues, 1, max(size(FractionsPerEmbryo)));
data.ydata = [X; FractionsPerEmbryo(:)'; TimeOnsPerEmbryo(:)']';

%%
saveStr = modelType;
if exitOnlyDuringOffStates
    saveStr = saveStr + "_exitOnlyOnOffStates";
end
% sims = load(dropboxfolder +  "\simulations\archive\" + "tf_paramsearch_"+saveStr+"_.mat", 'params', 'factive', 'mfpts');
sims = load(dropboxfolder +  "\simulations\" + "tf_paramsearch_"+saveStr+"_.mat", 'params', 'factive', 'mfpts');

%%
rng(1, 'combRecursive') %matlab's fastest rng. ~2^200 period
% options.drscale = 2; % a high value (5) is important for multimodal parameter spaces. 
options.drscale = 5;
options.waitbar = wb; %the waitbar is rate limiting sometimes
options.nsimu = nSteps; %should be between 1E3 and 1E6
options.updatesigma = 1; %honestly don't know what this does
%
names = ["c", "kd" , "nentrystates", "moffstates", "pentry", "pexit"];
% p0 = [1E5, 1E5, 5, 5, .1, 3];
p0 = [10, 1E3, 5, 5, 1, 1];
lb = [1E-1, 1E2, 0, 1, 1E-2, 0];
ub = [1E6, 1E6, 12, 12, 1E3, 1E1];

params = cell(1, length(p0));
pri_mu = NaN; %default prior gaussian mean
pri_sig = Inf; %default prior gaussian variance
localflag = 0; %is this local to this dataset or shared amongst batches?

for k = 1:length(names)
    
    targetflag = 1;
    
    if contains(names(k), "states") && fun=="table"
        targetflag = 0;
    end
    
    if modelType == "entry" && names(k) == "pexit"
        targetflag = 0;
        p0(k) = 0;
    end
    
    params{1, k}= {names(k), p0(k), lb(k), ub(k), pri_mu, pri_sig, targetflag, localflag};
end


if fun == "table"
    modelOpts = struct('sims', sims);
    modelOpts.model = modelType;
else
    modelOpts.exitOnlyDuringOffStates = true;
    modelOpts.nSims = 1E4;
    gpurng(1, "ThreeFry"); %fastest gpu rng
end

model = struct;

if modelType == "entryexit"
    if fun == "table"
        mdl = @(x, p) kineticFunForFits_table(x, p, modelOpts);
    elseif fun== "sim"
        mdl = @(x, p) kineticFunForFits_sim_gpu(x, p, modelOpts);
    end
elseif modelType == "entry"
    mdl = @(x, p) entryAnalytical(x, p, t_cycle);
end
% mdl = @(x, p) kineticFunForFits_sim(x, p, modelOpts);
model.modelfun   = mdl;  %use mcmcrun generated ssfun
model.ssfun = @(theta, data) sum( (data.ydata(:, 2:end) - mdl(data.ydata(:, 1), theta) ).^2, 'omitnan');

results = [];
% [results,~,~,~]=mcmcrun(model,data,params,options,results);
[results,~,~,~]=mcmcrun(model,data,params,options,results);
[results,chain,s2chain,~]=mcmcrun(model,data,params,options,results);

burnInTime = .25; %let's burn the first 25% of the chain just in case
chain = chain(round(burnInTime*nSteps):nSteps, :);
if ~isempty(s2chain)
    s2chain = s2chain(round(.25*nSteps):nSteps, :);
end

chainfig = figure(); clf
mcmcplot(chain,[],results,'chainpanel')

% Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.
%geweke is a measure of whether things converged between 0 and 1.
chainstats(chain,results)


%ideally, these guys look like ellipses. if certain parameters give weird
%shapes, it might mean those parameters should be removed from the model if
%possible
pairFig = figure; clf
% mcmcplot(chain,[],results,'pairs', .5);
mcmcplot(chain,[],results, 'pairs', 4);

% 
% figure;
% out = mcmcpred(results,chain,[],data, mdl);
% mcmcpredplot(out);

%%
if modelType == "entryexit"
    theta_mean = [results.mean(1:2), 5, 5, results.mean(3:4)];
elseif modelType == "entry"
    theta_mean = [results.mean(1:2), 5, 5, results.mean(3), 0];
end
yy = mdl(DorsalFluoValues, theta_mean); 

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

nexttile
scatter(sims.factive(:), sims.mfpts(:))

% nexttile;
% plot(data.ydata(1, :), data.ydata(2, :), data.ydata(2, :))
% ylabel('fraction')
% xlabel('dl')
% zlabel('onset')



disp('debug stop')
end