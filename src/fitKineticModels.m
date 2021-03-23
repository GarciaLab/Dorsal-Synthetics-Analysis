function results = fitKineticModels(varargin)

close all force;

wb = true;
nSims = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

datPath = "C:\Users\owner\Dropbox\DorsalSyntheticsDropbox\manuscript\window\basic\dataForFitting\";
load(datPath + "DorsalFluoValues.mat", "DorsalFluoValues");
load(datPath + "FractionsPerEmbryo.mat", "FractionsPerEmbryo");
load(datPath + "TimeOnsPerEmbryo.mat", "TimeOnsPerEmbryo");
FractionsPerEmbryo = nanmean(FractionsPerEmbryo,  1);
TimeOnsPerEmbryo = nanmean(TimeOnsPerEmbryo, 1);

%needs to be nobs x ny
data.ydata = [DorsalFluoValues; FractionsPerEmbryo; TimeOnsPerEmbryo]';

%%
model = "entryexit";
exitOnlyDuringOffStates = true;
[~, dropboxfolder] = getDorsalFolders;
saveStr = model;
if exitOnlyDuringOffStates
    saveStr = saveStr + "_exitOnlyOnOffStates";
end
sims = load(dropboxfolder + "\" + "tf_paramsearch_"+saveStr+"_.mat", 'params', 'factive', 'mfpts');
%
% dlbins = linspace(0,4500,20);
% simsbinneddl = BinData(sims.params.dls, dlbins);
% for k = 1:length(dlbins)
%     sims.factive_binned(k) =
%     sims.onset_binned(k) =
% end

%%
rng(1,'twister'); %set the rng seed so we get the same results every run of this function
options.drscale = 5; % a high value (5) is important for multimodal parameter spaces
options.waitbar = wb; %the waitbar is rate limiting sometimes
options.nsimu = nSims; %should be between 1E3 and 1E6
options.updatesigma = 1; %honestly don't know what this does
%
names = ["c", "kd" , "nentrystates", "moffstates", "pentry", "pexit"];
% p0 = [1E5, 1E5, 5, 5, .1, 3];
p0 = [10, 1E3, 5, 5, 1, 1];
lb = [1E-6, 1E2, 1, 1, 1E-2, 0];
ub = [1E6, 1E6, 12, 12, 1E3, 1E1];

params = cell(1, length(p0));
pri_mu = NaN; %default prior gaussian mean
pri_sig = Inf; %default prior gaussian variance
localflag = 0; %is this local to this dataset or shared amongst batches?


for k = 1:length(names)
    
    if contains(names(k), "states")
        targetflag = 0;
    else
        targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.
    end
    
    params{1, k}= {names(k), p0(k), lb(k), ub(k), pri_mu, pri_sig, targetflag, localflag};
end


modelOpts = struct('sims', sims);

model = struct;
mdl = @(x, p) kineticFunForFits(x, p, modelOpts);
model.modelfun   = mdl;  %use mcmcrun generated ssfun
model.ssfun = @(theta, data) sum( (data.ydata(:, 2:end) - mdl(data.ydata(:, 1), theta) ).^2, 'omitnan');

results = [];
[results,~,~,~]=mcmcrun(model,data,params,options,results);
[results,~,~,~]=mcmcrun(model,data,params,options,results);
[results,chain,s2chain,~]=mcmcrun(model,data,params,options,results);

burnInTime = .25; %let's burn the first 25% of the chain just in case
chain = chain(round(burnInTime*nSims):nSims, :);
if ~isempty(s2chain)
    s2chain = s2chain(round(.25*nSims):nSims, :);
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
mcmcplot(chain,[],results, 'pairs');

% 
% figure;
% out = mcmcpred(results,chain,[],data, mdl);
% mcmcpredplot(out);

figure;
scatter(sims.factive(:), sims.mfpts(:))


theta_mean = [results.mean(1:2), 5, 5, results.mean(3:4)];
yy = mdl(data.ydata(:, 1), theta_mean); 

figure;
tiledlayout('flow')
nexttile;
plot(data.ydata(:, 1), data.ydata(:, 2))
hold on
plot(data.ydata(:, 1), yy(:, 1))
ylabel('factive')
xlabel('dl')
legend('data', 'sim')
nexttile;
plot(data.ydata(:, 1), data.ydata(:, 3))
hold on
plot(data.ydata(:, 1), yy(:, 2))
ylabel('onset')
xlabel('dl')
legend('data', 'sim')
% nexttile;
% plot(data.ydata(1, :), data.ydata(2, :), data.ydata(2, :))
% ylabel('fraction')
% xlabel('dl')
% zlabel('onset')



disp('debug stop')
end