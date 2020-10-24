function [results,chain,s2chain]  = fitstuff_mcmc2glob(varargin)

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')

expmnt = "affinities";
md = "simpleweak";
metric = "fraction";
lsq = true;
noOff = false;
nSimu = 1E4;
minKD = 200;
maxKD = 1E4;
minw = 1E-2; %1E-2
maxw = 1E1; %1E2
minR = 10;
maxR = 1E3;
displayFigures = true;
wb = true;
%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if expmnt == "phaff"
    lsq = false;
    noOff = true;
end

if strcmpi(expmnt, 'affinities')
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
    scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
elseif strcmpi(expmnt, 'phases')
    enhancers = {'1Dg11', '1Dg-5', '1Dg-8D'};
    scores = [0, -5, -8]';
elseif expmnt=="phaff"
    enhancers_aff =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
    scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
    enhancers_ph = {'1Dg11', '1Dg-5', '1Dg-8D'};
    positions = [0, -5, -8]';
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW', '1Dg-5', '1Dg-8D'};
end


%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = getXRange(enhancers, expmnt);

nSets = length(enhancers);
xo = {};
yo = {};
xs = {};
ys = {};
dsid = [];
T = [];
Y = [];
for k = 1:nSets
    cond = strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k});
    xo{k} = dorsalResultsDatabase.dorsalFluoBins(cond);
    if metric == "fraction"
        yo{k} = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
    elseif metric == "fluo"
        yo{k} = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond);
    end
    
    [xs{k}, ys{k}]= processVecs(xo{k}, yo{k}, xrange(k, :));
    dsid = [dsid; k*ones(size(xs{k}))];
    
    if isnan(xrange(k, 1))
        xrange(k, 1) = min(xo{k});
    end
    if isnan(xrange(k, 2))
        xrange(k, 2) = max(xo{k});
    end
    
    assert(~isempty(xs{k}));
    
    T = [T; xs{k}];
    Y = [Y; ys{k}];
end

data.ydata = [T, Y];
data.dsid = dsid;
data.X =  [T dsid];

%%
% Refine the first guess for the parameters with fminseacrh and calculate residual variance as an estimate of the model error variance.


%rate, kd, hill, y offset
y_max = nanmax(Y(:));
x_max = max(T);

[p0, lb, ub] = getInits(expmnt, md, metric ,x_max, y_max, nSets,...
    'minKD', minKD, 'maxKD', maxKD, 'minw', minw, 'maxw', maxw, 'minR', minR, 'maxR', maxR);


if lsq
    [k0, mse] = globfit2('expmnt',expmnt, 'metric', metric, 'md', md, 'maxKD', maxKD, 'displayFigures', displayFigures);
else
    k0 = p0;
end


if noOff && metric=="fluo"
    p0 = p0(1:end-1);
    lb = lb(1:end-1);
    ub = ub(1:end-1);
    k0 = k0(1:end-1);
end

%
params = cell(1, 3);
for i = 1:length(k0)
    params{1, i} = {['k', num2str(i)],k0(i), lb(i), ub(i)};
end

model = struct;

if expmnt == "affinities" && md=="simpleweak" && metric=="fraction"
    model.ssfun = @funss_simpleweakfraction;
    mdl = @(x, p)subfun_simplebinding_weak_fraction_std2(x, p);
elseif expmnt == "affinities" && md=="simpleweak" && metric=="fluo"
    if noOff
        model.ssfun = @funss_simpleweakfluo_nooff;
        mdl = @(x, p)subfun_simplebinding_weak_fluo_std2_nooff(x, p);
    else
        model.ssfun = @funss_simpleweakfluo;
        mdl = @(x, p)subfun_simplebinding_weak_fluo_std2(x, p);
    end
elseif expmnt == "phases" && md=="simpleweak" && metric=="fraction"
    model.ssfun = @funss_simpleweakfraction_phases;
    mdl = @(x, p)subfun_simplebinding_weak_fraction_std2_phases(x, p);
elseif expmnt == "phases" && md=="simpleweak" && metric=="fluo"
    if noOff
        model.ssfun = @funss_simpleweakfluo_phases_nooff; %#ok<*UNRCH>
        mdl = @(x, p)subfun_simplebinding_weak_fluo_std2_phases_nooff(x, p);
    else
        model.ssfun = @funss_simpleweakfluo_phases;
        mdl = @(x, p)subfun_simplebinding_weak_fluo_std2_phases(x, p);
    end
elseif expmnt == "phaff" && md=="simpleweak" && metric=="fraction"
    model.ssfun = @funss_simpleweakfraction_phaff;
    mdl = @(x, p)subfun_simplebinding_weak_fraction_std2_phaff(x, p);
elseif expmnt == "phaff" && md=="simpleweak" && metric=="fluo"
    model.ssfun = @funss_simpleweakfluo_phaff_nooff;
    mdl = @(x, p)subfun_simplebinding_weak_fluo_std2_phaff_nooff(x, p);
end

if lsq
    model.sigma2 = mse;
end

options.drscale = 5;
options.waitbar = wb;
options.nsimu = nSimu;
options.updatesigma = 1;

rng(1,'twister'); %set the rng seed so we get the same results every run of this function

if lsq
    [results,chain,s2chain] = mcmcrun(model,data,params,options);
else
    results = [];
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
end

if displayFigures
    burnInTime = .25; %let's burn the first 25% of the chain
    chain = chain(round(burnInTime*nSimu):nSimu, :);
    if ~isempty(s2chain)
        s2chain = s2chain(round(.25*nSimu):nSimu, :);
    end
    
    chainfig = figure(); clf
    mcmcplot(chain,[],results,'chainpanel')
    
    % Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.
    chainstats(chain,results)
    
    
    
    %ideally, these guys look like ellipses. if certain parameters give weird
    %shapes, it might mean those parameters should be removed from the model if
    %possible
    pairFig = figure; clf
    mcmcplot(chain,[],results,'pairs', .5);
    %
    
    
    figure(4); clf
    til = tiledlayout(1, nSets);
    dsid2 = [];
    
    xx = (0:10:max(data.X(:,1)))';
    
    for k = 1:nSets
        dsid2 = [dsid2; k*ones(length(xx), 1)];
    end
    
    
    out = mcmcpred(results,chain,[],repmat(xx, nSets, 1), mdl);
    % mcmcpredplot(out);
    nn = (size(out.predlims{1}{1},1) + 1) / 2;
    plimi = out.predlims{1};
    yl = plimi{1}(1,:);
    yu = plimi{1}(2*nn-1,:);
    
    km = mean(chain);
    ks = std(chain);
    X2 = [repmat(xx, nSets, 1), dsid2];
    
    
    yf = plimi{1}(nn,:);
    
    for i = 1:nSets
        
        %     yy = yfit2(X2(:, 2)==i);
        yy = yf(X2(:, 2)==i);
        yyl = yl(X2(:, 2)==i);
        yyu = yu(X2(:, 2)==i);
        nexttile;
        fillyy(xx,yyl,yyu,[0.9 0.9 0.9]);
        hold on
        plot(xo{i}, yo{i}, 'o-', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'r');
        plot( xo{i}(xo{i} >= xrange(i, 1) & xo{i} <=xrange(i, 2) ), yo{i}(xo{i} >=  xrange(i, 1) & xo{i} <= xrange(i, 2)),'o-','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Color', 'r');
        plot(xx,yy,'-k')
        xlim([0,3500])
        
        if expmnt == "affinities"
            vartheta = 'KD = ';
            consttheta = ' \omega'' = ';
        elseif expmnt == "phases"
            consttheta = 'KD = ';
            vartheta = ' \omega'' = ';
        end
        
        
        if expmnt == "phaff"
            naff = 7;
            nph = 3;
            if i <= naff
                
                titleCell = {enhancers{i}, [' \omega'' = ', num2str(round2(km(1))), ' \pm ', num2str(round2(ks(1)))],...
                    [ 'KD = ', num2str(round2(km(nph+i))), ' \pm ', num2str(round2(ks(nph+i)))]};
            elseif i > naff
                titleCell = {enhancers{i}, [' \omega'' = ', num2str(round2(km(i - (naff-1) ))), ' \pm ', num2str(round2(ks(i - (naff-1))))],...
                    [ 'KD = ', num2str(round2(km(nph+1))), ' \pm ', num2str(round2(ks(nph+1)))]};
            end
        else
            titleCell = {enhancers{i}, [vartheta, num2str(round2(km(i+1))), ' \pm ', num2str(round2(ks(i+1)))],...
                [consttheta, num2str(round2(km(1))), ' \pm ', num2str(round2(ks(1)))]};
        end
        
        if metric=="fraction"
            ylim([0, 1])
        elseif metric=="fluo"
            titleCell = [titleCell,...
                [' amp = ', num2str(round2(km(nSets + 1))), ' \pm ', num2str(round2(ks(nSets + 1)))] ];
            if ~noOff
                titleCell = [titleCell, [' off = ', num2str(round2(km(end))), ' \pm ', num2str(round2(ks(end)))] ];
            end
            ylim([0, y_max])
        end
        title(titleCell);
    end
    %     title(til, 'global fit with mcmc')
    
    
    figure;
    if expmnt == "phaff"
        
        tilo = tiledlayout('flow');
        nexttile;
        errorbar(scores, km(nph+1:nph+naff), ks(nph+1:nph+naff));
        ylabel('KD (au)')
        xlabel('affinity (Patser score)')
        
        nexttile;
        errorbar(positions, km(1:nph), ks(1:nph));
        ylabel('w'' (au)')
        xlabel('position (bp)')
        xlim([-10, 2])
    else
        
        errorbar(scores, km(2:nSets+1), ks(2:nSets+1));
        if expmnt == "affinities"
            ylabel('KD (au)')
            xlabel('affinity (Patser score)')
        elseif expmnt == "phases"
            ylabel('w'' (au)')
            xlabel('position (bp)')
            xlim([-10, 2])
        end
        
    end
    
    
    covfig = figure;
    cv = @(x, y) sqrt(abs(x)) ./ sqrt((y'*y));
    imagesc(cv(results.cov, results.mean));
    colorbar;     
    ylabel('parameter 1')
    xlabel('parameter 2')
    title('Covariance matrix of the parameters');
   
    colormap(viridis);
    
    
    
end

end
%%
function ss = funss_hill(params, data)
% sum-of-squares


x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
nSets = max(dsid);
params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
n = params(nSets + 2);
offset = params(nSets + 3);

Amodel = amplitude.*(((x./KD(dsid)).^n)./(1+((x)./KD(dsid)).^n))+offset;

ss = sum(( data.ydata(:,2)-Amodel).^2);
end

%%
function ss = funss_simpleweak(params, data)
% sum-of-squares

x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
omegaDP = params(nSets + 2);
offset = params(nSets + 3);

Amodel = amplitude.*(((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP)))+offset;

ss = sum((data.ydata(:,2)-Amodel).^2);
end

function yfit = subfun_simple_weak(params,X)

%simplebinding in the weak promoter limit.

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
omegaDP = params(nSets + 2);
offset = params(nSets + 3);

yfit = amplitude.*(((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+((x./KD(dsid)).*omegaDP)))+offset;

end

%% simple weak fraction

function ss = funss_simpleweakfraction(params, data)
% sum-of-squares

x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
nSets = max(dsid);
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:nSets+1)';
Amodel = (((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid)) + ((x./KD(dsid)).*omegaDP)));
ss = sum((data.ydata(:,2)-Amodel).^2);
end

function yfit = subfun_simplebinding_weak_fraction_std2(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:max(dsid)+1)';
yfit = (((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP)));

end

%% simple weak fluo


function ss = funss_simpleweakfluo(params, data)
% sum-of-squares

x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
Amodel = amp.*(((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP))) + offset;
ss = sum((data.ydata(:,2)-Amodel).^2);
end


function yfit = subfun_simplebinding_weak_fluo_std2(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
yfit = amp.*(((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP))) + offset;

end

function ss = funss_simpleweakfluo_nooff(params, data)
% sum-of-squares

x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
Amodel = amp.*(((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid)) + ((x./KD(dsid)).*omegaDP)));
ss = sum((data.ydata(:,2)-Amodel).^2);
end


function yfit = subfun_simplebinding_weak_fluo_std2_nooff(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
yfit = amp.*(((x./KD(dsid)).*omegaDP)./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP)));

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%phases

function ss = funss_simpleweakfraction_phases(params, data)
% sum-of-squares

x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
nSets = max(dsid);
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:nSets+1)';
Amodel = (((x./KD).*omegaDP(dsid))./(1+ (x./KD)+ ((x./KD).*omegaDP(dsid))));
ss = sum((data.ydata(:,2)-Amodel).^2);
end

function yfit = subfun_simplebinding_weak_fraction_std2_phases(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
yfit = (((x./KD).*omegaDP(dsid))./(1+ (x./KD)+((x./KD).*omegaDP(dsid))));

end

%% simple weak fluo


function ss = funss_simpleweakfluo_phases(params, data)
% sum-of-squares

x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
Amodel = amp.*(((x./KD).*omegaDP(dsid))./(1+ (x./KD)+ ((x./KD).*omegaDP(dsid)))) + offset;
ss = sum((data.ydata(:,2)-Amodel).^2);
end


function yfit = subfun_simplebinding_weak_fluo_std2_phases(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
offset = params(max(dsid)+3);
yfit = amp.*(((x./KD).*omegaDP(dsid))./(1+ (x./KD)+ ((x./KD).*omegaDP(dsid)))) + offset;

end


function ss = funss_simpleweakfluo_phases_nooff(params, data)
% sum-of-squares

x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
Amodel = amp.*(((x./KD).*omegaDP(dsid))./(1+(x./KD)+ ((x./KD).*omegaDP(dsid))));
ss = sum((data.ydata(:,2)-Amodel).^2);
end


function yfit = subfun_simplebinding_weak_fluo_std2_phases_nooff(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
KD = params(1);
omegaDP = params(2:max(dsid)+1)';
amp = params(max(dsid)+2);
yfit = amp.*(((x./KD).*omegaDP(dsid))./(1+ (x./KD)+ ((x./KD).*omegaDP(dsid))));

end


function ss = funss_simpleweakfluo_phaff_nooff(params, data)
% sum-of-squares
nph = 3;
naff = 7;

x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
KD = [params(1:naff), params(1)*ones(1, nph-1)]';
amp = params(nph+naff+1);
Amodel = amp.*(((x./KD(dsid)).*omegaDP(dsid))./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP(dsid))));
ss = sum( (data.ydata(:,2)-Amodel).^2 );
end


function yfit = subfun_simplebinding_weak_fluo_std2_phaff_nooff(x, params)
%simplebinding in the weak promoter limit.
%note that this could probably be cleaner if dsid was a 2d matrix and
%params were referenced in that matrix like data set 1 has kd1, omega1,
%etc.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

nph = 3;
naff = 7;
dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
KD = [params(1:naff), params(1)*ones(1, nph-1)]';
amp = params(nph+naff+1);
yfit = amp.*(((x./KD(dsid)).*omegaDP(dsid))./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP(dsid))));

end


function ss = funss_simpleweakfraction_phaff(params, data)
% sum-of-squares

nph = 3;
naff = 7;
x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
KD = [params(1:naff), params(1)*ones(1, nph-1)]';
Amodel = (((x./KD(dsid)).*omegaDP(dsid))./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP(dsid))));
ss = sum((data.ydata(:,2)-Amodel).^2);
end

function yfit = subfun_simplebinding_weak_fraction_std2_phaff(x, params)
%simplebinding in the weak promoter limit.

nph = 3;
naff = 7;

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
KD = [params(1:naff), params(1)*ones(1, nph-1)]';
yfit = (((x./KD(dsid)).*omegaDP(dsid))./(1+ (x./KD(dsid))+ ((x./KD(dsid)).*omegaDP(dsid))));

end