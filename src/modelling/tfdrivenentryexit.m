
close all force;


dmax = 5000;
nPlots = 10;
dl = 500;
dls = logspace(1, log10(dmax)); %aus
t = linspace(0, 10)'; %mins
kd = 500;
kds = logspace(2, 4, nPlots);
cs = logspace(1, 3, nPlots);
pi1s = logspace(-2, 1, nPlots);
pi2s = logspace(-2, 1, nPlots);

R = 500;
c = 2;

t_cycle = 10; %min

nSteps = 6;
nSims = 1E3;
nOffStates = 5;
nEntryStates = 5;
firstoffstate = nEntryStates+1;
onstate = nEntryStates + nOffStates+1;
silentstate = onstate+1;
nStates = nEntryStates + nOffStates + 1 + 1;
occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );
pi0 = 1; %min
pi1 = 1; %min
pi2 = 10; %min

dls = 700; %au
mfpts = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));

% cmap = colormap(viridis(nPlots));
cmap = colormap(parula(nPlots));
% cs = 10;
% pi1s = [1, .1];
pi2s = 1;
pi1 = 1/3;

for n = 1:length(pi2s)
    for m = 1:length(cs)
        for k = 1:length(pi1s)
            for i = 1:length(dls)
                for j = 1:length(kds)
                    pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1
                    pi1 = pi1s(k); %min-1
                    pi2 = pi2s(n);
                    
                    [fpt_on_observed, factive_temp] = averagePaths(nSims, nSteps, pi0, pi1, pi2,onstate, silentstate, t_cycle, firstoffstate);
                    
                    mfpts(i, j, k, m, n) = mean(fpt_on_observed);
                    factive(i,j,k,m,n) = factive_temp;
                end
            end
        end
    end
end



figure;
tiledlayout('flow')
for j = 1:length(cs)
    nexttile;
    for k = 1:length(pi1s)
        %     plot(kds, mfpts(40, :, 1));
        plot(kds, mfpts(1, :, k, j, 1), 'LineWidth', 2, 'Color', cmap(k, :));
        hold on
    end
    xlabel('K_D (au)')
    ylabel('mean time to turn on (min)')
    %     leg =legend(num2str(round2(pi1s')));
    %     title(leg, '\pi_1');
    %     set(gca, 'XScale', 'log')
    title(['c = ', num2str(cs(j))])
    
end


%%
figure;
% tiledlayout('flow')
for j = 1:length(cs)
    %     nexttile;
    %     for k = 1:length(pi1s)
    %     plot(kds, mfpts(40, :, 1));
    plot(kds, mfpts(1, :, 6, j, 1), 'LineWidth', 2, 'Color', cmap(j, :));
    hold on
    %     end
    xlabel('K_D (au)')
    ylabel('mean time to turn on (min)')
    %     leg =legend(num2str(round2(pi1s')));
    %     title(leg, '\pi_1');
    %     set(gca, 'XScale', 'log')
    %     title(['c = ', num2str(cs(j))])
    
end
leg =legend(num2str(round2(cs')));
title(leg, 'c (min-1)');
title({'pientry=1', 'pisilent=.3', 'Dl=700'})
%%
%
% figure;
% plot(pi1s, mfpts(:, 5));
% xlabel('super off transition rate (min-1)')
% ylabel('mean time to turn on (min)')

%
% figure;
% plot(dls, mfpts(:, 5));
% xlabel('[Dl] (au)')
% ylabel('mean time to turn on (min)')

%%


dmax = 5000;
nPlots = 10;
dl = 500;
t = linspace(0, 10)'; %mins
dls = logspace(1, log10(dmax)); %aus
kd = 500;
kds = logspace(2, 4, nPlots);
cs = logspace(1, 3, nPlots);
pi1s = logspace(-2, 1, nPlots);
pi2s = logspace(-2, 1, nPlots);

R = 500;
c = 2;

t_cycle = 10; %min

nSteps = 6;
nSims = 1E3;
nOffStates = 5;
nEntryStates = 5;
firstoffstate = nEntryStates+1;
onstate = nEntryStates + nOffStates+1;
silentstate = onstate+1;
nStates = nEntryStates + nOffStates + 1 + 1;
occupancy = @(d, kd) ( (d./kd) ./ (1 + d./kd) );
pi0 = 1; %min
pi1 = 1; %min
pi2 = 10; %min


% cmap = colormap(viridis(nPlots));
cmap = colormap(parula(nPlots));
% cs = 10;
% pi1s = [1, .1];
% cs = 100;
kds = 1000;
pi2s = 1;
pi1s = 1/3;

mfpts = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));

for n = 1:length(pi2s)
    for m = 1:length(cs)
        for k = 1:length(pi1s)
            for i = 1:length(dls)
                for j = 1:length(kds)
                    pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1
                    pi1 = pi1s(k); %min-1
                    pi2 = pi2s(n);
                    
                    [fpt_on_observed, factive_temp] = averagePaths(nSims, nSteps, pi0, pi1, pi2,onstate, silentstate, t_cycle, firstoffstate);
                    
                    mfpts(i, j, k, m, n) = mean(fpt_on_observed);
                    factive(i,j,k,m,n) = factive_temp;
                end
            end
        end
    end
end



figure;
tiledlayout('flow')

for m = 1:length(cs)
    nexttile;
    
    yyaxis left
    plot(dls, factive(:,1,1,m,1), 'LineWidth', 2)
    ylabel('factive)')
    
    yyaxis right
    plot(dls, mfpts(:,1,1,m,1), 'LineWidth', 2)
    ylabel('t on (min)')
    
    xlabel('dl (au)')
    title(['c = ', num2str(cs(m))])
    
end




