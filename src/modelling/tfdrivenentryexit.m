
% close all force;

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
pi0 = 1; %min-1
pi1 = 1; %min-1
pi2 = 10; %min-1

dls = 700; %au

% cmap = colormap(viridis(nPlots));
cmap = colormap(parula(nPlots));
% cs = 10;
% pi1s = [1, .1];
pi2s = 1;
pi1 = 1/3;

onlyEntry = true;
if onlyEntry == true
%     dls = linspace(1, log10(dmax)); %aus
    dls = linspace(1, dmax, 20);
    kds = logspace(2, 4, nPlots);
%     cs = logspace(1, 3, nPlots);
    cs = 1;
%     pi1s = logspace(-2, 1, nPlots);
    pi1s = 0;
    pi2s = logspace(-2, 1, nPlots);
%     pi2s = 100;
end

mfpts = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));
factive = nan(length(dls), length(kds), length(pi1s), length(cs), length(pi2s));

for n = 1:length(pi2s)
    for m = 1:length(cs)
        for k = 1:length(pi1s)
            for i = 1:length(dls)
                for j = 1:length(kds)
                    pi0 = cs(m).*occupancy(dls(i), kds(j)); %min-1
                    pi1 = pi1s(k); %min-1
                    pi2 = pi2s(n); %min-1
                    
                    [fpt_on_observed, factive_temp] = averagePaths_entryexit(nSims, nSteps, pi0, pi1, pi2,onstate, silentstate, t_cycle, firstoffstate);
                    

                    mfpts(i, j, k, m, n) = nanmean(fpt_on_observed);

                    factive(i,j,k,m,n) = factive_temp;
                end
            end
        end
    end
end

dt = mfpts(:, 10, :, :, :) - mfpts(:, 4, :, :, :); %kd(10)=10k, kd(4)=400

dt = repmat(dt, [1 length(kds) 1 1 1]);



%%
% 
% figure;
% tiledlayout('flow')
% for j = 1:length(cs)
%     nexttile;
%     for k = 1:length(pi1s)
%         %     plot(kds, mfpts(40, :, 1));
%         plot(kds, mfpts(1, :, k, j, 1), 'LineWidth', 2, 'Color', cmap(k, :));
%         hold on
%     end
%     xlabel('K_D (au)')
%     ylabel('mean time to turn on (min)')
%     %     leg =legend(num2str(round2(pi1s')));
%     %     title(leg, '\pi_1');
%     %     set(gca, 'XScale', 'log')
%     title(['c = ', num2str(cs(j))])
%     
% end
% 
% 
% %%
% figure;
% % tiledlayout('flow')
% for j = 1:length(cs)
%     %     nexttile;
%     %     for k = 1:length(pi1s)
%     %     plot(kds, mfpts(40, :, 1));
%     plot(kds, mfpts(1, :, 6, j, 1), 'LineWidth', 2, 'Color', cmap(j, :));
%     hold on
%     %     end
%     xlabel('K_D (au)')
%     ylabel('mean time to turn on (min)')
%     %     leg =legend(num2str(round2(pi1s')));
%     %     title(leg, '\pi_1');
%     %     set(gca, 'XScale', 'log')
%     %     title(['c = ', num2str(cs(j))])
%     
% end
% leg =legend(num2str(round2(cs')));
% title(leg, 'c (min-1)');
% title({'pientry=1', 'pisilent=.3', 'Dl=700'})
% %%
% %
% % figure;
% % plot(pi1s, mfpts(:, 5));
% % xlabel('super off transition rate (min-1)')
% % ylabel('mean time to turn on (min)')
% 
% %
% % figure;
% % plot(dls, mfpts(:, 5));
% % xlabel('[Dl] (au)')
% % ylabel('mean time to turn on (min)')
% 
% 
% figure;
% tiledlayout('flow')
% 
% for m = 1:length(cs)
%     nexttile;
%     
%     yyaxis left
%     plot(dls, factive(:,1,1,m,1), 'LineWidth', 2)
%     ylabel('factive)')
%     
%     yyaxis right
%     plot(dls, mfpts(:,1,1,m,1), 'LineWidth', 2)
%     ylabel('t on (min)')
%     
%     xlabel('dl (au)')
%     title(['c = ', num2str(cs(m))])
%     
% end
% 
% 
% figure
% scatter(factive(:), mfpts(:))
% hold on
% rectangle('Position',[.1 3.5 .9 3.5], 'LineWidth', 2, 'EdgeColor', 'r', 'Curvature',0.2) %lower left corner x&y, width, height
% xlabel('factive')
% ylabel('mean turn on (min)')
% 
% %%
% 
% figure; 
% tiledlayout('flow')
% for k = 1:length(pi2s)
% %dl kd pi1 c pi2
%     nexttile;
%     factivetemp = factive(:,:, :, :, k);
%     mfptstemp = mfpts(:,:, :, :, k);
% 
%     scatter(factivetemp(:), mfptstemp(:))
%     hold on
%     fmin = .1;
%     fmax = 10;
%     tmin = 3.5;
%     tmax = 7;
%     rectangle('Position',[.1 3.5 .9 3.5], 'LineWidth', 2, 'EdgeColor', 'r', 'Curvature',0.2) %lower left corner x&y, width, height
%     xlabel('factive')
%     ylabel('mean turn on (min)')
%     
%     [in,on] = inpolygon(factivetemp(:),mfptstemp(:),[.1 .1 1 1],[3.5 7 7 3.5]);
%     title(pi2s(k))
% 
% end
% 
% 
% figure; 
% tiledlayout('flow')
% for k = 1:length(pi1s)
% %dl kd pi1 c pi2
%     nexttile;
%     factivetemp = factive(:,:,k, :, :);
%     mfptstemp = mfpts(:,:, k, :, :);
% 
%     scatter(factivetemp(:), mfptstemp(:))
%     hold on
%     fmin = .1;
%     fmax = 10;
%     tmin = 3.5;
%     tmax = 7;
%     rectangle('Position',[.1 3.5 .9 3.5], 'LineWidth', 2, 'EdgeColor', 'r', 'Curvature',0.2) %lower left corner x&y, width, height
%     xlabel('factive')
%     ylabel('mean turn on (min)')
%     
%     [in,on] = inpolygon(factivetemp(:),mfptstemp(:),[.1 .1 1 1],[3.5 7 7 3.5]);
%     title(pi1s(k))
% 
% end
% 
% figure; 
% tiledlayout('flow')
% for k = 1:length(cs)
% %dl kd pi1 c pi2
%     nexttile;
%     factivetemp = factive(:,:,:, k, :);
%     mfptstemp = mfpts(:,:, :, k, :);
% 
%     scatter(factivetemp(:), mfptstemp(:))
%     hold on
%     fmin = .1;
%     fmax = 10;
%     tmin = 3.5;
%     tmax = 7;
%     rectangle('Position',[.1 3.5 .9 3.5], 'LineWidth', 2, 'EdgeColor', 'r', 'Curvature',0.2) %lower left corner x&y, width, height
%     xlabel('factive')
%     ylabel('mean turn on (min)')
%     
%     [in,on] = inpolygon(factivetemp(:),mfptstemp(:),[.1 .1 1 1],[3.5 7 7 3.5]);
%     title(cs(k))
% 
% end
% 
% 
% figure; 
% tiledlayout('flow')
% for k = 1:length(pi2s)
% %dl kd pi1 c pi2
%     nexttile;
%     factivetemp = factive(:,:, :, :, k);
%     mfptstemp = mfpts(:,:, :, :, k);
%     scatter(factivetemp(:), mfptstemp(:))
%     hold on
%     fmin = .1;
%     fmax = 10;
%     tmin = 3.5;
%     tmax = 7;
%     rectangle('Position',[.1 3.5 .9 3.5], 'LineWidth', 2, 'EdgeColor', 'r', 'Curvature',0.2) %lower left corner x&y, width, height
%     xlabel('factive')
%     ylabel('mean turn on (min)')
%     
%     [in,on] = inpolygon(factivetemp(:),mfptstemp(:),[.1 .1 1 1],[3.5 7 7 3.5]);
%     title(pi2s(k))
% 
% end


% %%
% figure
% 
% scatter3(factive(:), mfpts(:), dt(:))
% hold on
% fmin = .1;
% fmax = 1;
% tmin = 3.5;
% tmax = 7;
% deltatmin = 0;
% deltatmax = 3;
% % rectangle('Position',[.1 3.5 .9 3.5], 'LineWidth', 2, 'EdgeColor', 'r', 'Curvature',0.2) %lower left corner x&y, width, height
% % plot3( [.1 .1 1 1 .1], [3.5 7 7 3.5 3.5], [0 0 3 3 0], 'LineWidth', 3, 'Color', 'r' )
% [X,Y,Z] = ellipsoid(.55, 5.25, 1.5, .45, 1.75, 1.5);
% s = surf(X,Y,Z,'FaceAlpha',0.5, 'FaceColor', 'r');
% s.EdgeColor = 'none';
% 
% xlabel('factive')
% ylabel('mean turn on (min)')
% zlabel('delta t')
% 
% [in,on] = inpolygon(factivetemp(:),mfptstemp(:),[.1 .1 1 1],[3.5 7 7 3.5]);

[~, dropboxfolder] = getDorsalFolders;

if onlyEntry
    save([dropboxfolder, 'tfentry_paramsearch.mat'])
else
    save([dropboxfolder, 'tfentryexit_paramsearch.mat'])
end
% % load([dropboxfolder, 'tfentryexit_paramsearch.mat'])
% load([dropboxfolder, 'tfentry_paramsearch.mat'])

try
    plotTFDrivenParams(factive, dt, mfpts, 'nPoints', 2E4)
catch
    plotTFDrivenParams(factive, dt, mfpts);
end