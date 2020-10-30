pdets = [.01:.01:1];
nTrials = 2;
for k = 1:length(pdets)
    [fractive(k), take3pactive(k)] = pactiveSim('pdetect', pdets(k), 'nTrials', nTrials);
end
figure; 

tiledlayout('flow')
nexttile;
plot(pdets, fractive, pdets, take3pactive);
xlabel('true p detect')
ylabel('pactive')
legend('true fraction active', 'inferred p active');
ylim([0, 1])
nexttile
plot(pdets, 1 - take3pactive ./ fractive);
xlabel('true p detect')
ylabel('error from estimate')

figure;
pdets = [.01:.01:1];
nTrials = [1:10];
palette = parula(length(nTrials));

for j = 1:length(nTrials)
    for k = 1:length(pdets)
        [fractive(k), take3pactive(k)] = pactiveSim('pdetect', pdets(k), 'nTrials', nTrials(j));
    end
    plot(pdets, 1 - take3pactive ./ fractive, 'Color', palette(j,:));
    hold on
end
xlabel('true p detect')
ylabel('error from estimate')
cbh = colorbar;
set(cbh,'YTickLabel', nTrials)
ylabel(cbh, 'nTrials');
