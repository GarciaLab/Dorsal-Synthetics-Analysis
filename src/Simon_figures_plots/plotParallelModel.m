function plotParallelModel(dorsalVals,theta)

%theta = [c,kd,nSwitches,DlIndependentK,tcycle]
% figure
% hold on
% ParallelModel(dorsalVals,theta,[])
% hold off

%%
close all
KdValues = logspace(2,3.5,7);
Palette = flip(viridis(length(KdValues)));
figure
hold on
for k = 1:length(KdValues)
    Color = Palette(k,:);
    kd = KdValues(k);
    theta(2) = kd;
    fraction_onset = ParallelModel(dorsalVals,theta,[]);
    plot(dorsalVals,fraction_onset(:,1),'Color',Color,'LineWidth',2)
    %plot(dorsalVals,fraction_onset(:,2),'Color',Color,'LineWidth',2)
end
hold off
ylim([0 1])
xlim([1 dorsalVals(end)])
set(gca,'XScale','log')

