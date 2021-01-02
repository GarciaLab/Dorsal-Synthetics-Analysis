function Comparing_enhancers

%% grab results
clear all
dorsalResultsDatabase = createDorsalResultsDatabase;

dlfluobins_1 = dorsalResultsDatabase.dorsalFluoBins; dlfluobins_1 = dlfluobins_1(1:37);
%cond1 = [dorsalResultsDatabase.DataType]=="1Dg11" & [dorsalResultsDatabase.nc]==12;
cond1 = contains([dorsalResultsDatabase.DataType],'1Dg11_2xDl_FFF') & [dorsalResultsDatabase.nc]==12;
frac_1 = dorsalResultsDatabase.meanFracFluoEmbryo(cond1);
sefrac_1 = dorsalResultsDatabase.seFracFluoEmbryo(cond1);
max_1 = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond1);
semax_1 = dorsalResultsDatabase.seAllMaxFluoEmbryo(cond1);
accum_1 = dorsalResultsDatabase.meanallmrnasEmbryo(cond1);
seaccum_1 = dorsalResultsDatabase.seallmrnasEmbryo(cond1);
turnOn_1 = dorsalResultsDatabase.meanTurnOnsEmbryo(cond1);
seTurnON_1 = dorsalResultsDatabase.seTurnOnsEmbryo(cond1);


convFactor = 1.2371; %with optogenetics settings Venus is this much brighter than in normal settings

cond2 = [dorsalResultsDatabase.DataType]=="1Dg11_noExport" & [dorsalResultsDatabase.nc]==12;
dlfluobins_2 = dorsalResultsDatabase.dorsalFluoBins(cond2);
dlfluobins_2 = dlfluobins_2 ./ convFactor;
frac_2 = dorsalResultsDatabase.meanFracFluoEmbryo(cond2);
sefrac_2 = dorsalResultsDatabase.seFracFluoEmbryo(cond2);
max_2 = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond2);
semax_2 = dorsalResultsDatabase.seAllMaxFluoEmbryo(cond2);
accum_2 = dorsalResultsDatabase.meanallmrnasEmbryo(cond2);
seaccum_2 = dorsalResultsDatabase.seallmrnasEmbryo(cond2);
turnOn_2 = dorsalResultsDatabase.meanTurnOnsEmbryo(cond2);
seTurnON_2 = dorsalResultsDatabase.seTurnOnsEmbryo(cond2);
close all

%% plot results
frac_1(isnan(frac_1)) = 0;
max_1(isnan(max_1))=0;
accum_1(isnan(accum_1))=0;
frac_2(isnan(frac_2)) = 0;
max_2(isnan(max_2))=0;
accum_2(isnan(accum_2))=0;

%results = dorsalResults{1,1}

figure
hold on
errorbar(dlfluobins_1+62.5,frac_1,sefrac_1,'ro-','CapSize',0,'MarkerFaceColor','r')
errorbar(dlfluobins_2+62.5,frac_2,sefrac_2,'bo-','CapSize',0,'MarkerFaceColor','b')
ylabel('fraction')
ylim([0 1])
xlim([0 3250])
hold off
legend('2X and 1X combined','opto Dorsal')
title('fraction active')


figure
hold on
errorbar(dlfluobins_1+62.5,max_1,semax_1,'ro-','CapSize',0,'MarkerFaceColor','r')
%errorbar(dlfluobins+10,maxFFF,semaxFFF,'bo-','CapSize',0,'MarkerFaceColor','b')
errorbar(dlfluobins_2+62.5,max_2,semax_2,'bo-','CapSize',0,'MarkerFaceColor','b')
ylabel('max fluo')
ylim([0 600])
xlim([0 3250])
hold off
legend('2X and 1X combined','opto Dorsal')
title('maximum fluorescence')


figure
hold on
errorbar(dlfluobins_1+62.5,accum_1,seaccum_1,'ro-','CapSize',0,'MarkerFaceColor','r')
errorbar(dlfluobins_1+62.5,accum_2,seaccum_2,'bo-','CapSize',0,'MarkerFaceColor','b')
ylabel('accumulated mRNA')
ylim([0 1200])
xlim([0 3250])
hold off
legend('2X and 1X combined','opto Dorsal')
title('accumulated fluorescence')


figure
hold on
errorbar(dlfluobins_1+62.5,turnOn_1,seTurnON_1,'ro-','CapSize',0,'MarkerFaceColor','r')
errorbar(dlfluobins_1+62.5,turnOn_2,seTurnON_2,'bo-','CapSize',0,'MarkerFaceColor','b')
ylabel('turn on time')
ylim([0 10])
xlim([0 3250])
hold off
legend('2X and 1X combined','opto Dorsal')
title('turn on time')



end
