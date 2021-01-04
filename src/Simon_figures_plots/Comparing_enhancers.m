function Comparing_enhancers

%% grab results
clear all
dorsalResultsDatabase = createDorsalResultsDatabase;
dlfluobins = dorsalResultsDatabase.dorsalFluoBins; dlfluobins = dlfluobins(1:37);
%cond1 = [dorsalResultsDatabase.DataType]=="1Dg11" & [dorsalResultsDatabase.nc]==12;
cond1 = contains([dorsalResultsDatabase.DataType],'1Dg11_2xDl_FFF') & [dorsalResultsDatabase.nc]==12;
frac2x = dorsalResultsDatabase.meanFracFluoEmbryo(cond1);
sefrac2x = dorsalResultsDatabase.seFracFluoEmbryo(cond1);
max2x = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond1);
semax2x = dorsalResultsDatabase.seAllMaxFluoEmbryo(cond1);
accum2x = dorsalResultsDatabase.meanallmrnasEmbryo(cond1);
seaccum2x = dorsalResultsDatabase.seallmrnasEmbryo(cond1);
turnOn2x = dorsalResultsDatabase.meanTurnOnsEmbryo(cond1);
seTurnON2x = dorsalResultsDatabase.seTurnOnsEmbryo(cond1);

% cond2 = [dorsalResultsDatabase.DataType]=="1Dg11_2xDl" & [dorsalResultsDatabase.nc]==12;
% dlfluobinsfff = dorsalResultsDatabase.dorsalFluoBins(cond2);
% fracFFF = dorsalResultsDatabase.meanFracFluoEmbryo(cond2);
% sefracFFF = dorsalResultsDatabase.seFracFluoEmbryo(cond2);
% maxFFF = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond2);
% semaxFFF = dorsalResultsDatabase.seAllMaxFluoEmbryo(cond2);
% accumFFF = dorsalResultsDatabase.meanallmrnasEmbryo(cond2);
% seaccumFFF = dorsalResultsDatabase.seallmrnasEmbryo(cond2);

close all

%% plot results
frac2x(isnan(frac2x)) = 0;
max2x(isnan(max2x))=0;
accum2x(isnan(accum2x))=0;

%results = dorsalResults{1,1}

figure
hold on
errorbar(dlfluobins+62.5,frac2x,sefrac2x,'ro-','CapSize',0,'MarkerFaceColor','r')
%errorbar(dlfluobins+10,fracFFF,sefracFFF,'bo-','CapSize',0,'MarkerFaceColor','b')
%errorbar(dlfluobins,results.meanFracFluoEmbryo,results.seFracFluoEmbryo,'r-','CapSize',0)
ylabel('fraction')
ylim([0 1])
xlim([0 3250])
hold off
legend('2X and FFF combined','FFF')
title('1Dg')


figure
hold on
errorbar(dlfluobins+62.5,max2x,semax2x,'ro-','CapSize',0,'MarkerFaceColor','r')
%errorbar(dlfluobins+10,maxFFF,semaxFFF,'bo-','CapSize',0,'MarkerFaceColor','b')
%errorbar(dlfluobins,results.meanAllMaxFluoEmbryo,results.seAllMaxFluoEmbryo,'r-','CapSize',0)
ylabel('max fluo')
ylim([0 600])
xlim([0 3250])
hold off
legend('2X and FFF combined','FFF')
title('1Dg')


figure
hold on
errorbar(dlfluobins+62.5,accum2x,seaccum2x,'ro-','CapSize',0,'MarkerFaceColor','r')
%errorbar(dlfluobins+10,accumFFF,seaccumFFF,'bo-','CapSize',0,'MarkerFaceColor','b')
%errorbar(dlfluobins,results.meanallmrnasEmbryo,results.seallmrnasEmbryo,'r-','CapSize',0)
ylabel('accumulated mRNA')
ylim([0 1200])
xlim([0 3250])
hold off
legend('2X and FFF combined','FFF')
title('1Dg')


figure
hold on
errorbar(dlfluobins+62.5,turnOn2x,seTurnON2x,'ro-','CapSize',0,'MarkerFaceColor','r')
%errorbar(dlfluobins+10,accumFFF,seaccumFFF,'bo-','CapSize',0,'MarkerFaceColor','b')
%errorbar(dlfluobins,results.meanallmrnasEmbryo,results.seallmrnasEmbryo,'r-','CapSize',0)
ylabel('turn on time')
ylim([0 10])
xlim([0 3250])
hold off
legend('2X and FFF combined','FFF')
title('1Dg')



end
