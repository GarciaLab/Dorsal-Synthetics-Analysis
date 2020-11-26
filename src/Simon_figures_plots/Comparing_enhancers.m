
dorsalResultsDatabase = createDorsalResultsDatabase;
dlfluobins = dorsalResultsDatabase.dorsalFluoBins; dlfluobins = dlfluobins(1:19);
cond1 = [dorsalResultsDatabase.DataType]=="TwiPEv5(7)_2xDl" & [dorsalResultsDatabase.nc]==12;
fracTwi = dorsalResultsDatabase.meanFracFluoEmbryo(cond1);
sefracTwi = dorsalResultsDatabase.seFracFluoEmbryo(cond1);
maxTwi = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond1);
semaxTwi = dorsalResultsDatabase.seAllMaxFluoEmbryo(cond1);
accumTwi = dorsalResultsDatabase.meanallmrnasEmbryo(cond1);
seaccumTwi = dorsalResultsDatabase.seallmrnasEmbryo(cond1);

cond2 = [dorsalResultsDatabase.DataType]=="1Dg11_2xDl" & [dorsalResultsDatabase.nc]==12;
frac1Dg = dorsalResultsDatabase.meanFracFluoEmbryo(cond2);
sefrac1Dg = dorsalResultsDatabase.seFracFluoEmbryo(cond2);
max1Dg = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond2);
semax1Dg = dorsalResultsDatabase.seAllMaxFluoEmbryo(cond2);
accum1Dg = dorsalResultsDatabase.meanallmrnasEmbryo(cond2);
seaccum1Dg = dorsalResultsDatabase.seallmrnasEmbryo(cond2);




%results = dorsalResults{1,1}

figure
hold on
errorbar(dlfluobins,fracTwi,sefracTwi,'b-','CapSize',0)
errorbar(dlfluobins,frac1Dg,sefrac1Dg,'r-','CapSize',0)
%errorbar(dlfluobins,results.meanFracFluoEmbryo,results.seFracFluoEmbryo,'r-','CapSize',0)
ylabel('fraction')
ylim([0 1])
hold off
legend('TwiPE','1Dg')


figure
hold on
errorbar(dlfluobins,maxTwi,semaxTwi,'b-','CapSize',0)
errorbar(dlfluobins,max1Dg,semax1Dg,'r-','CapSize',0)
%errorbar(dlfluobins,results.meanAllMaxFluoEmbryo,results.seAllMaxFluoEmbryo,'r-','CapSize',0)
ylabel('max fluo')
ylim([0 400])
hold off
legend('TwiPE','1Dg')


figure
hold on
errorbar(dlfluobins,accumTwi,seaccumTwi,'b-','CapSize',0)
errorbar(dlfluobins,accum1Dg,seaccum1Dg,'r-','CapSize',0)
%errorbar(dlfluobins,results.meanallmrnasEmbryo,results.seallmrnasEmbryo,'r-','CapSize',0)
ylabel('accumulated mRNA')
ylim([0 1200])
hold off
legend('TwiPE','1Dg')
