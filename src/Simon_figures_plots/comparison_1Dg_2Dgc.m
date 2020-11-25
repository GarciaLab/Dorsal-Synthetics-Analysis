cond = [dorsalResultsDatabase.DataType]=="1Dg11_2xDl" & [dorsalResultsDatabase.nc]==12;
frac1dg = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
sefrac1dg = dorsalResultsDatabase.seFracFluoEmbryo(cond);
max1dg = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond);
semax1dg = dorsalResultsDatabase.seAllMaxFluoEmbryo(cond);
accum1dg = dorsalResultsDatabase.meanallmrnasEmbryo(cond);
seaccum1dg = dorsalResultsDatabase.seallmrnasEmbryo(cond);


results = dorsalResults{1,1}

figure
hold on
errorbar(dlfluobins,frac1dg,sefrac1dg,'b-','CapSize',0)
errorbar(dlfluobins,results.meanFracFluoEmbryo,results.seFracFluoEmbryo,'r-','CapSize',0)
ylabel('fraction')
ylim([0 1])
hold off


figure
hold on
errorbar(dlfluobins,max1dg,semax1dg,'b-','CapSize',0)
errorbar(dlfluobins,results.meanAllMaxFluoEmbryo,results.seAllMaxFluoEmbryo,'r-','CapSize',0)
ylabel('max fluo')
ylim([0 400])
hold off

figure
hold on
errorbar(dlfluobins,accum1dg,seaccum1dg,'b-','CapSize',0)
errorbar(dlfluobins,results.meanallmrnasEmbryo,results.seallmrnasEmbryo,'r-','CapSize',0)
ylabel('accumulated mRNA')
ylim([0 1200])
hold off