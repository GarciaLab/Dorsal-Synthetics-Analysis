function exploreBasicModel(theta)



%theta = [c,kd,nInactive,nOff,piEntry,piExit,tCycle];
dorsalVals = linspace(0,3800,18);
modelOpts.TimeVariantDorsalValues = [];
modelOpts.nSims = 1;
fraction_onset = BasicModel_masterEq_DorsalTrace_AR(dorsalVals,theta, modelOpts);

yyaxis left
plot(dorsalVals,fraction_onset(:,1))
ylim([0 1])

yyaxis right
plot(dorsalVals,fraction_onset(:,2))
ylim([0 8])
