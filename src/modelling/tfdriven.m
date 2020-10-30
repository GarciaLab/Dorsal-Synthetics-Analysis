

%% setup parameters and stuff
close all;

d = logspace(1, 4); % range of dorsal concentrations in AUs
t = linspace(0, 10)'; % time in mins
kd = 500; % Kd of Dorsal binding in AUs
cs = logspace(-1, 2, 10);  % rate increase per unit dorsal
R = 1; % rate of transcription when the dorsal site is 100% occupied
nSteps = 5; % number of irreversible steps
t_nc13 = 10; % duration of the transcriptional window in minutes

% [D, T] = meshgrid(d,t); %?
% 
% figure;
% t_surf = tiledlayout('flow');

f2 = figure; ax2 = axes(f2);
f3 = figure; ax3 = axes(f3);
f4 = figure; ax4 = axes(f4);
f5 = figure; ax5 = axes(f5);
f7 = figure; ax7 = axes(f7);
f8 = figure; ax8 = axes(f8);

% initialize arrays containing the model outputs
dmrnadt = [];
paccessible  = [];
n = 0; %for indexing
fon = [];
ton2 = [];
cdfg = [];
temp8 = [];
cdfg2 = [];
ton3 = [];
ton4 = [];

% loop over values of c
for c = cs
    
    n = n + 1; %for indexing
    
%     %five irreversible steps, SA: I think this is a wrong equation
%     temp = (-1/24).*d.*exp(-c.*d.*t).*(d+kd).^(-1).*R.*(24+(-24).* ...
%       exp(c.*d.*t)+c.*d.*t.*(24+c.*d.*t.*(12+c.*d.*t.*(4+c.*d.*t))) ...
%       );
%      dmrnadt(n, :, :) = temp;
%      mrna(n,:) = trapz(t, temp, 1);
%  
%  temp2 = -((d R (-24 t + 
%     E^(-c d t) (-(120/(c d)) - 96 t - 36 c d t^2 - 8 c^2 d^2 t^3 - 
%        c^3 d^3 t^4)))/(24 (d + kd)))
%   

%     % probability of being in the last state as a function of [Dl] and time
%     % for a given c
%     temp3 = (1/24)*exp(-c.*d.*t).* (-24 + 24*exp(c.*d.*t) - 24.*c.*d.*t - 12.*(c.*d.*t).^2 - 4*(c.*d.*t).^3 - (c.*d.*t).^4);
%     paccessible(n, :, :) = temp3;
%     
%     %odds of reaching the end of the cycle without turning on
%     fon(n, :) = squeeze(paccessible(n, t_nc13, :)); 

    temp6 = (1/24).*d.*exp(1).^((-1).*c.*d.*(d+kd).^(-1).*t).*(d+kd).^(-1).* ...
    R.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+c.*d.*(d+kd).^(-4).* ...
    t.*((-24).*(d+kd).^3+(-12).*c.*d.*(d+kd).^2.*t+(-4).*c.^2.*d.^2.*( ...
    d+kd).*t.^2+(-1).*c.^3.*d.^3.*t.^3));

    dmrnadt(n, :, :) = temp6;
    mrna2(n,:) = trapz(t, temp6, 1);

    % probability of being in the last state as a function of [Dl] and time
    % for a given c
    temp7 = (1/24).*exp(1).^((-1).*c.*d.*(d+kd).^(-1).*t).*(d+kd).^(-4).*(24.* ...
      ((-1)+exp(1).^(c.*d.*(d+kd).^(-1).*t)).*kd.^4+24.*d.*kd.^3.*((-4)+ ...
      4.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t)+12.*d.^2.*kd.^2.*(( ...
      -12)+12.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t.*(6+c.*t))+4.* ...
      d.^3.*kd.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t.*( ...
      18+c.*t.*(6+c.*t)))+d.^4.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).* ...
      t)+(-1).*c.*t.*(24+c.*t.*(12+c.*t.*(4+c.*t)))));

    paccessible2(n, :, :) = temp7;   
    %odds of reaching the end of the cycle without turning on
    fon2(n, :) = squeeze(paccessible2(n, t_nc13, :));

    temp10 = 5.*c.^(-1).*(1+d.^(-1).*kd+2500.*c.^5.*d.^4.*(3.*d.^4+30.*c.*d.^4+ ...
      150.*c.^2.*d.^4+500.*c.^3.*d.^4+1250.*c.^4.*d.^4+12.*d.^3.*kd+90.* ...
      c.*d.^3.*kd+300.*c.^2.*d.^3.*kd+500.*c.^3.*d.^3.*kd+18.*d.^2.* ...
      kd.^2+90.*c.*d.^2.*kd.^2+150.*c.^2.*d.^2.*kd.^2+12.*d.*kd.^3+30.* ...
      c.*d.*kd.^3+3.*kd.^4+(-3).*exp(1).^(10.*c.*d.*(d+kd).^(-1)).*(d+ ...
      kd).^4).^(-1));

     ton4(n, :) = temp10;

     
    % plot results
    
    % 3D plots
%     nexttile(t_surf)
%     surf(D, T, squeeze(dmrnadt(n, :, :)));
%     xlim([0, 1E4]);
%     title(num2str(c));
    
%     % fraction on as a function of [Dl]
%     plot(ax2, d, fon(n, :))
%     xlabel(ax2,'Dorsal concentration (AU)')
%     ylabel(ax2,'fraction of active nuclei')
%     hold(ax2, 'on')

    % accumulated mRNA
    plot(ax4, d, mrna2(n, :))
    xlabel(ax4,'Dorsal concentration (AU)')
    ylabel(ax4,'acumulated mRNA (AU)')
    hold(ax4, 'on')
    
    % fraction active
    plot(ax5, d, fon2(n, :))
    xlabel(ax5,'Dorsal concentration (AU)')
    ylabel(ax5,'fraction of active nuclei')
    hold(ax5, 'on')
    
    % time on
    plot(ax8, d, ton4(n, :))
    xlabel(ax8,'Dorsal concentration (AU)')
    ylabel(ax8,'time on')
    hold(ax8, 'on')

end

%%


% 
% title(ax2, 'predicted fraction active')
% xlim(ax2, [0, 3500]);
% ylim(ax2,[0, 1]);
% leg2 = legend(ax2, num2str(round(cs', 2, 'significant')));
% title(leg2, 'c')
% xlabel(ax2,'[Dorsal] (au)')
% ylabel(ax2,'fraction active')


title(ax5, 'predicted fraction active')
xlim(ax5, [0, 3500]);
ylim(ax5,[0, 1]);
leg5 = legend(ax5, num2str(round(cs', 2, 'significant')));
title(leg5, 'c')
xlabel(ax5,'[Dorsal] (au)')
ylabel(ax5,'fraction active')


% title(ax3, 'predicted accumulated mRNA')
% xlim(ax3, [0, 3500]);
% leg3 = legend(ax3, num2str(round(cs', 2, 'significant')));
% title(leg3, 'c')
% xlabel(ax3,'[Dorsal] (au)')
% ylabel(ax3,'normalized accumulated mRNA')


title(ax4, 'predicted accumulated mRNA')
xlim(ax4, [0, 3500]);
leg4 = legend(ax4, num2str(round(cs', 2, 'significant')));
title(leg4, 'c')
xlabel(ax4,'[Dorsal] (au)')
ylabel(ax4,'normalized accumulated mRNA')


c = 5;
kd = logspace(0, 4, 100);
d = 1000;
n = 0;
temp10 = [];
ton4 = [];
for k = kd
    n = n + 1;
   temp10 = 5.*c.^(-1).*(1+d.^(-1).*kd+2500.*c.^5.*d.^4.*(3.*d.^4+30.*c.*d.^4+ ...
  150.*c.^2.*d.^4+500.*c.^3.*d.^4+1250.*c.^4.*d.^4+12.*d.^3.*kd+90.* ...
  c.*d.^3.*kd+300.*c.^2.*d.^3.*kd+500.*c.^3.*d.^3.*kd+18.*d.^2.* ...
  kd.^2+90.*c.*d.^2.*kd.^2+150.*c.^2.*d.^2.*kd.^2+12.*d.*kd.^3+30.* ...
  c.*d.*kd.^3+3.*kd.^4+(-3).*exp(1).^(10.*c.*d.*(d+kd).^(-1)).*(d+ ...
  kd).^4).^(-1));

 ton4(n, :) = temp10;
     
plot(ax7, kd, ton4(n,:));
hold(ax7, 'on');

end

title(ax7, 'predicted \langle T_{on} \rangleT_{on} (min)')
xlim(ax7, [0, 3500]);
% leg7 = legend(ax7, num2str(round(cs', 2, 'significant')));
% title(leg7, 'c')
xlabel(ax7,'KD (au)')
ylabel(ax7,'mean turn on time (min)')


title(ax8, 'predicted T_{on} (min)')
xlim(ax8, [0, 3500]);
leg8 = legend(ax8, num2str(round(cs', 2, 'significant')));
title(leg8, 'c')
xlabel(ax8,'[Dl] (au)')
ylabel(ax8,'mean turn on time (min)')

%% Compare with data

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])%, 'dorsalResultsDatabase')
AllStruct = combinedCompiledProjects_allEnhancers([combinedCompiledProjects_allEnhancers.cycle]==12);
b = load('S:\Simon\Dropbox\DorsalSyntheticsDropbox\dorsalResultsDatabase.mat');

Datasets = {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3','1DgVW'};
PatserScores = [6.23,5.81,5.39,5.13,4.8,4.73,4.29];
AllTurnOnData = struct('Embryos',[]);

%gather the data
for e = 1:length(Datasets)
    EnhancerStruct = AllStruct(contains({AllStruct.dataSet}, Datasets{e}));
    EnhancerPrefixes = unique({EnhancerStruct.prefix});   
    EnhancerEmbryos = struct('embryoTurnOnTimes',[]);   
   for p = 1:length(EnhancerPrefixes)
        PrefixTurnOnTimes = nan(19,20); % array of dorsal bins x nuclei containing turn on times
        PrefixStruct = EnhancerStruct(strcmpi({EnhancerStruct.prefix},EnhancerPrefixes{p}) &...
       ~isnan([EnhancerStruct.dorsalFluoBin]));
        onNucleiIdx = find(arrayfun(@(PrefixStruct)~isempty(PrefixStruct.particleTimeOn),PrefixStruct));
        PrefixOnNucleiStruct = PrefixStruct(onNucleiIdx);
        for n = 1:length(PrefixOnNucleiStruct)
            PrefixTurnOnTimes(PrefixOnNucleiStruct(n).dorsalFluoBin,n) = PrefixOnNucleiStruct(n).particleTimeOn;
        end
        EnhancerEmbryos(p).embryoTurnOnTimes = PrefixTurnOnTimes;
   end   
   AllTurnOnData(e).Embryos= EnhancerEmbryos;   
end

% plot stuff
figure
hold on
for enh = 1:length(Datasets)
    
    Embryos = AllTurnOnData(enh).Embryos;
    AllsingleEmbryoMeans = [];
    for emb = 1:length(Embryos)
        singleEmbryoTurnOns = Embryos(emb).embryoTurnOnTimes;
        singleEmbryoMean = nanmean(singleEmbryoTurnOns,2);
        AllsingleEmbryoMeans(emb,:) = singleEmbryoMean;
    end
    NEmbryos = sum(~isnan(AllsingleEmbryoMeans));
    mean = nanmean(AllsingleEmbryoMeans);
    SEM = nanstd(AllsingleEmbryoMeans,1)./NEmbryos;
    mean(mean>10) = nan;
    errorbar([0:250:4500],mean,SEM,'LineWidth',2,'CapSize',0)
    meanTurnOnTimePerKd(enh) = nanmean(mean);
    SEMturnOnTimePerKd(enh) = nanstd(mean)./sum(~isnan(mean));
end
hold off
legend(Datasets)
xlabel('Dorsal concentration')
ylabel('mean turn on time')
ylim([0 10])

figure
errorbar(PatserScores,meanTurnOnTimePerKd,SEMturnOnTimePerKd,'ko','MarkerFaceColor','k','CapSize',0)
%ylim([0 10])
xlabel('patser score')
ylabel('mean turn on time across [Dl]')
set(gca, 'XDir','reverse')


%     
%     
%     AllEmbryosCombined = []
%     for emb = 1:length(Embryos)
%         AllEmbryosCombined = [AllEmbryosCombined,Embryos.embryoTurnOnTimes];
%     end
%     
%     Nnuclei = sum(~isnan(AllEmbryosCombined),2);
%     meanTurnOn = nanmean(AllEmbryosCombined,2);
%     SETurnOn = nanstd(AllEmbryosCombined,2)./Nnuclei;
%     errorbar(meanTurnOn,SETurnOn)
%     ylim([0 10])

        
        
        
   
   












