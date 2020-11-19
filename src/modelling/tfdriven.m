function tfdriven()
close all;

nPlots = 10;
dmax = 5000;
d = logspace(1, log10(dmax)); %aus
t = linspace(0, 10)'; %mins
kd = 500;
kds = logspace(2, 4, nPlots);
cs = logspace(-1, 2, nPlots);
R = 500;
nSteps = 5;
t_nc13 = 10;
cmap = colormap(viridis(nPlots));


%%
%functions
tonfun  = @(d, c, kd) 5.*c.^(-1).*(1+d.^(-1).*kd+2500.*c.^5.*d.^4.*(3.*d.^4+30.*c.*d.^4+ ...
    150.*c.^2.*d.^4+500.*c.^3.*d.^4+1250.*c.^4.*d.^4+12.*d.^3.*kd+90.* ...
    c.*d.^3.*kd+300.*c.^2.*d.^3.*kd+500.*c.^3.*d.^3.*kd+18.*d.^2.* ...
    kd.^2+90.*c.*d.^2.*kd.^2+150.*c.^2.*d.^2.*kd.^2+12.*d.*kd.^3+30.* ...
    c.*d.*kd.^3+3.*kd.^4+(-3).*exp(1).^(10.*c.*d.*(d+kd).^(-1)).*(d+ ...
    kd).^4).^(-1));

[D, T] = meshgrid(d,t);

figure;
t_surf = tiledlayout('flow');


f2 = figure; ax2 = axes(f2);
f3 = figure; ax3 = axes(f3);
f4 = figure; ax4 = axes(f4);
f5 = figure; ax5 = axes(f5);


dmrnadt = [];
paccessible  = [];
n = 0;
fon = [];
ton2 = [];
cdfg = [];
temp8 = [];
cdfg2 = [];
ton3 = [];
ton4 = [];
for c = cs
    n = n + 1;
    %five irreversible steps
temp= (-1/24).*d.*exp(-c.*d.*t).*(d+kd).^(-1).*R.*(24+(-24).* ...
  exp(c.*d.*t)+c.*d.*t.*(24+c.*d.*t.*(12+c.*d.*t.*(4+c.*d.*t))) ...
  );
 dmrnadt(n, :, :) = temp;
%  
%  temp2 = -((d R (-24 t + 
%     E^(-c d t) (-(120/(c d)) - 96 t - 36 c d t^2 - 8 c^2 d^2 t^3 - 
%        c^3 d^3 t^4)))/(24 (d + kd)))
%        
   temp3 = (1/24)*exp(-c.*d.*t).* (-24 + 24*exp(c.*d.*t) - 24.*c.*d.*t - 12.*(c.*d.*t).^2 - 4*(c.*d.*t).^3 - (c.*d.*t).^4);
    
   paccessible(n, :, :) = temp3;
   
   fon(n, :) = squeeze(paccessible(n, t_nc13, :)); %odds of reaching the end of the cycle without turning on
   
 
  
   
  
   
 mrna(n,:) = trapz(t, temp, 1);
 
 temp6 = (1/24).*d.*exp(1).^((-1).*c.*d.*(d+kd).^(-1).*t).*(d+kd).^(-1).* ...
  R.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+c.*d.*(d+kd).^(-4).* ...
  t.*((-24).*(d+kd).^3+(-12).*c.*d.*(d+kd).^2.*t+(-4).*c.^2.*d.^2.*( ...
  d+kd).*t.^2+(-1).*c.^3.*d.^3.*t.^3));
 
 dmrnadt(n, :, :) = temp6;
 mrna2(n,:) = trapz(t, temp6, 1);

 

temp7 = (1/24).*exp(1).^((-1).*c.*d.*(d+kd).^(-1).*t).*(d+kd).^(-4).*(24.* ...
  ((-1)+exp(1).^(c.*d.*(d+kd).^(-1).*t)).*kd.^4+24.*d.*kd.^3.*((-4)+ ...
  4.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t)+12.*d.^2.*kd.^2.*(( ...
  -12)+12.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t.*(6+c.*t))+4.* ...
  d.^3.*kd.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).*t)+(-1).*c.*t.*( ...
  18+c.*t.*(6+c.*t)))+d.^4.*((-24)+24.*exp(1).^(c.*d.*(d+kd).^(-1).* ...
  t)+(-1).*c.*t.*(24+c.*t.*(12+c.*t.*(4+c.*t)))));
 

paccessible2(n, :, :) = temp7;   
fon2(n, :) = squeeze(paccessible2(n, t_nc13, :)); %odds of reaching the end of the cycle without turning on
 
nexttile(t_surf)
surf(D, T, squeeze(dmrnadt(n, :, :)));
xlim([0, 1E4]);
title(num2str(c));



plot(ax2, d, fon(n, :))
hold(ax2, 'on')


plot(ax3, d, mrna(n, :))
hold(ax3, 'on')


plot(ax4, d, mrna2(n, :), 'LineWidth', 2, 'Color', cmap(n, :))
hold(ax4, 'on')


plot(ax5, d, fon2(n, :),  'LineWidth', 2, 'Color', cmap(n, :))
hold(ax5, 'on')
end





title(ax2, 'predicted fraction active')
xlim(ax2, [0, 3500]);
ylim(ax2,[0, 1]);
leg2 = legend(ax2, num2str(round(cs', 2, 'significant')));
title(leg2, 'c')
xlabel(ax2,'[Dorsal] (au)')
ylabel(ax2,'fraction active')




title(ax3, 'predicted accumulated mRNA')
xlim(ax3, [0, 3500]);
leg3 = legend(ax3, num2str(round(cs', 2, 'significant')));
title(leg3, 'c')
xlabel(ax3,'[Dorsal] (au)')
ylabel(ax3,'normalized accumulated mRNA')


title(ax4, 'predicted accumulated mRNA.  pi \alpha occupancy')
xlim(ax4, [10, dmax]);
leg4 = legend(ax4, num2str(round(cs', 2, 'significant')));
title(leg4, 'c')
xlabel(ax4,'[Dorsal] (au)')
ylabel(ax4,'accumulated mRNA (au)')
set(ax4, 'XScale', 'log')


title(ax5, 'predicted fraction active. pi \alpha occupancy')
xlim(ax5, [10, dmax]);
ylim(ax5,[0, 1]);
% leg5 = legend(ax5, num2str(round(cs', 2, 'significant')));
% title(leg5, 'c')
xlabel(ax5,'[Dorsal] (au)')
ylabel(ax5,'fraction active')
set(ax5, 'XScale', 'log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


figure;
c = 10;

ton4 = [];
for k = 1:length(kds)
    ton4(k, :) = tonfun(d, c, kds(k));
 plot(d, ton4(k, :), 'LineWidth', 2, 'Color', cmap(k, :))
hold on;
end


title('predicted T_{on} (min)')
xlim([10, dmax]);
set(gca, 'XScale', 'log')
ylim([0, 10]);
leg = legend(num2str(round(kds', 2, 'significant')));
title(leg, 'K_D')
xlabel('[Dl] (au)')
ylabel('mean turn on time (min)')

end