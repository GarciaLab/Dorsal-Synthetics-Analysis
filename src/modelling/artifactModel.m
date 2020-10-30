v = .67;
L  = [.05, .1, .25, .9];
Ls = .25;
kd = 810;
kds = [250, 500, 1000, 5000];
R = 1;
d = [0.01:1:4000];
n = 1;

dmRNAdt = @(d, kd, n, R) (1+(d.*kd.^(-1)).^n).^(-1).*(d.*kd.^(-1)).^n.*R;

dmRNAdtDist = @(d, kd, n, v, R, L) (0.1E1+(d.*kd.^(-1)).^n).^(-1).*(d.*kd.^(-1)).^n.*R.*( ...
  0.1E1+erf(0.353553E0.*v.^(-1).*(v.^2+(-0.2E1).*log(L)+0.2E1.*log( ...
  R+(-0.1E1).*(0.1E1+(d.*kd.^(-1)).^n).^(-1).*R)))).*erfc( ...
  0.353553E0.*v.^(-1).*(v.^2+0.2E1.*log(L)+(-0.2E1).*log(R+(-0.1E1) ...
  .*(0.1E1+(d.*kd.^(-1)).^n).^(-1).*R))).^(-1);

dmRNAdtObs = @(d, kd, n, v, R, L) (0.1E1+(d.*kd.^(-1)).^n).^(-1).*(d.*kd.^(-1)).^n.*R.*( ...
  0.1E1+erf(0.353553E0.*v.^(-1).*(v.^2+(-0.2E1).*log(L)+0.2E1.*log( ...
  R+(-0.1E1).*(0.1E1+(d.*kd.^(-1)).^n).^(-1).*R)))).*erfc( ...
  0.353553E0.*v.^(-1).*(v.^2+0.2E1.*log(L)+(-0.2E1).*log(R+(-0.1E1) ...
 .*(0.1E1+(d.*kd.^(-1)).^n).^(-1).*R))).^(-1);

fractionActive = @(d, kd, n, v, R, L)  (1/2).*erfc((-1).*2.^(-1/2).*v.^(-1).*((-0.5E0).*v.^2+( ...
  -1).*log(L)+log((1+(d.*kd.^(-1)).^n).^(-1).*(d.*kd.^(-1)).^n.*R)));

close all force

figure; tiledlayout('flow')
for k = 1:length(Ls)
    nexttile;
    plot(d, dmRNAdt(d, kd, n, R), d, dmRNAdtObs(d, kd, n, v, R, Ls(k)), d, fractionActive(d, kd, n, v, R, Ls(k)), 'LineWidth', 2 );
    xlabel('[Dorsal] (au)')
    legend({'dmRNA/dt true', 'dmRNA/dt observed', 'fraction active'});
    title("L = " + Ls(k) );
end


figure; tiledlayout('flow')
for k = 1:length(kds)
    nexttile;
    plot(d, dmRNAdt(d, kds(k), n, R), d, dmRNAdtObs(d, kds(k), n, v, R, L), d, fractionActive(d, kds(k), n, v, R, L), 'LineWidth', 2 );
    xlabel('[Dorsal] (au)')
    legend({'dmRNA/dt true', 'dmRNA/dt observed', 'fraction active'});
    title("KD = " + kds(k) );
end