function plotGoodCurves(factive, dt, mfpts, params)

goodMatrixIndices = plotTFDrivenParams(factive, dt, mfpts, 'dim', 2);

%sort to get dorsals arrayed contiguously
goodMatrixIndices = sortrows(goodMatrixIndices,[5 4 3 2]);

%get locations of contiguous dorsal arrays
pos = findArray(goodMatrixIndices(:, 1), length(params.dls));

factive_theory = []; dt_theory = []; onset_theory = [];
figure;
tiledlayout('flow')

for j = 1:length(pos)
    
    inds = pos(j) : pos(j)+length(params.dls)-1;
    gmi = goodMatrixIndices( inds, : );
    
    for k = 1:length(params.dls)
        factive_theory(k) = factive(gmi(k, 1), gmi(k, 2), gmi(k, 3), gmi(k, 4), gmi(k, 5));
        dt_theory(k) = dt(gmi(k, 1), gmi(k, 2), gmi(k, 3), gmi(k, 4), gmi(k, 5));
        onset_theory(k) = mfpts(gmi(k, 1), gmi(k, 2), gmi(k, 3), gmi(k, 4), gmi(k, 5));
    end
    
    if factive_theory(1) < .2 && factive_theory(end) > .6 &&...
           factive_theory(3) < .8
        %only plot curves that span the full factive range
        
        nexttile(1)
        plot(params.dls, factive_theory, 'LineWidth', 2)
        xlim([0, 4000])
        xlabel('[Dorsal] (au)')
        ylabel('fraction of active nuclei')
        ylim([0, 1])
        hold on
        
        nexttile(2)
        plot(params.dls, dt_theory, 'LineWidth', 2)
        xlabel('[Dorsal] (au)')
        ylabel('Change in mean turn on time across large range of affinities (min)')
        xlim([0, 4000])
        hold on
        
        
        nexttile(3)
        plot(params.dls, onset_theory, 'LineWidth', 2)
        hold on
        xlim([0, 4000])
        set(gca, 'XScale', 'log');
        ylim([0, 10])
        ylabel('mean time to turn on (min)')
        xlabel('[Dorsal] (au)')
    end
    
end