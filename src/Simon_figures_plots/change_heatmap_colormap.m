
%% this is for the states heatmap
NumStates = 7;
OnColor = [.6 .96 .64];
Palette = (cbrewer('seq', 'OrRd',NumStates));
%Palette = flip(cbrewer('seq', 'YlGnBu',NumStates));
Palette(end,:) = OnColor;
colormap(Palette)

%% this is for the states heatmap
NumStates = 6;
OnColor = [.6 .96 .64];
Palette = flip(cbrewer('seq', 'OrRd',NumStates));
%Palette = flip(cbrewer('seq', 'YlGnBu',NumStates));
Palette = [Palette;OnColor];
colormap(Palette)
