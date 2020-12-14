 predicted_df  = readtable('C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\dat\\predicted_df.csv');
 datTempAve = readtable('C:\\Users\\owner\\Documents\\Dorsal-Synthetics-Analysis\\dat\\dat.csv');
 %%
 enhancers = table2array(unique(datTempAve(:, 4)));
 figure; tiledlayout(1, length(enhancers))
 
 
 cmap = colormap(viridis(length(enhancers)));

 for k = 1:size(enhancers, 1)
 nexttile
 
 hold on
 
 fillyy(predicted_df(predicted_df.dsid==k, :).x, predicted_df(predicted_df.dsid==k, :).Y_mean_cred_lower,...
     predicted_df(predicted_df.dsid==k, :).Y_mean_cred_upper,[0.9 0.9 0.9]);
  plot(predicted_df(predicted_df.dsid==k, :).x, predicted_df(predicted_df.dsid==k, :).Y_mean_mean,  'LineWidth', 2, 'Color', cmap(k, :))
 plot(datTempAve(datTempAve.dsid==k, :).x, datTempAve(datTempAve.dsid==k, :).Y, 'k')
set(gca, 'XScale', 'log')
 title(enhancers{k})
%  ylabel('Max trace fluorescence (au)')
%  xlabel('Dorsal concentration (au)')
xlim([0, 4500]); 

 end
 %%