try
    %this has a function that overloads "mad" used in statistics toolbox
    rmpath('.\lib\dependencies\mcmcstat-master\')
    rmpath('P:\Armando\Dorsal-Synthetics-Analysis\lib\mcmcstat-master')
end

prefix1s ={ '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_1',...
  '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_2',...
  '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_3',...
  '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_5',...
  '2020-11-03-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_6',...
  '2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_8',...
  '2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_9',...
  '2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_10',...
  '2020-11-04-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_11',...
  '2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_18',...
  '2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_19',...
  '2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_20',...
  '2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_21',...
  '2020-12-06-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_22',...
  '2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_26',...
  '2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_27',...
  '2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_28',...
  '2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_29',...
  '2020-12-07-2xIntB2-1Dg456_parB2GFP-2xnosMCPmCh_30' };

N = length(prefix1s);

max_intensities_non = cell(N, 1);
max_intensities_transcription = cell(N, 1);
mean_intensities_non = cell(N, 1);
mean_intensities_transcription = cell(N, 1);

for k = 1:length(prefix1s)
    
    prefix1 = prefix1s{k};
    prefix2 = strcat(prefix1, 'ch2filt');
    
    [max_intensities_non{k},max_intensities_transcription{k},mean_intensities_non{k}, ...
        mean_intensities_transcription{k}] = MeanMax2Channel(prefix1,prefix2);
    
end

New_ParB_1Dg_max_trans = [max_intensities_transcription{:}];
New_ParB_1Dg_max_trans_50 = New_ParB_1Dg_max_trans + 50;

New_ParB_1Dg_max_non = [max_intensities_non{:}];
New_ParB_1Dg_max_non_50 = New_ParB_1Dg_max_non + 50;
