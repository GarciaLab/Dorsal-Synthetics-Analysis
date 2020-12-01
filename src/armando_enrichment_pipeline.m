
dataType = '1Dg11_2xDl';
dataType = '1DgS2_2xDl';

addpath(genpath('S:\Armando\nick_enrichment_fork\tf_enrichment_pipeline'));
% addpath(genpath('S:\Armando\tf_enrichment_pipeline'));

% main01_compile_traces(dataType,'S:\Armando\Dropbox\DorsalSyntheticsDropbox\','firstNC', 12)
main01_compile_traces(dataType,'firstNC', 12)

main02_sample_local_protein(dataType,'ignoreQC', true, 'max_nucleus_radius_um', 6,'segmentNuclei', 1, 'NumWorkers', 1);

main04_make_exploratory_figs(dataType,'S:\Armando\Dropbox\',...
   'dlv', 'mcpmch')
