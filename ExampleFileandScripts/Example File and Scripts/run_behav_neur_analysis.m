%% =====  =====  =====  =====  ===== 
% Run full analysis for: NN, NB, NNB
%  =====  =====  =====  =====  ===== 

% fix to set unix search path
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/brycechung/Documents/GitHub/ContinuousMIEstimation/kraskovStoegbauerGrassberger']);

%%

unit1 = 'B';
unit2 = 'D';

standard_fileInfo = '20200422_bl21lb21_example_data';

behav_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));
neural_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));

nFs = 30000;
bFs = 30000;

verbose_level = 5;

%% MAKE PRESSURE DATA OBJECT

% load('/Users/brycechung/Documents/GitHub/ContinuousMIEstimation/ExampleFileandScripts/Example File and Scripts/20200422_bl21lb21_example_data.mat');
load('/Users/brycechung/Google Drive/My Drive/__Research/__SOBER/__PROJECTS/_Analysis/Mutual Information/Validation/Transfer of Knowledge/Raw Data Files/20200422_bl21lb21-spiketimes-peaks.mat');

% unit1_dat = eval(['unit' unit1]);
% unit2_dat = eval(['unit' unit2]);

unit1_dat = spike_times_adjusted.(unit1);
unit2_dat = spike_times_adjusted.(unit2);


b = mi_data_pressure(behav_fileInfo, 'verbose', 4);

b.set_data_files('ext', '.rhd');

b.set_calibration('slope', 15.4467, 'first_file_mean', 1.5336, 'pressure_offset', 1.5336);

b.build_behavior();

% b.add_cycleTimes(final_cycle_times', strcat(standard_fileInfo, datestr(date, 'mmddyyyy')), 30000);

% Adjust time offsets for sampling error with DataView
% dtvw_sampling_int = 0.0333000011742115;
% samp_n = 27014400;
% 
% dat_ts = [];
% dat_ts(end+1:samp_n+1) = ([1:samp_n]*dtvw_sampling_int/1000.) + 0.0;
% dat_ts(end+1:samp_n+1) = ([1:samp_n]*dtvw_sampling_int/1000.) + 900.4800;
% dat_ts(end+1:samp_n+1) = ([1:samp_n]*dtvw_sampling_int/1000.) + 1800.9600;
% dat_ts(end+1:samp_n+1) = ([1:samp_n]*dtvw_sampling_int/1000.) + 2701.4400;
% dat_ts(end+1:samp_n+1) = ([1:samp_n]*dtvw_sampling_int/1000.) + 3601.9200;
% 
% b.data.dat_ts.data = dat_ts;


b.process_behavior();


%% MAKE NEURAL DATA OBJECT

% Set up info strings for each unit
str_Unit1 = strcat(standard_fileInfo, ['/unit' unit1]);
% str_Unit2 = strcat(standard_fileInfo, ['/unit' unit2]);

% Construct neural data object
d = mi_data_neural(neural_fileInfo, 'verbose', verbose_level);

% Add spikes to neural data object
d.add_spikes(unit1_dat, str_Unit1, nFs, ['unit' unit1]);
% d.add_spikes(unit2_dat, str_Unit1, nFs, ['unit' unit2]);


%% make raster plot to confirm cycles

unitName = unit1;

cycleTimes = b.get_data('cycleTimes');
spikeTimes = d.get_data(['unit' unitName]);

mat_raster_ixs = nan(size(spikeTimes,2),1);
mat_raster_ts = nan(size(spikeTimes,2), 1);

for i=1:size(cycleTimes,1)
    ixs = find(spikeTimes >= cycleTimes(i,1) & spikeTimes <= cycleTimes(i,2));
    mat_raster_ts(ixs) = spikeTimes(ixs) - cycleTimes(i,1);
    mat_raster_ixs(ixs) = i;
end

figure();
scatter(mat_raster_ts, mat_raster_ixs, 'k.');


%%

a_cb = calc_count_behav(d, b, ...
    {['unit' unit1]}, ...
    'b_timeBase', 'phase', ...
    'feature', 'pca', ...
    'start', 0, ...
    'dur', 2*pi, ...
    'nSamp', 10, ...
    'nPC', 3, ...
    'verbose', verbose_level);

a_cb.buildMIs();


a_cb.arrMIcore{1}.x = a_cb.arrMIcore{1}.x';
a_cb.arrMIcore{1}.y = a_cb.arrMIcore{1}.y';


a_cb.calcMIs();

% a_cb.returnMIs();

%%

% make_kSelectionPlots

a_cb.auditKs();

a_cb.returnMIs();