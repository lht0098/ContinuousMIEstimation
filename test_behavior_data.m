% setup behavior object

standard_fileInfo = '20200422_bl21lb21_example_data';

behav_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));

load('/Users/brycechung/Documents/GitHub/ContinuousMIEstimation/ExampleFileandScripts/Example File and Scripts/20200422_bl21lb21_example_data.mat', 'final_cycle_times');


b = mi_data_pressure(behav_fileInfo, 'verbose', 4);

b.set_data_files('ext', '.rhd');

b.set_calibration('slope', 15.4467, 'first_file_mean', 1.5336, 'pressure_offset', 1.5336);

b.build_behavior();

% b.add_cycleTimes(final_cycle_times', strcat(standard_fileInfo, datestr(date, 'mmddyyyy')), 30000);

b.process_behavior();


%% ===== ===== ===== ===== =====
% Test behavioral metrics
%  ===== ===== ===== ===== =====

%% Time | raw
start = 1;
dur= 200;
samps = 100;

x = b.get_behavior('time', 'raw', start, dur, samps);
start_samples = ceil(start*30000/1000.);
dur_samples = ceil(dur*30000/1000.);
stop_samples = start_samples + dur_samples;

ix = 1;
tmp_dat = b.data.wav_pressure.data{ix}(start_samples:stop_samples);
newSamples = round(linspace(1,length(tmp_dat),samps));

figure();
plot(x');
title('Resampled Waveforms - Time: Raw')

figure();
plot(b.data.wav_pressure.data{1}, 'k-', 'LineWidth', 2);
hold on;
scatter(newSamples+start_samples, x(1,:), 'ro');
title('Resampling verification - Time Raw');


%% Time | pca

start = 1;
dur= 300;
samps = 100;

x = b.get_behavior('time', 'raw', start, dur, samps);


hf_wavs = figure();
% hp_wavs = plot(x', 'tag', 'rollover');
hp_wavs = plot(x', 'tag', 'rollover', 'LineWidth', 1);
title('Resampled Waveforms - Time: PCA');


[~,score,~] = pca(x);
y = b.get_behavior('time', 'pca', start, dur, samps);

hf_pcas = figure();
scatter(score(:,1), score(:,2), 'k.');
hold on;
hp_pcas = scatter(y(:,1), y(:,2), 'ro');
title({'PCA Comparison - Time: PCA' 'k. - manual | ro - library'});


%% Time | residual

start = 1;
dur= 300;
samps = 100;

x = b.get_behavior('time', 'raw', start, dur, samps);

figure();
plot(x');
hold on;
plot(mean(x), 'k-', 'LineWidth', 5);
title('Resampled Waveforms - Time: Residual');

y = b.get_behavior('time', 'residual', start, dur, samps);

figure();
plot(y');
hold on;
plot(mean(y), 'r-', 'LineWidth', 5);
title('Resampled Residuals - Time: Residual');


%% Phase | raw
start = 0;
dur= 2*pi;
samps = 100;

x = b.get_behavior('phase', 'raw', start, dur, samps);
start_samples = ceil(start*30000/1000.);
dur_samples = ceil(dur*30000/1000.);
stop_samples = start_samples + dur_samples;

figure();
plot(x');
title('Resampled Waveforms - Phase Raw');


%% Phase | pca
start = 0;
dur= 2*pi;
samps = 100;

x = b.get_behavior('phase', 'raw', start, dur, samps);

figure();
plot(x');
title('Resampled Waveforms - Phase: PCA')

[~,score,~] = pca(x);
y = b.get_behavior('phase', 'pca', start, dur, samps);

hf_pcas = figure();
scatter(score(:,1), score(:,2), 'k.');
hold on;
hp_pcas = scatter(y(:,1), y(:,2), 'ro');
title({'PCA Comparison - Phase: PCA' 'k. - manual | ro - library'});

%% Phase | residual

start = 0;
dur= 2*pi;
samps = 100;

x = b.get_behavior('phase', 'raw', start, dur, samps);

figure();
plot(x');
hold on;
plot(mean(x), 'k-', 'LineWidth', 5);
title('Resampled Waveforms - Phase: Residual');

y = b.get_behavior('phase', 'residual', start, dur, samps);

figure();
plot(y');
hold on;
plot(mean(y), 'r-', 'LineWidth', 5);
title('Resampled Residuals - Phase: Residual');
