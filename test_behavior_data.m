% setup behavior object

standard_fileInfo = '20200422_bl21lb21_example_data';
behav_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));

load('/Users/brycechung/Documents/GitHub/ContinuousMIEstimation/ExampleFileandScripts/Example File and Scripts/20200422_bl21lb21_example_data.mat', 'final_cycle_times');


b = mi_data_pressure(behav_fileInfo, 'verbose', 5);

b.set_data_files('ext', '.rhd');

b.set_calibration('slope', 15.4467, 'first_file_mean', 1.5336, 'pressure_offset', 1.5336);

b.build_behavior();

b.add_cycleTimes(final_cycle_times', strcat(standard_fileInfo, datestr(date, 'mmddyyyy')), 30000);

b.process_behavior();