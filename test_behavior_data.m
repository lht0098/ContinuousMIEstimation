% setup behavior object

behav_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));

b = mi_data_pressure(behav_fileInfo, 'verbose', 5);

b.set_data_files('ext', '.rhd');

b.set_calibration('slope', 15.4467, 'first_file_mean', 1.5336, 'pressure_offset', 1.5336);

b.build_behavior();