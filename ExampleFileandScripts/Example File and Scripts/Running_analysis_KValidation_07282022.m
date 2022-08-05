% CONSTRUCT DATA OBJECT 

% LOAD DATA
load('20200422_bl21lb21_example_data.mat')

% SET RELEVANT VARIABLES
% Note- for your data, you make have different info for neurons and
% pressure data
standard_fileInfo = '20200422_bl21lb21_example_data';

% Set behav sample frequency
bFs = 30000;

% Set neural sample frequency
nFs = 30000;

behav_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));

neural_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));

verbose_level = 5;

% MAKE PRESSURE DATA OBJECT

% Construct pressure behavioral object
b = mi_data_pressure(behav_fileInfo, 'verbose', verbose_level);

% Add cycle times to pressure object
b.add_cycleTimes(final_cycle_times', strcat(standard_fileInfo, datestr(date, 'mmddyyyy')), bFs);

% MAKE NEURAL DATA OBJECT

% Set up info strings for each unit
str_UnitG = strcat(standard_fileInfo, '/unitG');
str_UnitD = strcat(standard_fileInfo, '/unitD');

% Construct neural data object
d = mi_data_neural(neural_fileInfo, 'verbose', verbose_level);

% Add spikes to neural data object
d.add_spikes(unitG, str_UnitG, nFs, 'unitG');
d.add_spikes(unitD, str_UnitD, nFs, 'unitD');

% 20220727 LHT: 
fileName = strcat(standard_fileInfo, datestr(date, 'mmddyyyy'), 'dataObjs.mat');

% Save Data Objects
save(fileName, 'd', 'b')

%% 
% CONSTRUCT ANALYSIS OBJECT

% Load Data objects
load('20200422_bl21lb21_example_data07272022dataObjs.mat')

% Set analysis variables of interest
unit1_name = 'unitD';
unit2_name = 'unitG';

% NOTE- SAVE A COPY OF THE SCRIPT WITH A NEW DATA AND UNIT NAMES AFTER
% REVISING THE VARIABLES ABOVE

verbose_level = 4;

% Construct CC analysis object and run estimates
% Construct mi_analysis object
a_cc = calc_count_count(d, b, {unit1_name , unit2_name}, 'verbose', verbose_level);

a_cc.buildMIs();

% For K Selection Training, set rng seed:
a_cc.arrMIcore{1,1}.set_seed = true;

a_cc.calcMIs();

a_cc.getMIs();

save(strcat('20200427_bl21lb21_', datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_analysis_c1c2.mat'), 'a_cc')

% Construct mi_analysis object
a_tc = calc_timing_count(d, b, {unit1_name , unit2_name}, 'verbose', verbose_level, 'discard_omittedData', true);

a_tc.buildMIs();

% For K selection training, set rng seed:
for i = 1:size(a_tc.arrMIcore,1)
    a_tc.arrMIcore{i,1}.set_seed = true;
end

a_tc.calcMIs();

a_tc.getMIs();

save(strcat('20200427_bl21lb21_', datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_analysis_t1c2.mat')', 'a_tc')

% Construct TC analysis object and run estimates

% Construct mi_analysis object
a_tc = calc_timing_count(d, b, {unit2_name , unit1_name}, 'verbose', verbose_level, 'discard_omittedData', true);

a_tc.buildMIs();

% For K selection training, set rng seed:
for i = 1:size(a_tc.arrMIcore,1)
    a_tc.arrMIcore{i,1}.set_seed = true;
end

a_tc.calcMIs();

a_tc.getMIs();

save(strcat('20200427_bl21lb21_', datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_analysis_t2c1.mat'), 'a_tc')

% Construct TT analysis object and run estimates

% Construct mi_analysis object
a_tt = calc_timing_timing(d, b, {unit1_name , unit2_name}, 'verbose', verbose_level, 'discard_omittedData', true);

a_tt.buildMIs();

% For K selection training, set rng seed:
for i = 1:size(a_tt.arrMIcore,1)
    a_tt.arrMIcore{i,1}.set_seed = true;
end

a_tt.calcMIs();

a_tt.getMIs();

save(strcat('20200427_bl21lb21_', datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_analysis_t1t2.mat'), 'a_tt')

% Construct TT analysis object and run estimates

% Construct mi_analysis object
test = calc_timing_timing(d, b, {unit1_name , unit2_name}, 'verbose', verbose_level, 'discard_omittedData', true);

test.buildMIs();

% For K selection training, set rng seed:
for i = 1:size(test.arrMIcore,1)
    test.arrMIcore{i,1}.set_seed = true;
end

test.calcMIs();

test.getMIs();