%% Load the data object and set the units to analyze 

clear all

% Load the data object
load('20200427_bl21lb21_06032020data.mat')

% Set analysis variables of interest
unit1_name = 'unitB';
unit2_name = 'unitG';  

verbose_level = 4;

%% Set the parameters for the noise

% Set the mean and standard deivation parameter for the noise
mu = 0;
sigma = 1.5; % CHANGE

% Change the mu and sigma of the gaussian noise to a string
str_mu = num2str(mu);
str_sigma = num2str(sigma);

% Ensure that there are no dots in the name of the file
str_mu(strfind(str_mu, '.')) = [];
str_sigma(strfind(str_sigma, '.')) = [];

save_str = ['20200427_bl21lb21_06032020data' '_' unit2_name '_mu' str_mu '_sig' str_sigma '_'];

%% Construct the analysis object with the subgroup jitter

% Construct mi_analysis object 

% t1c2
% a_tc = calc_timing_count(d, b, {unit1_name , unit2_name}, 'verbose', verbose_level);

% t2c1
a_tc = calc_timing_count(d, b, {unit2_name , unit1_name}, 'verbose', verbose_level);

a_tc.buildMIs();

% Determine which unit is x and y (audit)
units = a_tc.varNames;
unitX_time = units(1) % spike time - x
unitY_count = units(2) % spike count - y

% Set the unit of interest 
unit_of_interest = 'x'; % CHANGE 

% Initate an array to hold: subgroup number, if subgroup has spikes, number
% of rows of subgroup, number of columns of subgroup (audit)
subgroup_info = NaN(length(a_tc.arrMIcore), 4);

for i = 1:length(a_tc.arrMIcore)

    % Get the subgroup number
    subgroup_info(i, 1) = i;
    % Determine if there are spikes in the subgroup  
    subgroup_info(i, 2) = isempty(a_tc.arrMIcore{i,1}.(unit_of_interest));
    % Get the number of rows and the number of columns 
    subgroup_info(i, 3:4) = size(a_tc.arrMIcore{i,1}.(unit_of_interest));
    % If the unit in the subgroup does have spike times
    
    if ~isempty(a_tc.arrMIcore{i,1}.(unit_of_interest));
        % Create a guassian noise array of the same dimensions as the
        % unit's spike times
        noise = normrnd(mu, sigma, subgroup_info(i,3), subgroup_info(i,4));
        % Add the noise to the spike time
        noisySpikeTime = noise + a_tc.arrMIcore{i, 1}.(unit_of_interest);     
        % Replace the spike time with the jittered spike time 
        a_tc.arrMIcore{i, 1}.(unit_of_interest) = noisySpikeTime;
    end 
    
end 

a_tc.calcMIs();

a_tc.getMIs();

% t1c2
% save(strcat(save_str, datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_subgroupJitter_analysis_t1c2.mat')', 'a_tc')

% t2c1
save(strcat(save_str, datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_subgroupJitter_analysis_t2c1.mat'), 'a_tc')
