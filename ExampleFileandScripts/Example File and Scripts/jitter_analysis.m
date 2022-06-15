% order of neurons: unit D, unit B, unit E, and unit G
% change unit D (corrupted), and then run it with B, E, and G (which are
% not corrupted)
%% Create the corrupted neuronal unit

clear all

% Load the data object that contains the units:
load('20200427_bl21lb21_06032020data.mat')

% Set the unit of interest
unit = 'unitD';

% Get the unit's spike train
spikeTrain = getfield(d, 'data', unit, 'data');

% Generate guassian noise array for the same length with mean value of [1e-4,
% 0.005, 0.02] s - at a time:

% Set the mean and standard deivation parameter
mu = 0;
sigma = 0.02;

% Create the guassian noise array
noise = normrnd(mu, sigma, size(spikeTrain, 1), size(spikeTrain, 2));

% Add that gaussian nosie to the spike train 

noisy_spikeTrain = spikeTrain + noise;

figure()
plot(noise)

%% Create a new name for the copy of the data object that will contain the 
% corrupted unit: 

% Change the mu and sigma of the gaussian noise to a string
str_mu = num2str(mu);
str_sigma = num2str(sigma);

% Ensure that there are no dots in the name of the file
str_mu(strfind(str_mu, '.')) = [];
str_sigma(strfind(str_sigma, '.')) = [];

% Create the new name for the corrupted data object file
corrupt_file = ['20200427_bl21lb21_06032020data' '_' unit '_mu' str_mu '_sig' str_sigma '.mat' ];

% Create a copy of the source data object 
copyfile('20200427_bl21lb21_06032020data.mat', corrupt_file);

%% Input the corrupted unit into the copy of the data object

% Load the copy of the data object
load(corrupt_file);

% Delete the uncorrupted spike train
d.data.(unit) = rmfield(d.data.(unit), 'data');

% Input the corrupted spike train into the copy of the data object
d.data.(unit).data = noisy_spikeTrain;

% Save the changes to the field
save(corrupt_file);

%% Use that jitter spike train to rerun all analyses with the other uncorrupted spike
%trains

% Get the mutual information for all of the analysis objects 

a_cc.returnMIs

a_tc.returnMIs

a_tt.returnMIs

%% Result: series of bar plots that are data unjittered and data not
% jittered: 
