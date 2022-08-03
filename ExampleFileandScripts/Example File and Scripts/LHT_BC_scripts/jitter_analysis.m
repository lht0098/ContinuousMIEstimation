% order of neurons: unit D, unit B, unit E, and unit G
% change unit D (corrupted), and then run it with B, E, and G (which are
% not corrupted)
%% Create the corrupted neuronal unit

clear all

% Load the data object that contains the units:
load('20200427_bl21lb21_06032020data.mat')

% Set the unit of interest 
unit = 'unitG'; % CHANGE

% Get the unit's spike train
spikeTrain = getfield(d, 'data', unit, 'data');

% Generate guassian noise array for the same length with mean value of [1e-4,
% 0.005, 0.02] s - at a time:

% Set the mean and standard deivation parameter
mu = 0;
sigma = 1800; % CHANGE

% Create the guassian noise array
noise = normrnd(mu, sigma, size(spikeTrain, 1), size(spikeTrain, 2));

% Add that gaussian nosie to the spike train 
noisy_spikeTrain = spikeTrain + noise;

% Apply a randome permutation to the spike train
% noisy_spikeTrain = spikeTrain(randperm(length(spikeTrain)));

figure()

histogram(noise)

%% Create a new name for the copy of the data object that will contain the 
% corrupted unit: 

% Change the mu and sigma of the gaussian noise to a string
str_mu = num2str(mu);
str_sigma = num2str(sigma);

% Ensure that there are no dots in the name of the file
str_mu(strfind(str_mu, '.')) = [];
str_sigma(strfind(str_sigma, '.')) = [];

% Create the new name for the corrupted data object file
jitter_file = ['20200427_bl21lb21_06032020data' '_' unit '_mu' str_mu '_sig' str_sigma '_' '.mat' ];

% Create a copy of the source data object 
copyfile('20200427_bl21lb21_06032020data.mat', jitter_file);

%% Input the corrupted unit into the copy of the data object

% Load the copy of the data object
load(jitter_file);

% Delete the uncorrupted spike train
d.data.(unit) = rmfield(d.data.(unit), 'data');

% Input the corrupted spike train into the copy of the data object
d.data.(unit).data = noisy_spikeTrain;

% Clear all the variables 
% clearvars -except b d jitter_file noise

% Save the changes made to the file
save(jitter_file, 'b', 'd');
