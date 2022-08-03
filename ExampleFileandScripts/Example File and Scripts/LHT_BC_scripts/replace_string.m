%% Replace the sigma string in the k plot files 

clear all

% Get all the png files in the folder 
files = dir('*png') ; 

% List of sigma strings that need to be changed 
sigma_string = {'sig_05', 'sig_15', 'sig_30', 'sig_60', 'sig_90'};

new_sigma_string = {'sig05', 'sig15', 'sig30', 'sig60', 'sig90'};

% SET THE SIGMA THAT WILL BE CHANGED 
sigma = num2str(cell2mat(sigma_string(1))); 

% SET THE SIGMA THAT WILL BE CHANGED 
new_sigma = num2str(cell2mat(new_sigma_string(1))); 

% Go through all of the files and change only the sigma string
for i = 1:length(files)
    front_string = extractBefore(files(i).name, sigma);
    back_string = extractAfter(files(i).name, sigma);
    newName = [front_string, new_sigma, back_string]; 
    movefile(files(i).name, newName)
end

%% Replace the sigma string in the mat files 

clearvars -except sigma new_sigma

% Get all the png files in the folder 
files = dir('*mat') ; 

% Go through all of the files and change only the sigma string
for i = 1:length(files)
    front_string = extractBefore(files(i).name, sigma);
    back_string = extractAfter(files(i).name, sigma);
    newName = [front_string, new_sigma, back_string]; 
    movefile(files(i).name, newName)
end