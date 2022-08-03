%% Sort the file names from lowest to greatest sigma

clear all

cd % get the current directoy
myDir = cd; % set the current dir as the dir of interest
myFiles = dir(fullfile(myDir,'*.mat')); % puts all of the files into a struct
fileNames = repmat({'NaN'}, 1, length(myFiles)); % saves files for audit

% Create an array that will hold the mi and err 
mi_and_err = NaN(2, length(myFiles));

% List of all possible sigma strings
sigs = {'sig0_', 'sig00001_', 'sig0005_', 'sig002_', 'sig05_', 'sig15_', 'sig30_', 'sig60_', 'sig90_', 'sig900_'};  

% Save the file names 
for i = 1:length(myFiles)
    fileNames(1,i) = cellstr(myFiles(i).name); % save the file name for auditing
end 

% Initiate an array to save the sorted file names
sort_fileNames = repmat({'NaN'}, 1, length(fileNames));

% Sort the files from the lowest to the greatest sigma
for i = 1:length(sort_fileNames)
    % check if the index of the file name and the simga matches
    if contains(fileNames(1, i), sigs(i))
        sort_fileNames(1, i) = fileNames(1, i);
    else 
        % save the file name if they do not match
        file_of_interest = fileNames(1, i);
        % compare each of the sigmas with the one in the file name 
        for i = 1:length(sigs)
            % once the matching sigma has been found, save the file name
            % into the respective position in the new file name array
            if contains(file_of_interest, sigs(i));
            sort_fileNames(1, i) = file_of_interest;
            end 
        end 
    end 
end 

%% Plot the spike count comparison

% Get the spike count of sigma 0

fileLoad_sig0 = num2str(cell2mat(sort_fileNames(1, 1))); % converts the file name for loading
load(fileLoad_sig0);

spikeCount_sig0 = a_cc.arrMIcore{1, 1}.x;  

% Possible y labels for the sigmas 
sig_labels = {'sigma 0 s', 'sigma 0.0001 s', 'sigma 0.005 s', 'sigma 0.02 s', 'sigma 0.5 s', 'sigma 1.5 s', 'sigma 30 s', 'sigma 60 s', 'sigma 90 s' , 'sigma 900 s'};

% List of all possible sigma strings
sig_names = {'sig0', 'sig00001', 'sig0005', 'sig002', 'sig05', 'sig15', 'sig30', 'sig60', 'sig90', 'sig900'};  

for i = 1:length(sort_fileNames)

    fileLoad = num2str(cell2mat(sort_fileNames(1, i))); 
    load(fileLoad);

    spikeCount_sigx = a_cc.arrMIcore{1, 1}.x;  

    x = spikeCount_sig0;
    x1 = spikeCount_sigx; 

    sig_label = num2str(cell2mat(sig_labels(1, i))); 

    figureName = ['spikeCount_comparison_' 'sig0_' num2str(cell2mat(sig_names(1, i))) '.png'];

    figure()

    scatter(x+normrnd(0, 0.1, size(x, 1), size(x, 2)), ...
    x1+normrnd(0, 0.1, size(x1, 1), size(x1, 2)), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.5); 

    title('Spike Count Comparison'); 
    xlabel('sigma 0 s'); 
    ylabel(sig_label);
    xlim([-10 50])
    ylim([-10 50])

    saveas(gcf, figureName)

end 

