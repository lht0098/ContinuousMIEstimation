%% Sort the file names from lowest to greatest sigma

clear all

cd % get the current directoy
myDir = cd; % set the current dir as the dir of interest
myFiles = dir(fullfile(myDir,'*.mat')); % puts all of the files into a struct
fileNames = repmat({'NaN'}, 1, length(myFiles)); % saves files for audit

% Create an array that will hold the mi and err 
mi_and_err = NaN(2, length(myFiles));

% List of all possible unit pairings 
unit_pairings = {'unitDunitE_analysis_t1c2', 'unitDunitE_analysis_t2c1', 'unitEunitG_analysis_t2c1', 'unitEunitG_analysis_t1c2', 'unitBunitD_analysis_t1c2', 'unitBunitD_analysis_t2c1', 'unitBunitE_analysis_t1c2', 'unitBunitE_analysis_t2c1', 'unitBunitG_analysis_t1c2', 'unitBunitG_analysis_t2c1', 'unitDunitG_analysis_t2c1', 'unitDunitG_analysis_t1c2'};

% Save the file names 
for i = 1:length(myFiles)
    fileNames(1,i) = cellstr(myFiles(i).name); % save the file name for auditing
end 

% Initiate an array to save the sorted file names
sort_fileNames = repmat({'NaN'}, 1, length(fileNames));

% Sort the files from the lowest to the greatest sigma
for i = 1:length(sort_fileNames)
    % check if the index of the file name and the simga matches
    if contains(fileNames(1, i), unit_pairings(i))
        sort_fileNames(1, i) = fileNames(1, i);
    else 
        % save the file name if they do not match
        file_of_interest = fileNames(1, i);
        % compare each of the sigmas with the one in the file name 
        for i = 1:length(unit_pairings)
            % once the matching sigma has been found, save the file name
            % into the respective position in the new file name array
            if contains(file_of_interest, unit_pairings(i));
            sort_fileNames(1, i) = file_of_interest;
            end 
        end 
    end 
end 

%% Find the mi and err for the analysis objects

% Go through each of the analysis objects and save the mutual information 
% and err of each one
for i = 1:length(sort_fileNames)

    fileLoad = num2str(cell2mat(sort_fileNames(1,i))); % converts the file name for loading
    load(fileLoad);

        if ~isempty(strfind(fileLoad,'c1c2'))
            data = cell2mat(struct2cell(a_cc.returnMIs));
            mi_and_err(1, i) = data(1,1); % input the mi into the 1st row 
            mi_and_err(2, i) = data(2,1); % input the err into the 2nd row
        elseif ~isempty(strfind(fileLoad,'t1c2'))
            data = cell2mat(struct2cell(a_tc.returnMIs));
            mi_and_err(1, i) = data(1,1); % input the mi into the 1st row 
            mi_and_err(2, i) = data(2,1); % input the err into the 2nd row
        elseif ~isempty(strfind(fileLoad,'t2c1'))
            data = cell2mat(struct2cell(a_tc.returnMIs));
            mi_and_err(1, i) = data(1,1); % input the mi into the 1st row 
            mi_and_err(2, i) = data(2,1); % input the err into the 2nd row
        elseif ~isempty(strfind(fileLoad,'t1t2'))
            data = cell2mat(struct2cell(a_tt.returnMIs));
            mi_and_err(1, i) = data(1,1); % input the mi into the 1st row 
            mi_and_err(2, i) = data(2,1); % input the err into the 2nd row
        end 

end 

%% Get the names of the units and the analysis object for naming

% Get the name of the current folder 
str = pwd; % get the str of the path
idx = strfind(str,'\'); % find all of the slashes of the path
folderName = str(idx(end)+1:end); % only get the string of the current folder
figureName = [folderName, '.png']; % set the name of the current folder as the name of the figure

% Get part of the folder name
part_folderName = extractAfter(folderName, 'jittered_');

% Units analyzed 
unit_jittered = extractBefore(folderName, '_jittered'); % extract the name of the unit that was jittered
units = {'unitB', 'unitD', 'unitE', 'unitG'};
units_of_interest = repmat({'NaN'}, 1, length(units)); % array to save the analyzed units

for i = 1:length(units)

    unit_str = num2str(cell2mat(units(i))); % get the unit 
    
    % Extract the unit pair that was analyzed 
    if ~isempty(strfind(part_folderName, unit_str));
       if unit_str == unit_jittered
          unit_str = append(unit_str, '*');
       end 
       units_of_interest(1, i) = cellstr(unit_str);
    end

end

% Find the indexes of the units
units_idx = ~(cellfun(@(units_of_interest) strcmp(units_of_interest, 'NaN'), units_of_interest)); 
units_title = num2str(cell2mat(units_of_interest(units_idx)));

% Possible analysis objects
pairs = {'C1C2', 'T1C2', 'T2C1', 'T1T2'};
pairs_of_interest = repmat({'NaN'}, 1, length(pairs));

for i = 1:length(pairs)

    pair_str = num2str(cell2mat(pairs(i))); % get the pair
    
    % Extract the pair that was analyzed 
    if ~isempty(strfind(part_folderName, pair_str));
       pairs_of_interest(1, i) = cellstr(pair_str);
    end

end

% Find the indexes of the pairs
pairs_idx = ~(cellfun(@(pairs_of_interest) strcmp( pairs_of_interest, 'NaN'), pairs_of_interest)); 
pairs_title = num2str(cell2mat(pairs_of_interest(pairs_idx)));

%% Plot the mi and err

x = 1:length(myFiles);

% List of all possible unit pairings 
units_label = {'DE', 'ED', 'GE', 'EG', 'BD', 'DB', 'BE', 'EB', 'BG', 'GB', 'GD', 'DG'};  

figure()

bar(x, mi_and_err(1,:))

hold on

errorbar(x, mi_and_err(1,:), mi_and_err(2,:), 'bp')

title('Timing-Timing')
xticklabels(units_label)
ylabel('Mutual Information (bits)')
ylim([0 0.4])

hold off

saveas(gcf, figureName)
