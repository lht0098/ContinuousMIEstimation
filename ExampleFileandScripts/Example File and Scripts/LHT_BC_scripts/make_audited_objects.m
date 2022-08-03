%% Create a new analsis obejct with audited Ks

% Instuctions: RUN ONE SECTION AT A TIME!

clear all 

% get the name of the mat file that will be audited 
str_object = '20200427_bl21lb21_06032020data_unitB_mu0_sig1800_07222022unitBunitG_analysis_t2c1.mat';

% new string that will be added to new mat file
str_audit = '-audited_LHT';

% list of the count and timing pairs 
pairs = {'c1c2', 't1c2', 't2c1', 't1t2'};

% list of possible objects
objects = {'a_cc', 'a_tc', 'a_tc', 'a_tt'};

% add the new string to the name of the mat file right after the name of
% the pairing
if ~isempty(strfind(str_object, num2str(cell2mat(pairs(1)))))
    str_newObject = insertAfter(str_object, num2str(cell2mat(pairs(1))), str_audit)   
    object = num2str(cell2mat(objects(1)));
elseif ~isempty(strfind(str_object, num2str(cell2mat(pairs(2)))))
    str_newObject = insertAfter(str_object, num2str(cell2mat(pairs(2))), str_audit)   
    object = num2str(cell2mat(objects(2)));
elseif ~isempty(strfind(str_object, num2str(cell2mat(pairs(3)))))
    str_newObject = insertAfter(str_object, num2str(cell2mat(pairs(3))), str_audit)   
    object = num2str(cell2mat(objects(3)));
elseif ~isempty(strfind(str_object, num2str(cell2mat(pairs(4)))))
    str_newObject = insertAfter(str_object, num2str(cell2mat(pairs(4))), str_audit)   
    object = num2str(cell2mat(objects(4)));
end 

% load the old mat file 
load(str_object);

% make a copy of the old mat file with the new string name
copyfile(str_object, str_newObject);

% load the copy mat file for auditing 
load(str_newObject)

%% Save the changes made to the object

save(str_newObject, object);

