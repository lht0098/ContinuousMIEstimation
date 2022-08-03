clear all

load 20200427_bl21lb21_04052022unitDunitG_analysis_t1c2.mat

% flags to test different parts of the procedure;
% comment is rachel's initial values
f_dataFrac_triu = true; % true
f_dataFrac_thresh = 1; % 1
f_kStab_triu = true; % true
f_kStab_thresh = 1; % 1
f_stableKs_only = false; % false, comparing against all k values 
f_zeroStd_thresh = 1; % 1

core_object = a_tc.arrMIcore{3, 1};  

x = core_object.data_frac_stab_mat; 

z = core_object.k_val_stab_mat;  

ks = core_object.k_values;

weighted_k = zeros(1,size(x,3));

for i = 1:size(x,3)

    dataFrac_stab_matrix = x(:,:,i);

    if f_dataFrac_triu
        stab_mat = triu(dataFrac_stab_matrix);
    else
        stab_mat = dataFrac_stab_matrix; 
    end 
    
    % First, assess stability of each data frac comparison
    % 20220331 LHT: stability_boolean = stab_mat < 1;
    stability_boolean = stab_mat < f_dataFrac_thresh;
    
    % Determine if the frist four data fracs are stable
    if all(stability_boolean(1:4,1:4))
        stability_weight = 1;
    else
        stability_weight = 0;
    end
    
    % Add to the weight based on the total number of stable data
    % fracs
    n_stableFracs = sum(sum(stability_boolean));
    
    % Find total number of comparisons
    nFracPairs = numel(dataFrac_stab_matrix);
    
    % Add stability number to weight
    stability_weight = stability_weight + n_stableFracs/nFracPairs;

    weighted_k(i) = stability_weight;

end 

notBad_Ks = find(weighted_k >= 1);

stable_Ks = notBad_Ks; 

final_MIs = zeros(size(ks));
final_stds = zeros(size(ks));

% Only check for stability across ks if there are more than 1
% stable k values. Otherwise, audit is necessary
if length(notBad_Ks) > 1

    % Find matrix to describe stability across good k values
    if all(weighted_k < 1)
        k_stab = [];
    else
        k_stab = z(:,:,1);
    end
    % Get stability metric value for k stability
    if all(weighted_k < 1)
        k_stability_weights = zeros(size(ks));
    else
        % 20220406 LHT and Bryce: k_stability_weights = test_k_stab(obj, ks, notBad_Ks, k_stab);
        k_stability_matrix = k_stab;
        
        if f_kStab_triu
            % Focus only on the upper diagonal to compare forwards only
            stab_mat = triu(k_stability_matrix);
        else 
            stab_mat = k_stability_matrix;
        end
        
        if ~f_stableKs_only
            % First assess stability of each k value comparison
            stability_boolean = stab_mat < f_kStab_thresh; % Should this be only stable ks? 
        else 
             stab_mat_mask = zeros(length(notBad_Ks));
             for i = 1:length(notBad_Ks)
                 for j = 1:length(notBad_Ks)
                     stab_mat_mask(i,j) = stab_mat(notBad_Ks(i), notBad_Ks(j));
                 end
             end 
             stability_boolean = stab_mat_mask < f_kStab_thresh;
        end 

        % Make k stability weight placeholder vector
        k_stability_weights = zeros(size(ks));
        
        % stable_Ks = notBad_Ks; 

        % compare condition of requiring stability across all k's versus
        % versus requiring stability across stable k's 

        % Determine whether all ks are within error of each other
        if all(stability_boolean) % Should this be only stable ks? 
            k_stability_weights(stable_Ks) = 1;
        else
            % If all ks are not within error of each other, check if
            % they are within error of zero
            %final_MIs = zeros(size(stable_Ks));
            %final_stds = zeros(size(stable_Ks));
            %for ik = 1:length(stable_Ks)
            for ik = 1:length(ks)
                % Run getMIs to return the raw estimated values for all possible k-values
                k = ks(ik);
                %k = stable_Ks(ik); %20220408 LHT and BC: correcting the index
                r = core_object.get_singlek_mi(k);
                ifinal_MI = r.mi;
                final_MIs(ik) = ifinal_MI;
                ifinal_std = r.err;
                final_stds(ik) = ifinal_std;
            end
            
            if ~f_stableKs_only
                if all(final_MIs - final_stds * f_zeroStd_thresh <= 0) 
                    k_stability_weights(stable_Ks) = 1;
                end
            else 
                % if all(final_MIs - final_stds <= 0)
                if all(final_MIs([stable_Ks]) - final_stds([stable_Ks]) * f_zeroStd_thresh <= 0) 
                    k_stability_weights(stable_Ks) = 1;
                end
            end
        end
        
        % Return k stability weight vector
        % r = k_stability_weights;
    end
    % Get new weights for ks according to k stability matrix
    weighted_k = weighted_k + k_stability_weights;
    
else
    k_stab = [];
    
end



% Return all values to core object
opt_k = {final_MIs, final_stds, weighted_k, ['k < 1: Bad; 1 <= k < 2: Not Bad (Ks have data fraction stability, but do not match each othebr, see matrix); k > 2: Good! Ks in this range match and have consistent ' ...
                        'data fractions!'], k_stab};

valid_ks = opt_k{1,3};
MIs = opt_k{1,1};
errs = opt_k{1,2};

% 20220413 LHT/BC
% if sum(valid_ks > 2) == 0
if sum(valid_ks >= 2) == 0
   % 20220413 LHT/BC
   % if sum(valid_ks > 2) == 0
   if sum(valid_ks >= 1) == 0
       %error('Error: No k values have stable data fractions. Please manually select a k')
       % 20220413 LHT/BC
       % warning('NO k values have stable data fractions. Please manually select a k. FOR NOW- selecting minimum k with max stability that minimizes error')
       warning('The first 4 data fractions are not stable for any k. Please manually select a k. FOR NOW- selecting minimum k with max stability that minimizes error')
       best_weight = max(valid_ks);
       min_errIdx = find(errs == min(errs(valid_ks == best_weight)));
       % RC 20200310: FOR NOW, we just return an NaN value
       if isempty(min_errIdx)
           % 20220413 LHT/BC
           % best_kIdx = 1;
           best_kIdx = 0;
           % Place holder to get past test
           
           % 20220413 LHT/BC
           % MI = MIs(1);
           MI = NaN; 

           err = NaN;
           warning('Cannot minimize error with maximum stability.')
       else
           best_kIdx = min(min_errIdx); % Need to check that ks are monotonically increasing
           MI = MIs(best_kIdx);
           err = errs(best_kIdx);
       end

   else
       % Choose k with maximum stability that minimizes stdev
       % 20220413 LHT/BC
       % warning('Warning: K values are stable, but do not have consistent estimates. Selecting minimum k with maximum stability that minimizes error. Audit recommended.')
       warning('Warning: At least one K value has stable data fractions, but MI is not consistent across stable K values. Selecting minimum k with maximum stability that minimizes error. Audit recommended.')
       % 20220413 LHT/BC
       % best_weight = max(valid_ks(valid_ks > 1));
       best_weight = max(valid_ks(valid_ks >= 1));
       min_errIdx = find(errs == min(errs(valid_ks == best_weight)));
       % RC 20200310: FOR NOW, we just return an NaN value
       if isempty(min_errIdx)
           % 20220413 LHT/BC
           % best_kIdx = 1;
           best_kIdx = 0;
           % Place holder to get past test
           
           % 20220413 LHT/BC
           % MI = MIs(1);
           MI = NaN; 

           err = NaN;
           warning('Cannot minimize error with maximum stability.')
       else
           best_kIdx = min(min_errIdx);
           MI = MIs(best_kIdx);
           err = errs(best_kIdx);
       end

   end
else
   disp('At least one k value has stable data fractions and is consistent with other ks')
   % Choose minimum k with maximum stability that minimizes stdev
   % 20220413 LHT/BC
   % best_weight = max(valid_ks(valid_ks > 2));
   best_weight = max(valid_ks(valid_ks >= 2));
   min_errIdx = find(errs == min(errs(valid_ks == best_weight)));
   % RC 20200310: FOR NOW, we just return an NaN value
   if isempty(min_errIdx)
           % 20220413 LHT/BC
           % best_kIdx = 1;
           best_kIdx = 0;
           % Place holder to get past test
           
           % 20220413 LHT/BC
           % MI = MIs(1);
           MI = NaN; 

           err = NaN;
           warning('Cannot minimize error with maximum stability.')
   else
       best_kIdx = min(min_errIdx);
       MI = MIs(best_kIdx);
       err = errs(best_kIdx);
       
       %k_vals = [obj.mi_data{:,4}];
       %ks = unique(k_vals);
       %if ks(best_kIdx) == 1 %& MI - err >= 0
%        if best_kIdx == 0 
%            warning('Warning: K values are stable and consistent, but k = 1 was selected. Audit recommended')
%        else valid_ks(best_kIdx) == 1  %& MI - err >= 0
%            warning('Warning: K values are stable and consistent, but k = 1 was selected. Audit recommended')
%        end
   end     

end

% Output k value that was selected
% 
%k_vals = [obj.mi_data{:,4}];
%ks = unique(k_vals);

ks(best_kIdx)
MI 
err  

% Purpose: Select the final k value

