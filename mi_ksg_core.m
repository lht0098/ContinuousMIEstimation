classdef mi_ksg_core < handle
    % MI_KSG_core is used to set up sets of simulations to determine an
    % optimum k-value for a mutual information estimate and also calculates
    % mutual information and estimation error.
    
    properties
        verbose % for debugging purposes
        
        x % matrix of x data
        
        y % matrix of y data
        
        k_values % array of k-values
        
        mi_data % MI value, error estimate, data fraction, k-value
        
        opt_k % optimized k value; if -1, only runs MI calculation without any error estimate
        
        data_fracs = 10 % number of data fractions

        %%% 20220404 LHT 
        data_frac_stab_mat % data fraction stability matrix for auditing

        %%% 20220404 LHT
        k_val_stab_mat % k value stability matrix for auditing 
        
        append % Taken from analysis object. 
               % Specify whether to re-run all analysis or to just run analysis for k-values that have not been previously included. 
        
        sim_obj % sim_manager object
        
        set_seed = false % Can manually set to true for troubleshooting if want. 
        
        analysis_failure = 'NONE' % Defaults to NONE. Value changed if NaN value is returned for overall MI. 
    end
    
    methods
        function obj = mi_ksg_core(sim_obj, x, y, varargin)
            % This function generates the core object with the x and y
            % values, desired k values to run, and whether to estimate the error. 
            
            % Instantiate input parser
            p = inputParser;
            
            % Set up required inputs
            validate_sim_obj = @(x_var) assert(isa(x_var, 'mi_ksg_sims'), 'sim_obj must be a valid sim object');
            p.addRequired('sim_obj', validate_sim_obj);
            
            validate_x = @(x_var) assert(ismatrix(x_var), 'x must be a matrix');
            p.addRequired('x', validate_x);
            
            validate_y = @(x_var) assert(ismatrix(x_var), 'y must be a matrix');
            p.addRequired('y', validate_y);
            
            % Add optional
            
            % ks_arr
            default_ks_arr = 1:9;
            validate_ks_arr = @(x_var) assert(ismatrix(x_var), 'ks_arr must be a vector of integers');
            p.addParameter('ks_arr', default_ks_arr, validate_ks_arr);
            
            % opt_k
            valid_opt_k = [0 1 -1];
            default_opt_k = 1;
            validate_opt_k = @(x_var) assert(ismember(x_var, valid_opt_k), 'opt_k must be 0, 1 or -1');
            p.addParameter('opt_k', default_opt_k, validate_opt_k);
            
            % append
            default_append= true;
            validate_append = @(x_var) assert(islogical(x_var), 'append must be a logical value');
            p.addParameter('append', default_append, validate_append);
            
            % verbose
            default_verbose = 1; 
            validate_verbose = @(x_var) assert(isnumeric(x_var) && (rem(x_var,1) == 0), 'verbose must be an integer');
            p.addParameter('verbose', default_verbose, validate_verbose);
            
            % Parse the inputs
            p.KeepUnmatched = 1;
            p.parse(sim_obj, x, y, varargin{:});
            
            obj.x = p.Results.x;
            obj.y = p.Results.y;            
            obj.sim_obj = p.Results.sim_obj;
            obj.k_values = p.Results.ks_arr;
            obj.opt_k = p.Results.opt_k;
            obj.append = p.Results.append;
            obj.verbose = p.Results.verbose;
            
            add_sim(sim_obj, obj); % add this core obj to sim_manager list
        end
        
        % These methods are used to interface with other classes for data
        % analysis and visualization        
        function r = get_core_dataset(obj)
            
            
            % get cell array of data for MI calculation
            r = cell(0,4);
            if obj.opt_k < 0
                % only run MI calculation without error estimate
                for i=1:length(obj.k_values)
                    while 1
                        % generate unique key to track each simulation
                        key = num2str(dec2hex(round(rand(1)*100000)));
                        if ~any(strcmp(r(:,4), key))
                            break;
                        end
                    end
                    r = cat(1, r, {obj.x obj.y obj.k_values(i) key});
                end
            else
                % run MI calculation with error estimates
                 if obj.append
                    if isempty(obj.mi_data)
                        % Run estimates for all k values if none have been
                        % run yet
                        %Temporary to make the data fraction set
                        %predictable across runs. 
                        % RC 20200406: We will need to decide if we want to
                        % set the seed prior to running estimates and where
                        % we want to do that. 
                        if obj.set_seed
                            rng(0);
                        end
                        for i = 1:length(obj.k_values)
                            % create datasets for data fractions with unique key
                            % to track each simulation
                            r = cat(1, r, fractionate_data(obj, obj.k_values(i)));
                        end
                    else
                        mi_data = cell2mat(obj.mi_data);
                        k_finished = unique(mi_data(:,4));
                        for i = 1:length(obj.k_values)
                            if ismember(obj.k_values(i), k_finished)
                                continue
                            else
                                % create datasets for data fractions with unique key
                                % to track each simulation
                                r = cat(1, r, fractionate_data(obj, obj.k_values(i)));
                            end
                        end
                    end
                else
                    for i=1:length(obj.k_values)
                        % create datasets for data fractions with unique key
                        % to track each simulation
                        r = cat(1, r, fractionate_data(obj, obj.k_values(i)));
                    end
                end
            end 
        end
        
        function set_core_data(obj, dataset)
            % take sim_manager MI calculations and process results
            
            data_keys = unique([dataset(:,3)]); % extract simulation keys
            
            tmp_mi_data = cell(0,4);

            for key_ix = 1:length(data_keys) % iterate through each MI error estimation set
                tmp_match = strcmp([dataset(:,3)], data_keys(key_ix)); % find MI calculations that correspond to same data fractions
                count = sum(tmp_match); % determine number of data fractions
                data_ixs = find(tmp_match == 1); % identify which simulations to include
                               
                % Check to see if there are multiple k-values for this set
                ks = unique(cell2mat(dataset(data_ixs,2)));
                if length(ks) > 1
                    % We have duplicate IDs and need to resolve
                    for ik = 1:length(ks)
                        k_idx = cell2mat(dataset(data_ixs,2)) == ks(ik);
                        k_data_ixs = data_ixs(k_idx);
                        % Find MI and k for ID specific to this k value
                        mi = [dataset{k_data_ixs,1}];
                        k = dataset{k_data_ixs(1),2};
                        
                        % Find new count for ID specific to this k value
                        count = sum(k_idx);
                        
                        % Append for each k within the duplicate IDs
                        tmp_mi_data = cat(1, tmp_mi_data, {mean(mi) var(mi) count k}); % append MI with error estimation
                        
                    end

                else
                    % Set mi and k for this ID
                    mi = [dataset{data_ixs,1}];
                    k = dataset{data_ixs(1),2};
                    tmp_mi_data = cat(1, tmp_mi_data, {mean(mi) var(mi) count k}); % append MI with error estimation
                end

            end
            
            if obj.append
                tmp_mi_data = cat(1, tmp_mi_data, obj.mi_data);
            end
            
            obj.mi_data = sortrows(tmp_mi_data,[4,3]);
            
        end
        
        function r = get_mi(obj, errThreshold, varargin)
            
            p = inputParser;
            
            % Set required inputs
            validate_obj = @(x) assert(isa(x, 'mi_ksg_core'), 'obj must be a core object');
            p.addRequired('obj', validate_obj);
            
            validate_errThreshold = @(x) assert(isnumeric(x), 'errThreshold must be a numeric quantity');
            p.addRequired('errThreshold', validate_errThreshold);
            
            % Set Optional Input
            default_k = 0;
            validate_k = @(x) assert(rem(x,1) == 0, 'k value must be a whole number, divisable by 1');
            p.addParameter('k', default_k, validate_k);
            
            % Parse inputs
            p.parse(obj, errThreshold, varargin{:});
            
            % Define validated inputs
            obj = p.Results.obj;
            errThreshold = p.Results.errThreshold;
            k = p.Results.k;

            if isempty(obj.mi_data)
                disp('Data was not able to be analyzed. Returning NaN value for MI');
                r.mi = NaN;
                r.err = NaN;
                r.k = NaN;
                
                return

            else
                disp('Analyzing data');
                 % Find MIs differently depending on value of k
                if k == 0
                    % Find MI for optimized k.
                    r = obj.get_opk_mi();

                   % Sanity check that our MI value matches what it would be if
                   % we had inputted the optimized k value

                   %GET VALUES FOR SANITY CHECK
                   % Find MI calcs with k-value
                    test_data_ixs = cell2mat(obj.mi_data(:,4)) == r.k;
                    if length(test_data_ixs) == 0
                        r.mi = NaN;
                        r.err = NaN;
                        r.k = NaN;
                        return;
                    end

                    if isnan(r.mi) | r.mi == 0; return; end      

                    % calculate estimated error
                    test_listSplitSizes = cell2mat(obj.mi_data(test_data_ixs,3));
                    test_MIs = cell2mat(obj.mi_data(test_data_ixs,1));
                    test_listVariances = cell2mat(obj.mi_data(test_data_ixs,2));

                    if sum(isnan(test_MIs) ~= 0)
                        if sum(~isnan(test_MIs)) >= 4
                            test_listSplitSizes = test_listSplitSizes(~isnan(test_MIs));
                            test_MIs = test_MIs(~isnan(test_MIs));
                            test_listVariances = test_listVariances(~isnan(test_listVariances));
                        end
                        % Check that sizes are all still consistent
                        if size(test_listSplitSizes) ~= size(test_MIs) | size(test_listSplitSizes) ~= size(test_listVariances)
                            error('Error: Sizes of vectors without NaN values do not match')
                        end
                    end


                    test_listVariances = test_listVariances(2:end);

                    test_k_err = test_listSplitSizes(2:end);
                    test_variancePredicted = sum((test_k_err-1)./test_k_err.*test_listVariances)./sum((test_k_err-1));

                    test_MI = test_MIs(1);
                    test_err = test_variancePredicted.^0.5;
                    r.mi
                    test_MI
                    r.err
                    test_err

                    % ACTUAL SANITY CHECK
                    if ~isequaln(r.mi,test_MI) | ~isequaln(r.err, test_err) 
                        keyboard
                        error('Optimized K MI does not match that stored in mi_data'); 
                    end


                else
                    % Find MI for specified k.
                    r = obj.get_singlek_mi(k);   
                    

                end
% -------------THIS CALCULATES THE ERROR OF THE ERROR-----------------------------------            
%             N = size(obj.x,2);
%             Sml=variancePredicted*N;
%             varS = 2*Sml^2/sum((k-1)); %Estimated variance of the variance of the estimate of the mutual information at full N
%             stdvar = sqrt(varS/N^2); %the error bars on our estimate of the variance
% --------------------------------------------------------------------------------------

                % 2019107 BC
                % Adding hack to filter mutual information results that are
                % within three S.D. from 0
                if errThreshold < 0
                    r.mi = r.mi;
                    r.err = r.err;

                elseif errThreshold == 0
                    if MI >= 0
                        r.mi = r.mi;
                    else
                        r.mi = 0;
                    end
                    r.err = r.err;

                elseif errThreshold > 0
                    if (MI - errThreshold*(err)) > 0 || ((errThreshold == 0) && (MI > 0))
                        r.mi = r.mi;
                        r.err = r.err;
                    else
                        r.mi = 0;
                        r.err = r.err;
                    end      
                end  
            end            
        end
        
        function find_k_value(obj)
            % determine best k-value to use
            data = cell2mat(obj.mi_data);
            
            k_vals = [obj.mi_data{:,4}];
            ks = unique(k_vals);
            weighted_k = zeros(size(ks)); 

            % 20220404 LHT L331-L338: 
            % initialize matrix to store all data fraction 
            % stability matrices
            num_dataFracs = data(find(data(:,4) == 1),3);
            obj.data_frac_stab_mat = zeros(size(num_dataFracs,1), size(num_dataFracs,1), size(ks,2));
             
            for ik = 1:length(ks)
                % Find data fraction stability matrix for each k
                dataFrac_stab = get_stabMat_dataFrac(obj,ks(ik));

                % put the matrix in the next page
                obj.data_frac_stab_mat(:,:,ks(ik)) = dataFrac_stab;
                
                % Get stability metric value for each k
                weighted_k(ik) = test_dataFrac_stab(obj, dataFrac_stab);               
            end  

            % Identify k values that satisfy the minimum criteria for
            % stability
            notBad_Ks = find(weighted_k >= 1);
            
            %%% get MI and std for k neighbor eval and to return to core
            %%% object
            final_MIs = zeros(size(ks));
            final_stds = zeros(size(ks));

            for ik = 1:length(ks)
                k = ks(ik);
                r = get_singlek_mi(obj, k);
                final_MIs(ik) = r.mi;
                final_stds(ik) = r.err;
            end

            % Only check for stability across ks if there are more than 1
            % stable k values. Otherwise, audit is necessary
            if length(notBad_Ks) > 1

                % Find matrix to describe stability across good k values
                %%% initiate array of k stab values for all ks  
                k_stab = get_stabMat_kvals(obj, ks);

                %%% 20220404 LHT L376-384: K stability matrices for auditing
                %%% initialize matrix to store all k stability matrices
                obj.k_val_stab_mat = zeros(size(ks,2), size(ks,2));
                %%% save the matrix
                obj.k_val_stab_mat = k_stab;

                % Get stability metric value for k stability
                %%% initiate array of k stability weights 
                k_stability_weights = test_k_stab(obj, ks, notBad_Ks, k_stab);

                % Get new weights for ks according to k stability matrix
                weighted_k = weighted_k + k_stability_weights;
                
            end

            %%% k neighbor stability evaluation
            if all(weighted_k < 2)
                k_stab_neighbor_i = find(weighted_k > 1); % get the k vals that have stable data fractions
                k_stab_neighbor = zeros(2,length(k_stab_neighbor_i)); % initiate an array to hold the std of both neighbors
                for i=1:length(k_stab_neighbor_i) 
                    if k_stab_neighbor_i(i)+1 <= length(final_MIs)
                        k_stab_neighbor(1,i) = abs(final_MIs(k_stab_neighbor_i(i)+1) - final_MIs(k_stab_neighbor_i(i))) / final_stds(k_stab_neighbor_i(i)); 
                    end
                    if k_stab_neighbor_i(i)-1 >= 1
                        k_stab_neighbor(2,i) = abs(final_MIs(k_stab_neighbor_i(i)-1) - final_MIs(k_stab_neighbor_i(i))) / final_stds(k_stab_neighbor_i(i));
                    end 
                end 

                % BC: check that stability is within threshold criteria
                k_stab_neighbor_comp = zeros(1,size(k_stab_neighbor,2)); % initate an array to hold thresh eval
                k_stab_neighbor_op = zeros(1,size(k_stab_neighbor,2)); % initiate an array to hold stability eval 
                for i=1:size(k_stab_neighbor,2)
                    k_stab_neighbor_comp(i) = k_stab_neighbor(1,i) < 1.5 & k_stab_neighbor(2,i) < 1.5; % check that both neighbors are within the threshold
                    % BC: find k value that optimizes stability with neighbors
                    if k_stab_neighbor_comp(i) == 1 % check that both neighbors are within thresh
                        k_stab_neighbor_op(i) = mean(k_stab_neighbor(:,i)); % evaluate the average std
                    else 
                        k_stab_neighbor_op(i) = NaN;
                    end 
                end 

                if all(isnan(k_stab_neighbor_op))
                    best_neighStab = 0;
                else 
                    [~, lowest_err] = min(k_stab_neighbor_op); % get the index of the k that has the pair with the lowest avg std
                    best_neighStab = k_stab_neighbor_i(lowest_err); % get the k val that has the highest stability with neighbors
                end 
            %%% set the best_neighStab to 0 if there is a k weight => 2
            else
                best_neighStab = 0;
            end
            
            % Return all values to core object including best_neighStab
            obj.opt_k = {final_MIs, final_stds, weighted_k, ['k < 1: Bad; 1 <= k < 2: Not Bad (Ks have data fraction stability, but do not match each othebr, see matrix); k > 2: Good! Ks in this range match and have consistent ' ...
                                    'data fractions!'], best_neighStab};
        end
        
        function r = get_opk_mi(obj)
                % Find MIs differently depending on value of k
                % We have either already optimized k, or we need to
                % optimize k  
                if ~iscell(obj.opt_k)
                    % k has not been optimized yet. We need to run
                    % find_k_value
                    obj.find_k_value();
                end
               valid_ks = obj.opt_k{1,3};
               MIs = obj.opt_k{1,1};
               errs = obj.opt_k{1,2};
               %%% get best_neighStab
               best_neigh = obj.opt_k{1,5};

               %%% k evaluation
               best_weight = max(valid_ks);
               min_errIdx = find(errs == min(errs(valid_ks == best_weight)));
               best_kIdx = min(min_errIdx);
               MI = MIs(best_kIdx);
               err = errs(best_kIdx);

               %%% k selection
               if any(valid_ks >= 2) 
                   disp('At least one k value has stable data fractions and is consistent with other ks')
               elseif any(valid_ks >= 1) && best_neigh == 0
                   warning('Warning: At least one K value has stable data fractions, but MI is not consistent across stable K values. Selecting minimum k with maximum stability that minimizes error. Audit recommended.')
                   MI = NaN;
                   err = NaN; 
               elseif any(valid_ks >= 1) && best_neigh > 0    
                   best_kIdx = best_neigh;
                   MI = MIs(best_neigh);
                   err = errs(best_neigh);
               else    
                   warning('The first 4 data fractions are not stable for any k. Please manually select a k. FOR NOW- selecting minimum k with max stability that minimizes error')
                   MI = NaN;
                   err = NaN;
                   %%% if all k's do not have stable data fractions, select
                   %%% the k with the max stab and min err
               end 
               
               % Output k value that was selected
               % 
               k_vals = [obj.mi_data{:,4}];
               ks = unique(k_vals);
               
               % Set return values
               r.k = ks(best_kIdx);
               r.mi = MI;
               r.err = err;
              
        end
        
        function r = get_singlek_mi(obj,k)
                % get mutual information and error estimates
                data_ixs = cell2mat(obj.mi_data(:,4)) == k; % find MI calcs with k-value

                % calculate estimated error
                listSplitSizes = cell2mat(obj.mi_data(data_ixs,3));
                MIs = cell2mat(obj.mi_data(data_ixs,1));
                listVariances = cell2mat(obj.mi_data(data_ixs,2));

                
                if sum(isnan(MIs) ~= 0)
                    if sum(~isnan(MIs)) >= 4
                        listSplitSizes = listSplitSizes(~isnan(MIs));
                        MIs = MIs(~isnan(MIs));
                        listVariances = listVariances(~isnan(listVariances));
                    end
                    % Check that sizes are all still consistent
                    if size(listSplitSizes) ~= size(MIs) | size(listSplitSizes) ~= size(listVariances)
                        error('Error: Sizes of vectors without NaN values do not match')
                    end
                end
                
                listVariances = listVariances(2:end);
                k_err = listSplitSizes(2:end);
                variancePredicted = sum((k_err-1)./k_err.*listVariances)./sum((k_err-1));
                
                MI = MIs(1);
                err = variancePredicted.^0.5;
                
                % Set return values
                r.mi = MI;
                r.k = k;
                r.err = err;
        end
       
        
        function r = get_stabMat_dataFrac(obj,k)
            
                data = cell2mat(obj.mi_data);
                
                dataFracs = data(find(data(:,4) == k),3);
                
                dataFrac_data = data(find(data(:,4) == k), 1:end);
                
                stability_matrix = zeros(length(dataFracs), length(dataFracs));
            for iFrac_row = 1:length(dataFracs)
                % Find MI dataFrac row number
                iMI = dataFrac_data(iFrac_row, 1);
                
                % Iterate through each column for this row
                for jFrac_column = 1:length(dataFracs)
                    jMI = dataFrac_data(jFrac_column,1);
                    jvar = dataFrac_data(jFrac_column, 2);
                    jstd = jvar ^ 0.5; 
                    if jFrac_column == 1
                        % Fill in with zero automatically
                        stability_matrix(iFrac_row, jFrac_column) = 0;
                    else                    
                        % Test that MI is within 1 stdev of others
                        stability_matrix(iFrac_row, jFrac_column) = ((iMI - jMI)/jstd)^2;
                    end
                end   
            end
            % Output stability matrix
            r = stability_matrix;                     
        end
        
        function r = test_dataFrac_stab(obj, dataFrac_stab_matrix)
            % For now we are keeping the 4 stable data frac minimum and adding to the weight based on the number of consistent data fracs
            % This function returns a data fraction stability weight for
            % a single k value
            
            % Focus only on the upper diagonal to compare forwards only
            stab_mat = dataFrac_stab_matrix;
            
            % First, assess stability of each data frac comparison
            stability_boolean = stab_mat < 1.5;
            
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
            
            % Return stability weight
            r = stability_weight;
        end
        
        function r = get_stabMat_kvals(obj, ks)

            % Fill in stability matrix for all data fraction stable k
            % values
            for ik_row = 1:length(ks)
                % Run getMIs to return the raw estimated values for all possible k-values
                ik = ks(ik_row);
                r = get_singlek_mi(obj, ik);
                ifinal_MI = r.mi;
                for jk_column = 1:length(ks)
                    jk = ks(jk_column);
                    r = get_singlek_mi(obj, jk);
                    jfinal_MI = r.mi;
                    jfinal_std = r.err;
                    
                    % Input values into matrix
                    k_stability_matrix(ik_row, jk_column) = ((ifinal_MI - jfinal_MI)/jfinal_std)^2;
                end
            end
            r = k_stability_matrix;
        end
        
        function r = test_k_stab(obj, ks, stable_Ks, k_stability_matrix)
            
            % Focus only on the upper diagonal to compare forwards only
            %%% look at the full matrix
            stab_mat = k_stability_matrix;
            
            %%% only evaluate stable k's 
            stab_mat_mask = zeros(length(stable_Ks));
            for i = 1:length(stable_Ks)
                for j = 1:length(stable_Ks)
                    stab_mat_mask(i,j) = stab_mat(stable_Ks(i), stable_Ks(j));
                end
            end 

            % First assess stability of each k value comparison
            %%% assess the stability of the stable k's at higher k stability threshold
            stability_boolean = stab_mat_mask < 1.5;
            
            % Make k stability weight placeholder vector
            k_stability_weights = zeros(size(ks));

            %%% initiate quality k's var for evaluation
            quality_Ks = stable_Ks(any(stability_boolean.*stability_boolean'-eye(size(stability_boolean))));

            %%% determine if the stable k's are within error of each other
            if length(quality_Ks) >= 1
                disp('K values have stable data fractions and are stable across ks')
                k_stability_weights(quality_Ks) = 1;
            end 
            
            % Return k stability weight vector
            r = k_stability_weights;
        end
        
        function r = fractionate_data(obj, k)
            % return cell array of fractionated datasets with x-data,
            % y-data, k-value, and ix
            
            n = length(obj.x);
            r = cell(sum(1:obj.data_fracs),4);
            
            for frac_n = 1:obj.data_fracs
                % determine length of subsample
                a = randperm(n);
                l = round(linspace(0,n,frac_n+1));
                
                % generate unique key to track each simulation
                while 1
                    key = num2str(dec2hex(round(rand(1)*100000))); 
                    if ~any(strcmp(r(:,4), key))                       
                        break;
                    end
                end
                
                % select subsample of data and assign data and params to data cell array
                for j=1:frac_n 
                    xT = obj.x(a(l(j)+1:l(j+1)));
                    yT = obj.y(a(l(j)+1:l(j+1)));
                    r(sum(1:(frac_n-1))+j,:) = {xT yT k key};
                end
            
            end
        end
        
    end
end
