classdef mi_analysis < handle
    %  MI_KSG_data_analysis is a parent class used to set up a separate object for each pair of variables to feed into the MI calculations
    % 
    properties
        
        varNames % indicates which data(s) to pull
        
        objData % Reference to which data object to pull from
        
        objBehav % Reference to which behavior object to pull from (optional)
        
        % BC: This will be a list/cell array of objMIcore instances (may need to index)
        % cell array with structure:{{objMICore} {coeff} {k-value} {MIestimate} {Error} {coreID}}
        arrMIcore % Reference to MIcore object 
	
        sim_manager % Sim manager reference object
        
        append % Specify whether to re-run analysis or just for k-values that have not been previously included
        
        verbose % level of output for progress and troubleshooting/debugging
        
        notes %Indicates how much data has been omitted (optional)
        
        reparam %Do you reparameterize
        
        k_audited = 'No' % Have the k-values been manually audited? Automatically set to yes once auditing function is run.
 
        
        cycle_select
    end

    methods

        function obj = mi_analysis(objData, objBehav, varNames, varargin)
            % This funtion inputs the data object reference and variable references
            
            % Instantiate input parser
            p = inputParser;
            
            % Set up required inputs
            validate_objData = @(x) assert(isa(x, 'mi_data'), 'objData must be a valid data object');
            p.addRequired('objData', validate_objData);
            
            validate_objBehav = @(x) assert(isa(x, 'mi_data'), 'objBehav must be a valid behavior object');
            p.addRequired('objBehav', validate_objBehav);
            
            validate_varNames = @(x) assert(iscell(x), 'varNames must be a cell array of strings');
            p.addRequired('varNames', validate_varNames);
            
            % append
            default_append = true;
            validate_append = @(x) assert(islogical(x), 'append must be a logical value');
            p.addParameter('append', default_append, validate_append);
            
            % verbose
            default_verbose = 1;
            validate_verbose = @(x) assert(isnumeric(x) && rem(x,1) == 0, 'verbose must be an integer');
            p.addParameter('verbose', default_verbose, validate_verbose);
            
            % reparam
            default_reparam = false;
            validate_reparam = @(x) assert(islogical(x), 'reparam must be a logical value');
            p.addParameter('reparam', default_reparam, validate_reparam);
            
            
            % cycle_select
            default_cycle_select = -1;
            validate_cycle_select = @(x) assert(isnumeric(x), 'cycle_selection must be array of integers');
            p.addParameter('cycle_select', default_cycle_select, validate_cycle_select);
            
            % Parse the inputs
            % Set up InputParser to handle extra inputs from subclasses
            p.KeepUnmatched = 1;
            p.parse(objData, objBehav, varNames, varargin{:});
            
            obj.objData = p.Results.objData;
            obj.varNames = p.Results.varNames;
            obj.objBehav = p.Results.objBehav;
            obj.append = p.Results.append;
            obj.verbose = p.Results.verbose;
            obj.reparam = p.Results.reparam;
            obj.cycle_select = p.Results.cycle_select;
            
            % Temporarily set arrMIcore and instantiate a sim_manager
            % object
            obj.arrMIcore = {};
            obj.sim_manager = mi_ksg_sims(1,3);
            
            if obj.verbose > 0; disp([newline 'mi_analysis instantiated']); end
        end
        
        function set_parameters(obj,varargin)
            p = inputParser;           

            % append
            default_append = true;
            validate_append = @(x) assert(islogical(x), 'append must be a logical value');
            p.addParameter('append', default_append, validate_append);
            
            % verbose
            default_verbose = 1;
            validate_verbose = @(x) assert(isnumeric(x) && rem(x,1) == 0, 'verbose must be an integer');
            p.addParameter('verbose', default_verbose, validate_verbose);
            
            % reparam
            default_reparam = true;
            validate_reparam = @(x) assert(islogical(x), 'reparam must be a logical value');
            p.addParameter('reparam', default_reparam, validate_reparam);
            
            % cycle_select
            default_cycle_select = -1;
            validate_cycle_select = @(x) assert(isnumeric(x), 'cycle_selection must be numeric array');
            p.addParameter('cycle_select', default_cycle_select, validate_cycle_select);
            
            p.KeepUnmatched = 1;
            p.parse(varargin{:});
            
            obj.append = p.Results.append;
            obj.verbose = p.Results.verbose;
            obj.reparam = p.Results.reparam;
            obj.cycle_select = p.Results.cycle_select;
        end
        
        function buildMIs(obj, mi_data)

            % Set up empty array for obj.arrMIcore
            obj.arrMIcore = cell(size(mi_data,1),4);
            
            % Define groups to fill in arrMIcore
            xGroups = mi_data{1};
            yGroups = mi_data{2};
            coeffs = mi_data{3};
            
            v = obj.verbose;
            
            if v > 1; disp([newline 'Assigning ID to each subgroup...']); end
            
            for iGroup = 1:size(xGroups,1)
                x = xGroups{iGroup,1};                                                                                                                                                                                                                                    
                y = yGroups{iGroup,1};
                
                %Reparametrize data in each subgroup
                if obj.reparam == true
                    if ~all(rem(x,1) == 0) %Check for continuous data
                        x = reparameterize_matrix(x);
                    end
                    if  ~all(rem(y,1) == 0)
                        y = reparameterize_matrix(y);   
                    end
                end

              
                while 1 % generate random key to keep track of which MI calculations belong together
                    key = num2str(dec2hex(round(rand(1)*100000)));
                    % break the while loop if the key has not already been
                    % assigned.
                    if iGroup == 1
                        break
                    elseif ~ismember({obj.arrMIcore{1:end,4}}, key)
                        break
                    end
                end
                
                
                % RC: Why do we set the k values in the core object and in
                % the arrMIcore?
                core1 = mi_ksg_core(obj.sim_manager, x, y, 'ks_arr', 1:9, 'opt_k',1, 'append', obj.append, 'verbose', obj.verbose);
                if v > 2; disp([newline '--> Core object instantiated for group: ' num2str(iGroup)]); end
                
                if v > 2; disp([newline '--> Group ' num2str(iGroup) ' has ' num2str(max(size(x))) ' data points']); end
                
	            obj.arrMIcore(iGroup,:) = {core1 coeffs{iGroup,1} 0 key};
                
                if v > 2; disp([newline '--> arrMIcore assigned']); end
	            % BC: The obj.findMIs function basically calls run_sims
            end
            
            if v > 4; mi_ksg_viz.audit_plots(obj); end
            
	    % Sets up an MIcore object to calculate the MI values, and pushes the
	    % data from this object to the MIcore process. 
        
            if v > 0; disp(['COMPLETE: Added all data to arrMIcore']); end
        end
        
        function calcMIs(obj)
            v = obj.verbose;
            if v > 0; disp('Calculating mutual information...'); end
            run_sims(obj.sim_manager);
        end

        function getMIs(obj)
            for iCores = 1:size(obj.arrMIcore,1)
                core = obj.arrMIcore{iCores,1};
                disp(['MI core index: ' num2str(iCores)]);

                % Set the warning message to empty
                lastwarn('')

                % Call core function to find k value and output MI and error
                r = core.get_mi(-1);
                

                % Grab last warning message
                w = lastwarn;

                if isempty(w)
                    obj.arrMIcore{iCores,3} = r.k;
                else
                    obj.arrMIcore{iCores,3} = w;
                end

                % Set temporary MI and error values
                obj.arrMIcore{iCores,4} = r.mi;
                obj.arrMIcore{iCores,5} = r.err;
                
                % Log any reasons for analysis failures
                obj.arrMIcore{iCores,6} = core.analysis_failure;
                    
            end           
            
        end
        
        function r = returnMIs(obj)
            % NOTE WE NEED TO EDIT THIS AS WE MAKE DESIGN DECISIONS!!
            
            % Set all negative MI values to zero
            MIs = [];
            for iSubgroup = 1:size(obj.arrMIcore,1)
                MI = obj.arrMIcore{iSubgroup,4};
                if MI < 0
                    MI = 0;
                end
                MIs = [MIs; MI];
            end
            

            
            % Find Error Using Error Propagation Equations
            % We need to iterate through each subgroup to propagate error
            var_sum_vec = zeros(size(obj.arrMIcore, 1), 1);
            test_val = 0;
            if obj.discard_omittedData
                
                discard_reasons = unique(obj.arrMIcore(:,6));
                % NONE is not a discard reason because there was no
                % analysis failure for these cases.
                discard_reasons = discard_reasons(~strcmp(discard_reasons, 'NONE'));
                discard_reasons{end + 1} = 'Do NOT discard any data';
                
                % User prompt to select reason for 'NONE' response
                promptMessage = 'Please select which data you wish to discard from the probability term (press shift to select multiple options)';
               
                button = listdlg('PromptString', promptMessage, 'ListString', discard_reasons);

                if ~any(strcmp(discard_reasons(button), 'Do NOT discard any data')) 
                    
                    % Identify reasons to discard data:
                    discard = discard_reasons(button,:);
                    keep_flag = ~ismember(obj.arrMIcore(:,6),discard);
                    
                    obj.arrMIcore(:,7) = mat2cell(keep_flag, ones(size(obj.arrMIcore,1),1));
                    
                    MIcore = obj.arrMIcore(keep_flag,:); 
                    MIs = MIs(keep_flag);
                else
                    if length(button) > 1
                        warning('Do not discard data selected with other discard reasons. Proceeding without discarding data.')
                    end
                    MIcore = obj.arrMIcore;
                end
            else
                MIcore = obj.arrMIcore;
            end
            
            % Find weighted sum of MIs
            prob_vec = cell2mat(MIcore(:,2))/sum(cell2mat(MIcore(:,2)));
            % BC EDIT 20220803: matlab2022a deprecated
            % r.mi = sum(MIs.*prob_vec, 'omitnan');
            r.mi = nansum(MIs.*prob_vec);
            
            for iSubgroup = 1:size(MIcore, 1)
                iMI = MIcore{iSubgroup, 4};
                iprob = MIcore{iSubgroup, 2};
                %if obj.discard_omittedData
                    iprob = iprob/sum(cell2mat(MIcore(:,2)));
                %end
                iErr = MIcore{iSubgroup, 5};
                iVar = iErr^2;
                % Find total number of data points in this subgroup
                inx = size(MIcore{iSubgroup,1}.x, 2);
                ip_Var = ((1 - iprob)/(inx*iprob))*(iprob^2);
                
                var_sum_vec(iSubgroup, 1) = ( (iVar/iMI^2) + (ip_Var/iprob^2) )*iMI^2*iprob^2;
                test_val = test_val + iprob;

            end
            
	    % Sanity check
            if obj.discard_omittedData
                if ~ismembertol(test_val, 1, 1e-12); error('Error: Re-weighted probability terms do not sum to one.'); end
            end
            
            % LHT EDIT 20220803: matlab2022a deprecated
            % var_tot = sum(var_sum_vec, 'omitnan');
            var_tot = nansum(var_sum_vec);
            r.err = var_tot^0.5;
        end
        
        function auditKs(obj)
            % Go through each subgroup and determine if a k value was
            % automatically selected
            for iSubgroup = 1:size(obj.arrMIcore,1)
                k_val = obj.arrMIcore{iSubgroup,3};
                
                if ischar(k_val)
                    input_str = strcat('What K value is best for subgroup ', num2str(iSubgroup), ' (if unstable, input NONE?');
                    new_k_val = input(input_str);
                    
                    % Check for unstable data fractions
                    if ischar(new_k_val)
                        % Return NaN for k val, MI, and err
                        obj.arrMIcore{iSubgroup,3} = NaN;
                        obj.arrMIcore{iSubgroup,4} = NaN;
                        obj.arrMIcore{iSubgroup,5} = NaN;
                        
                        % User prompt to select reason for 'NONE' response
                        promptMessage = 'Please select your reason for inputting NONE';
                        titleBarCaption = 'Insufficient Data or Unstable Ks';
                        button = questdlg(promptMessage, titleBarCaption, 'Insufficnet Data', 'Unstable Ks', 'Unknown', 'Unknown');
                        
                        obj.arrMIcore{iSubgroup,6} = button;
                        
                    else
                        % Modify k value entry in core object
                        obj.arrMIcore{iSubgroup,3} = new_k_val;

                        % Recalculate MI and Error for k-value
                        core = obj.arrMIcore{iSubgroup,1};
                        r = core.get_mi(-1, 'k',new_k_val);

                        % Reset MI and error values in arrMIcore
                        obj.arrMIcore{iSubgroup,4} = r.mi;
                        obj.arrMIcore{iSubgroup,5} = r.err;
                    end
                end
                
            end
            obj.k_audited = 'Yes';
        end
        
    end
end
