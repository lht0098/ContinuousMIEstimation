classdef mi_data_pressure < mi_data_behavior
    %    
    % Class defined to load and process data from Intan RHD files for
    % air pressure assuming that the air pressure is recorded in the
    % board_adc channel
    %
   properties
%        omitOutliers % boolean to indicate whether to omit cycles with outlier lengths from analysis
%        outliers
       method       %   method for determining pressure cycles
       threshold    %   threshold value if using method == threshold

   end

   methods
       function obj = mi_data_pressure(ID, varargin)
            % Required arguments: ID

            obj@mi_data_behavior(ID,varargin{:}); 

            obj.method = '';
            obj.threshold = nan;
       end
       
       %function set_data_files(obj, arrDataFiles, varargin)
       function set_data_files(obj, varargin)
           % Input args:
           % folder     specify folder where data files are saved
           % files      specify cell array of file names (concatenated to
           %            folder argin if specified)
           % ext        specify file extension to use (only one ext
           %            allowed)
           
           v = obj.verbose;

           default_folder = '';
           default_files = {};
           default_ext = '.csv';

           if nargin > 1
               p = inputParser;
               
               validate_folder = @(x) assert(ischar(x) && isfolder(x), 'Parameter folder must be a folder name as a string');
               addParameter(p, 'folder', default_folder, validate_folder);
    
               validate_files = @(x) assert(iscell(x), 'Parameter files must be cell of strings');
               addParameter(p, 'files', default_files, validate_files);
    
               validate_ext = @(x) assert(ischar(x) && (strfind(x, '.') == 1), 'Parameter ext must be a string that begins with a period');
               addParameter(p, 'ext', default_ext, validate_ext);
    
               p.parse(varargin{:});

               fnames = p.Results.files;
               fldr = p.Results.folder;
               fext = p.Results.ext;
           else
               fnames = default_files;
               fldr = default_fldr;
               fext = default_ext;
                
           end

           % if no folder and no files specified, user needs to choose folder for data
           if length(fnames) == 0 && length(fldr) == 0
               if v>0; disp([newline 'Please select the folder where your data are saved or select one or more data files...']); end
               
               strOpts = {'Select file...', 'Select folder...'};
               s = listdlg('PromptString', 'Select a file:', 'SelectionMode','single','ListString',strOpts);
               if s==1
                   [fnames,fldr] = uigetfile('*', 'MultiSelect', 'on');

                   [filepath,name,ext] = fileparts(fnames{1});
                   if ~strcmp(fext, ext)
                           disp([newline 'File extension does not match specified/default data file extension']);
                           disp(['Specified/default data file extension: ' fext]);
                           disp(['Data file extension of selected data file: ' ext]);
                       s = input('Do you want to use the extension of the selected data file? (y,[n])', 's');
                       if strcmp(s, 'y') || strcmp(s,'Y'); fext = ext; end
                   end

                   for i=1:length(fnames)
                       [filepath,name,ext] = fileparts(fnames{i});
                       if ~strcmp(fext, ext)
                           error('File extenions does not match specified data file extension!');
                       end
                   end
               else
                   fldr = uigetdir('*');
               end
           end
           
           % if no files are specified, but a folder is specified, confirm
           % number of data files to process
           if length(fnames) == 0 && length(fldr) > 0
               fs = dir([fldr '/*' fext]);
               if v>0
                   disp([newline 'Searching in folder: ' fldr]);
                   disp(['File extension: ' fext]);
                   disp(['Data files: ' num2str(length(fs))]);
               end
               if length(fs) > 0; fnames = {fs.name}; end
           end

           % if files are specified, check that files exist with folder and
           % all files have same extension
           if length(fnames) > 0
               if v>0
                   disp([newline 'Data folder: ' fldr]);
                   disp(['# data files: ' num2str(length(fnames))]);
               end
               for i=1:length(fnames)
                    if isfile(fullfile(fldr, fnames{i}))
                        if v>0; disp(['File: ' fnames{i}]); end
                    else
                        error(['File does not exist: ' fullfile(fldr, fnames{i})]);
                    end

                    if ~logical(strfind(fnames{i}, fext))
                        error(['File does not match specified data file extension: ' fnames{i}]);
                    end
               end
           else
               error('No files found in folder that match data file extension!');
           end


            obj.arrFiles = fnames;
            obj.strFldr = fldr;
            obj.strExt = fext;

            if v>0; disp('COMPLETE: Data files set!'); end
       end
       
%        function r = outlierScrub(obj)
%            
%            % Find cycle intervals
%            cycleTimes = obj.data.cycleTimes.data;
%            cycleIntervals = diff(cycleTimes,1,2);
%            
%            % Find mean and stdev of cycle intervals
%            m = mean(cycleIntervals);
%            sd = std(cycleIntervals);
%            
%            % Find upper and lower bound for outliers using 3*stdev rule
%            upBound = m + 3*sd;
%            lowBound = m - 3*sd;
%            
%            % Get cycle indices for non-outliers
%            cycleIdx_rmOutliers = find(cycleIntervals >= lowBound & cycleIntervals <= upBound);
%            obj.outliers = find(cycleIntervals < lowBound | cycleIntervals > upBound);
%            
%            % Get cycle times without outliers
%            cycleTimes_rmOutliers = cycleTimes(cycleIdx_rmOutliers, :);
%            r = cycleTimes_rmOutliers;
%        end
       
        function build_behavior(obj, varargin)
            % INPUTS
            % cycle_times       array of cycle onsets and offsets
            % method            hilbert - use hilbert transform for cycles
            %                   threshold - use threshold crossing
            % threshold         if method == threshold, specify value
            %
            % Process behavior pulls data to build waveform matrix
            v = obj.verbose;
            
            default_cycleTimes = [];
            default_method = 'hilbert';
            default_threshold = 0;

            if nargin > 1
                p = inputParser;

                validate_cycle_times = @(x) assert(ismatrix(x), 'Parameter cycleTimes must be an array.');
                addParameter(p, 'cycleTimes', default_ext, validate_ext);
    
                validate_method = @(x) assert(ischar(x) && (strcmp(x,'hilbert') | strcmp(x,'threshold')), 'Parameter method must be either hilbert or threshold.');
                addParameter(p, 'mtehod', default_method, validate_method);

                validate_threshold = @(x) assert(isnumeric(x), 'Parameter threshold must be a number.');
                addParameter(p, 'threshold', default_threshold, validate_threshold);

                p.parse(varargin{:});

                cycleTimes = p.Results.cycleTimes;
                method = p.Results.method;
                threshold = p.Results.threshold;
            else
                if any(strcmp(fields(obj.data), 'cycleTimes'))
                    if length(obj.data.cycleTimes) < 1
                        cycleTimes = default_cycleTimes;
                    else
                        cycleTimes = obj.data.cycleTimes;
                    end
                else
                    cycleTimes = default_cycleTimes;
                end

                if length(obj.method) < 1
                    method = default_method;
                else
                    method = obj.method;
                end

                if isnan(obj.threshold)
                    threshold = default_threshold;
                else
                    threshold = obj.threshold;
                end
            end
            
            % check for consistency of parameters
            assert(size(cycleTimes,2) ~= 2, 'Parameter cycleTimes must have only 2 columns: [cycle onset, cycle offset]');
            assert(strcmp(method,'hilbert') || strcmp(method,'threshold'), 'Parameter method must be either hilbert or threshold');
            if strcmp(method,'threshold'); assert(~isnan(threshold), 'Parameter threshold must be numeric.'); end

            
            if v>1; disp([newline '--> Building behavioral data...']); end
            

%             % Make a cell array to hold cycle data
%             nCycles = size(cycle_times,1);
%             behaviorCycles = cell(nCycles,1);
            
            
%             behavOffset = 0;
            % Iterate and load data file info
%             for i=1:length(obj.arrFiles)
            dat_pressure = [];
            dat_ts = [];

            for i=1:length(obj.arrFiles)
                if v > 1
                    disp('===== ===== ===== ===== =====');
                    disp(['Processing file ' num2str(i) ' of ' num2str(length(obj.arrFiles))]);
                    disp(['File: ' obj.strFldr '\' obj.arrFiles{i}]);
                end
                [pressure_ts, pressure_wav] = read_Intan_RHD2000_nongui_adc(fullfile(obj.strFldr, obj.arrFiles{i}), v);

                
                if v>2
                    disp(['Start time: ' num2str(pressure_ts(1)*1000.) ' ms']); 
                    disp(['End time: ' num2str(pressure_ts(end)*1000) ' ms']); 
                end
                
                dat_pressure(end+1:end+length(pressure_ts)) = pressure_wav;
                dat_ts(end+1:end+length(pressure_ts)) = pressure_ts;

                % Filter Pressure Waves
%                 filterData = obj.filterBehavior(pressure_wav, obj.Fs, filterFreq); % This will change once we update filterBehavior func
%                 if v>3; disp('--> --> Filtering data...');end
%                 filterData = pressure_wav;

                
                % Find cycles that occur within the limits of the current
                % data file
%                 cycleIxs = cycle_times(:,1) >= pressure_ts(1)*1000. & cycle_times(:,2) <= pressure_ts(end)*1000.;
%                 validCycles = cycle_times(cycleIxs,:);
%                 tStart = pressure_ts(1)*1000.; % beginning of data file in ms
                
%                 if v>3; disp(['# Cycles: ' num2str(sum(cycleIxs))]); end
                
                % Assign pressure waves to cell array
                % Consider alternative ways to save speed here?
%                 for iCycle = 1:size(validCycles,1)
%                    % Find cycle start index
%                    cycleStart = ceil((validCycles(iCycle,1)-tStart)*obj.Fs/1000.);
%                    % Find cycle end index
%                    cycleEnd = floor((validCycles(iCycle,2)-tStart)*obj.Fs/1000.);
%                    behaviorCycles{iCycle+behavOffset} = filterData(cycleStart:cycleEnd);
%                 end
                
%                 behavOffset = behavOffset + size(validCycles,1);
                
                % Check that number of stored cycles matches the offset
                % being added to cycle indices. 
%                 if ~any(behavOffset == size(find(~cellfun('isempty',behaviorCycles)))); keyboard; error('Number of stored cycles does not match the offset cycle index'); end
                
            end
            
%             obj.rawBehav = behaviorCycles;
            % Check that number of stored cycles matches the final value of the offset. 
%             if ~any(behavOffset == size(find(~cellfun('isempty',obj.rawBehav)))); error('Total number of stored cycles does not match the final offset of cycle index'); end
            


            if v>0; disp('COMPLETE: Behavioral data loaded!'); end
       end
       
       
       function [filterData] = filterBehavior(obj, behavior, cycleFreq, filterFreq)
           
            if v>1; disp([newline '--> Filtering behavioral data...']); end 
           
            % This function prepares the raw behavioral data for analysis
            % Convert cycle freq to length of gaussian in samples
            cycleLengthSeconds = 1/cycleFreq;
            cycleLengthSamples = cycleLengthSeconds * obj.bFs;
            % Convert filter freq to width of gaussian in samples
            filterWinSeconds = 1/filterFreq;
            filterWinSamples = filterWinSeconds * obj.bFs;

            % Find alpha value for input to gaussian window function.
            alpha = (cycleLengthSamples - 1)/(2*filterWinSamples);

            % Generate the gaussian window for filter
            g = gausswin(cycleLengthSamples, alpha);

            filterData = conv(behavior,g,'same');

            if v>0; disp('COMPLETE: Behavioral data filetered!'); end
       end        
        
       function r = get_behavior(obj, timeBase, feature, start, dur, nSamp, varargin)
            % Implemented to return different features of behavior cycles after processing raw waveform matrix data
            % i.e., raw data, PCA, residual, area

            v = obj.verbose;
            
            p = inputParser;
            
            % Required: timeBase
            % 'time' or 'phase' to indicate how to calculate time
            valid_timeBase = {'time' 'phase'};
            validate_timeBase = @(x) assert(ischar(x) && ismember(x, valid_timeBase), 'timeBase must be a char/string: time, phase');
            p.addRequired('timeBase', validate_timeBase);
            
            % Required: feature
            % 'raw', 'pca', 'residual'
            valid_feature = {'raw' 'pca' 'residual'};
            validate_feature = @(x) assert(ischar(x) && ismember(x, valid_feature), 'feature must be a char/string: raw, pca, residual');
            p.addRequired('feature', validate_feature);
            
            % Required: start
            % in ms or rad/phase based on timeBase
            validate_start = @(x) assert(isnumeric(x), 'start must be numeric corresponding to timeBase');
            p.addRequired('start', validate_start);
            
            % Required: dur
            % in ms or rad/phase based on timeBase
            validate_dur = @(x) assert(isnumeric(x), 'dur must be numeric indicating the duration of the vector corresponding to timeBase');
            p.addRequired('dur', validate_dur);
            
            % Required: nSamp
            % in samples/data points
            validate_nSamp = @(x) assert(isnumeric(x) && mod(x,1) == 0, 'nSamp must be an integer indicating number of data points in sample');
            p.addRequired('nSamp', validate_nSamp);
            
            % Optional: nPC
            % integer to indicate number of principal components to return
            default_nPC = 3;
            validate_nPC = @(x) assert(isnumeric(x) && mod(x,1) == 0, 'nPC must be an integer indicating number of principal components to return');
            p.addOptional('nPC', default_nPC, validate_nPC);
            
            p.parse(timeBase, feature, start, dur, nSamp, varargin{:});
            
            timeBase = p.Results.timeBase;
            feature = p.Results.feature;
            start = p.Results.start;
            dur = p.Results.dur;
            nSamp = p.Results.nSamp;
            nPC = p.Results.nPC;
            
            
            if v>1; disp([newline '--> Getting formatted behavioral data...']); end
            
            switch(timeBase)
                case('time')
                    if v>2; disp([newline '--> --> Formatting behavior: time']); end
                    % Find the lengths of the cycles in samples
                    cycleLengths_samples = cell2mat(cellfun(@(x) length(x), obj.rawBehav, 'UniformOutput', false));

                    % Find the sample associated with the start time for each
                    % cycle in ms
                    start_samples = ceil(start*obj.Fs/1000.);

                    % Find the number of samples that encompases the window of
                    % interest for the cycles. 
                    % Convert windowOFInterest from ms to seconds
                    dur_samples = ceil(dur*obj.Fs/1000.);

                    % Find the stop sample of the window of interest for each
                    % cycle
                    stop_samples = start_samples + dur_samples;

                    if v>3
                        disp(['Cycle Start Time: ' num2str(start) ' ms']);
                        disp(['Cycle Sample Duration: ' num2str(dur) ' ms']);
                        disp(['Sampled Data Points: ' num2str(nSamp)]);
                    end
                    
                    % Set up empty matrix to store pressure data.
                    nCycles = length(cycleLengths_samples);
                    cycle_behavior = zeros(nCycles, nSamp);

                    tic
                    for cycle_ix = 1:nCycles
                       % Document all of the data points for the window of
                       % interest
                       dat = obj.rawBehav{cycle_ix};
                       if length(dat) >= stop_samples
                           cycle_data = dat(start_samples:stop_samples);
                           % Resample to get only the desired number of points
                           newSamples = round(linspace(1,length(cycle_data),nSamp));

                           resampled_cycle_data = cycle_data(newSamples);
                           cycle_behavior(cycle_ix,:) = resampled_cycle_data;   
                       else
                           cycle_behavior(cycle_ix,:) = nan(1,nSamp);
                       end
                    end
                    toc
                    
                    if v>2; warning(['WARNING: Ommitted ' num2str(sum(isnan(cycle_behavior(:,1)))) ' empty cycles']); end

                case('phase')
                    if v>2; disp([newline '--> --> Formatting behavior: phase']); end
                    % Find the lengths of the cycles in samples
                    cycleLengths_samples = cell2mat(cellfun(@(x) length(x), obj.rawBehav, 'UniformOutput', false));

                    % Find the sample associated with the start phase for each
                    % cycle
                    start_samples = ceil(cycleLengths_samples.*(start/(2*pi)));

                    % Find the number of samples that encompases the window of
                    % interest for the cycles.
                    windowOfInterest_samples = ceil((cycleLengths_samples).*(dur./(2*pi)));

                    % Find the stop sample of the window of interest for each
                    % cycle
                    stop_samples = start_samples + windowOfInterest_samples;

                    if v>3
                        disp(['Cycle Start Phase: ' num2str(start) ' rad']);
                        disp(['Cycle Sample Duration: ' num2str(dur) ' rad']);
                        disp(['Sampled Data Points: ' num2str(nSamp)]);
                    end
                    
                    % Set up empty matrix to store pressure data.
                    nCycles = length(cycleLengths_samples);
                    cycle_behavior = zeros(nCycles, nSamp);

                    tic
                    for cycle_ix = 1:nCycles
                        % Document all of the data points for the window of
                        % interest
                        if stop_samples(cycle_ix) > start_samples(cycle_ix)
                            dat = obj.rawBehav{cycle_ix};
                            cycle_data = dat(start_samples(cycle_ix):stop_samples(cycle_ix));

                            newSamples = round(linspace(1,length(cycle_data),nSamp));
                            % Resample to get only the desired number of points
                            resampled_cycle_data = cycle_data(newSamples);
                            cycle_behavior(cycle_ix, 1:nSamp) = resampled_cycle_data;   
                        else
                            cycle_behavior(cycle_ix,:) = nan(1,nSamp);
                        
                        end
                    end              
                    toc
                    
                    if v>2; warning(['WARNING: Ommitted ' num2str(sum(isnan(cycle_behavior(:,1)))) ' empty cycles']); end
            end
            
            switch(feature)
                case('raw')
                    if v>2; disp([newline '--> --> Processing behavior: raw']); end
                    r = cycle_behavior;
                case('pca')
                    if v>2; disp([newline '--> --> Processing behavior: PCA (' num2str(nPC) ') PCs']); end
                    [~,score,~] = pca(cycle_behavior);
                    r = score(:,1:nPC);
                case('residual')
                    if v>2; disp([newline '--> --> Processing behavior: residual']); end
                    r = cycle_behavior - mean(cycle_behavior,1,'omitnan');     
            end
            
            if v>0; disp('COMPLETE: Behavioral data retrieved!'); end
        end       
   end    
end

%            obj.b_startPhase = {};