classdef mi_data_behavior < mi_data
    properties
        rawBehav % struct of raw behavior data and information
        
        arrFiles % list of files with raw behavioral data
        strFldr % path to where the data files are saved
        strExt % extension of data file
        
        cycle_select
    end
    
    methods
        function obj = mi_data_behavior(ID, varargin)
            % Required arguments: ID
            obj@mi_data(ID,varargin{:});
        end
        
        function add_cycleTimes(obj, data, dataInfo, Fs, varargin)
            add_data(obj, data, dataInfo, Fs, 'cycleTimes', varargin{:});
        end
        
        function build_behav(obj)
            % Pull waveform data from data files according to
            % data/cycleTimes
            warning('NOT IMPLEMENTED ERROR');
        end
        
        function r = get_behav(obj)
            % Return behavior data as matrix of raw waveform/PCA/ICA/etc.
            warning('NOT IMPLEMENTED ERROR');
        end
        
        function r = get_cycleTimes(obj, varargin)
            % Return behavior cycle times as matrix of on/off pairs
            p = inputParser;

            default_cycle_select = -1;
            validate_cycle_select = @(x) assert(isnumeric(x), 'cycle_selection must be numeric array');
            p.addOptional('cycle_select', default_cycle_select, validate_cycle_select)
            p.KeepUnmatched = 1;
            p.parse(varargin{:});
            
            obj.cycle_select = p.Results.cycle_select;

            all_cycles = get_data(obj, 'cycleTimes');
            
            if obj.cycle_select == -1
                r = all_cycles;
            else
                r = all_cycles(obj.cycle_select,:);
            end
        end
        
    end
end