% Raster Plot for Auditing

clear all

load('20200427_bl21lb21_06032020data.mat')

% Set the unit of interest 
unit = 'unitD'

% Get the cycle times
cycleTimes = b.data.cycleTimes.data;

% Get the spike times
spikeTimes = d.data.(unit).data;

load('20200427_bl21lb21_06032020data_unitD_mu0_sig00001.mat')

corrupt_spikeTimes = d.data.(unit).data;

% Create an empty array of NaNs 
timeDiff = NaN(size(cycleTimes,1), 1);

% Find the length of time for each cycle 
for i = 1:size(timeDiff, 1)
    timeDiff(i, 1) = cycleTimes(i, 2) - cycleTimes(i, 1);
end 

% Create an empty array of NaNs
spikes_per_cycle = NaN(size(cycleTimes, 1), 8);
corrupt_spikes_per_cycle = NaN(size(cycleTimes, 1), 8);

for i = 1:size(spikes_per_cycle, 1)
    
    % Set the range of interest 
    cycleStart = cycleTimes(i, 1);
    cycleEnd = cycleTimes(i, 2);

    idx = spikeTimes >= cycleStart & spikeTimes <= cycleEnd;
    corrupt_idx = corrupt_spikeTimes >= cycleStart & corrupt_spikeTimes <= cycleEnd;
    
    % Find the column indexes of the spike times in the cycle of interest
    column_idx = find(idx);
    corrupt_column_idx = find(corrupt_idx)

    spikes_cycle = spikeTimes(column_idx);
    corrupt_cycle = corrupt_spikeTimes(corrupt_column_idx);

    column_diff = size(spikes_per_cycle, 2) - size(column_idx, 2);
    corrupt_diff = size(corrupt_spikes_per_cycle, 2) - size(corrupt_column_idx, 2);

    spikes_per_cycle(i,:) = horzcat(spikes_cycle, NaN(1, column_diff));
    corrupt_spikes_per_cycle(i,:) = horzcat(corrupt_cycle, NaN(1, corrupt_diff));
end

offset_spikes = NaN(size(cycleTimes, 1), 8);
corrupt_offset_spikes = NaN(size(cycleTimes, 1), 8);

for i = 1:size(offset_spikes, 1)
    offset_spikes(i,:) = spikes_per_cycle(i,:) - cycleTimes(i, 1);
    corrupt_offset_spikes(i,:) = corrupt_spikes_per_cycle(i,:) - cycleTimes(i, 1);
end

figure()

for i = 1:size(offset_spikes, 1)
    
    % Find the times of the spikes
    spike_idx = find(~isnan(offset_spikes(i,:)));
    corrupt_spike_idx = find(~isnan(corrupt_offset_spikes(i,:)));
    
    % Create a circle for every spike time in the cycyle
    plot([offset_spikes(i, spike_idx)], 2*i.* ones(1, length(spike_idx)), linestyle = 'none', marker = '.', color = 'b')

    % Keep each line from each trial
    hold on 

    % Create a circle for every spike time in the cycyle - corrupted 
    plot([corrupt_offset_spikes(i, corrupt_spike_idx)], (2*i-1).* ones(1, length(corrupt_spike_idx)), linestyle = 'none', marker = '.', color = 'g')

    hold on

end

ylim([0 10])
xlim([0 0.6])

% Reverse the y axis to make a closer replica of Rachel's plots
set(gca, 'YDir','reverse')
