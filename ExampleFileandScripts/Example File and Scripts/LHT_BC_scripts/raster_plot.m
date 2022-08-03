% Raster Plot for Auditing

clear all

load('20200427_bl21lb21_06032020data.mat')

% Set the unit of interest 
unit = 'unitB';

% Get the cycle times
cycleTimes = b.data.cycleTimes.data;

% Get the spike times
spikeTimes = d.data.(unit).data;

% Create an empty array of NaNs 
timeDiff = NaN(size(cycleTimes,1), 1);

% Find the length of time for each cycle 
for i = 1:size(timeDiff, 1)
    timeDiff(i, 1) = cycleTimes(i, 2) - cycleTimes(i, 1);
end 

% Initiate an array to hold the number of max spikes per cycle
cycleMax = NaN(size(cycleTimes, 1), 1);

for i = 1:size(cycleMax, 1)
    
    % Set the range of interest 
    cycleStart = cycleTimes(i, 1);
    cycleEnd = cycleTimes(i, 2);

    % Find spike times that are less than the start and
    % the end of the cycle
    idx = spikeTimes >= cycleStart & spikeTimes <= cycleEnd;
    
    % Find the column indexes of the spike times in the cycle of interest
    column_idx = find(idx);
    spikes_cycle{i,1} = spikeTimes(column_idx);
    
    % Get the number of max spikes per cycle 
    cycleMax(i) = size(spikes_cycle{i}, 2);

end

% Find the greatest number of spikes a cycle can have
maxSpikes = max(cycleMax);

% Create an empty array of NaNs
spikes_per_cycle = NaN(size(cycleTimes, 1), maxSpikes);

total_spikes = NaN(size(cycleTimes, 1), 1);

for i = 1:size(spikes_per_cycle, 1)

    % Get the number of spikes in the cycle
    num_spikes = size(spikes_cycle{i}, 2);

    total_spikes(i, 1) = num_spikes;

    % Replace the NaNs with spikes times
    spikes_per_cycle(i, 1:num_spikes) = cell2mat(spikes_cycle(i,:));

end

% Create an array to hold the offset spike times
offset_spikes = NaN(size(spikes_per_cycle, 1), size(spikes_per_cycle, 2));

for i = 1:size(offset_spikes, 1)
    offset_spikes(i,:) = spikes_per_cycle(i,:) - cycleTimes(i, 1);
end

figure()

for i = 1:size(offset_spikes, 1)
    
    % Find the times of the spikes
    spike_idx = find(~isnan(offset_spikes(i,:)));
    
    % Create a circle for every spike time in the cycyle
    plot([offset_spikes(i, spike_idx)], i.* ones(1, length(spike_idx)), linestyle = 'none', marker = '.', color = 'm')

    % Keep each line from each trial
    hold on 

end

ylim([0 8000])
xlim([0 0.6])

% Reverse the y axis to make a closer replica of Rachel's plots
set(gca, 'YDir','reverse')