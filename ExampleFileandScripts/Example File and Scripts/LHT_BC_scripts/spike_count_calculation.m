% get the spike count of the unit 
spk_ct = a_cc.arrMIcore{1, 1}.x;

% get the average spike count of the unit
disp(['average JITTERED spike count = ' num2str(mean(spk_ct))])

% get the total number of spike counts 
disp(['total jittered spike count = ' num2str(sum(spk_ct))])

% get the cycle times 
cycl_tm = a_cc.objBehav.data.cycleTimes.data; 

% get the beginning of the first cycle and the end of the last cycle
cycl_lims = [cycl_tm(1, 1) cycl_tm(end, 2)];

% get the spike times of unit 
spk_tm = a_cc.objData.data.unitB.data;

% get the number of spikes that is within the beginning end of the original
% total cycles 
valid_spk = sum(spk_tm >= cycl_lims(1) & spk_tm <= cycl_lims(2));

% get the average spike count of each cycle 
disp([newline 'EXPECTED average spike count = ' num2str(valid_spk/size(cycl_tm, 1))])

% get the total number of spikes for the unit  
disp(['total expected spike count = ' num2str(size(spk_tm, 2))])

% get the total number of spikes that were dropped as a result of the
% jitter
disp([newline 'dropped spikes = ' num2str(size(spk_tm, 2) - sum(spk_ct)) ' (' num2str( ((size(spk_tm, 2) - sum(spk_ct))/size(spk_tm, 2)) * 100 ) '%)'])

