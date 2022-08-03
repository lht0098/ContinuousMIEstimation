% Raster Plot for Auditing

clear all

load('20200427_bl21lb21_06032020data.mat')

% Set the unit of interest 
unit1 = 'unitD'
unit2 = 'unitB'

% Get the spike times
spikeTimes1 = d.data.(unit1).data;
spikeTimes2 = d.data.(unit2).data;

load('20200427_bl21lb21_06032020data_unitD_mu0_sig_15_.mat')

spikeTimes_jitter = d.data.(unit1).data;

figure()

plot(spikeTimes1, ones(size(spikeTimes1)), linestyle = 'none', marker = '.', color = 'b')

hold on 

plot(spikeTimes2, 2*ones(size(spikeTimes2)), linestyle = 'none', marker = '.', color = 'g')

hold on 

plot(spikeTimes_jitter, 3*ones(size(spikeTimes_jitter)), linestyle = 'none', marker = '.', color = 'r')

ylim([0 4])