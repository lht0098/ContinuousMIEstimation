%% Count-Count 
cc_units = [1, 2, 3, 4, 5, 6];
mi = [0, 0.0932, 0.3482, 0.0238, 0.1002, 0.0553];
std = [0, 0.0106, 0.0153, 0.0192, 0.0131, 0.0111];

cc_fig = figure()

bar(cc_units, mi) 

hold on

errorbar(cc_units, mi, std, 'bp')
xticklabels({'DE', 'GE', 'BD', 'BE', 'BG', 'DG'})
ylim([0 0.4])

hold off

%% Timing-Count 
tc_units = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
mi = [0.0063, 9.5919e-04, 0.0059, 0.0206, 0.0971, 0.0345, 0.0223, 0.0830, 0.0742, 0.0808, 0.0602, 0.0729];
std = [0.0166, 0.0430, 0.0133, 0.0475, 0.0342, 0.0815, 0.0146, 0.1241, 0.0328, 0.0564, 0.0518, 0.0476];

tc_fig = figure()

bar(tc_units, mi) 

hold on

errorbar(tc_units, mi, std, 'bp')
xticklabels({'DE', 'ED', 'GE', 'EG', 'BD', 'DB', 'BE', 'EB', 'BG', 'GB', 'GD', 'DG'})
ylim([0 0.25])

hold off

%% Timing-Timing 
tt_units = [1, 2, 3, 4];
mi = [0, 0.0345, 0.0320, 0.0477];
std = [0, 0.0815, 0.0296, 0.0705];

tt_fig = figure()

bar(tt_units, mi) 

hold on

errorbar(tt_units, mi, std, 'bp')
xticklabels({'GE', 'BD', 'BG', 'DG'})
ylim([0 0.25])

hold off

%% Count-Count 
cc_units = [1, 2, 3, 4, 5, 6];
mi = [0, nan, 0.3482, nan, 0.0553, nan];
std = [0, nan, 0.0153, nan, 0.0111, nan];

cc_fig = figure()

bar(cc_units, mi) 

hold on

errorbar(cc_units, mi, std, 'bp')
title('Count Count: Mu = 0, Sigma = 0.02')
xticklabels({'DE', 'DE*', 'BD', 'BD*', 'DG', 'DG*'})
ylim([0 0.4])

hold off

%% Timing-Count 
tc_units = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
mi = [0.0063, 0.0072, 9.5919e-04, 0, 0.0971, 0.1043, 0.0345, 0.1397, 0.0602, 0.0697, 0.0729, 0.0622];
std = [0.0166, 0.0142, 0.0430, 0.0137, 0.0342, 0.0399, 0.0815, 0.1069, 0.0518, 0.0503, 0.0476, 0.0149];

tc_fig = figure()

bar(tc_units, mi) 

hold on

errorbar(tc_units, mi, std, 'bp')
title('Timing Count: Mu = 0, Sigma = 1e-04')
xticklabels({'DE', 'DE*', 'ED', 'ED*', 'BD', 'BD*', 'DB', 'DB*', 'GD', 'GD*', 'DG', 'DG*'})
ylim([0 0.25])

hold off

%% Timing-Count 
tc_units = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
mi = [0.0063, nan, 9.5919e-04, nan, 0.0971, 0.1043, 0.0345, nan, 0.0602, nan, 0.0729, nan];
std = [0.0166, nan, 0.0430, nan, 0.0342, 0.0399, 0.0815, nan, 0.0518, nan, 0.0476, nan];

tc_fig = figure()

bar(tc_units, mi) 

hold on

errorbar(tc_units, mi, std, 'bp')

xticklabels({'DE', 'DE', 'ED', nan, 'BD', 'BD*', 'DB', nan, 'GD', nan, 'DG', nan})
ylim([0 0.25])

hold off
