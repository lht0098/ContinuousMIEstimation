clear all 

load('20200427_bl21lb21_05192022unitDunitG_analysis_c1c2.mat')

%%% test mi_ksg_core.m >> test_dataFrac_stab

a_cc.arrMIcore{1, 1}.data_frac_stab_mat  

a_cc.arrMIcore{1, 1}.test_dataFrac_stab(a_cc.arrMIcore{1, 1}.data_frac_stab_mat)

%%% test mi_ksg_core >> test_k_stab

% get the stable k's for the input argument 
stable_k_vals = a_cc.arrMIcore{1, 1}.k_values(find((a_cc.arrMIcore{1, 1}.opt_k{1,3})>=1))

a_cc.arrMIcore{1, 1}.test_k_stab(a_cc.arrMIcore{1, 1}.k_values, stable_k_vals, a_cc.arrMIcore{1, 1}.k_val_stab_mat)

load('20200427_bl21lb21_05192022unitDunitG_analysis_t2c1.mat')

%%% test mi_ksg_core.m >> find_k_value neighbor eval 

% get best_neighborStab column in core object
a_tc.arrMIcore{3, 1}.find_k_value

% check if opt_k has the best_neighborStab column
a_tc.arrMIcore{3, 1}.opt_k

% test if the script is working; outputs 
assert(size(a_tc.arrMIcore{3, 1}.opt_k, 2) == 6, 'not working')

% check that the correct value is being ouputted 
best_neighStab



