%% a_cc.arrMIcore{1, 1}

clear all 

load('20200427_bl21lb21_05192022unitDunitG_analysis_c1c2.mat')

%%% test mi_ksg_core.m >> test_dataFrac_stab

a_cc.arrMIcore{1, 1}.data_frac_stab_mat  

a_cc.arrMIcore{1, 1}.test_dataFrac_stab(a_cc.arrMIcore{1, 1}.data_frac_stab_mat)

%%% test mi_ksg_core >> test_k_stab

% get the stable k's for the input argument 
stable_k_vals = a_cc.arrMIcore{1, 1}.k_values(find((a_cc.arrMIcore{1, 1}.opt_k{1,3})>=1))

a_cc.arrMIcore{1, 1}.test_k_stab(a_cc.arrMIcore{1, 1}.k_values, stable_k_vals, a_cc.arrMIcore{1, 1}.k_val_stab_mat)

%%% test mi_ksg_core.m >> find_k_value neighbor eval 

% get best_neighborStab column in core object
a_cc.arrMIcore{1, 1}.find_k_value

% check if opt_k has the best_neighborStab column
a_cc.arrMIcore{1, 1}.opt_k

% test if the script is working; outputs 
assert(size(a_cc.arrMIcore{1, 1}.opt_k, 2) == 5, 'not working')

% check that the correct value is being ouputted; it should be 0
a_cc.arrMIcore{1, 1}.opt_k{1,5}

%%% test mi_ksg_core.m >> get_opk_mi

a_cc.arrMIcore{1, 1}.get_opk_mi

% check that the correct k value is being outputted 
a_cc.arrMIcore{1, 3}

%% 

clear all

load('20200427_bl21lb21_05192022unitDunitG_analysis_t1c2.mat')

%%% test mi_ksg_core.m >> find_k_value neighbor eval 

% get best_neighborStab column in core object
a_tc.arrMIcore{1, 1}.find_k_value

% check if opt_k has the best_neighborStab column
a_tc.arrMIcore{1, 1}.opt_k

% test if the script is working; outputs 
assert(size(a_tc.arrMIcore{1, 1}.opt_k, 2) == 5, 'not working')

% check that the correct value is being ouputted 
a_tc.arrMIcore{1, 1}.opt_k{1,5}

%%% test mi_ksg_core.m >> get_opk_mi

a_tc.arrMIcore{1, 1}.get_opk_mi

% check that the correct k value is being outputted 
a_tc.arrMIcore{1, 3}

%%

clear all

load('20200427_bl21lb21_05192022unitDunitG_analysis_t1c2.mat')

%%% test mi_ksg_core.m >> find_k_value neighbor eval 

% get best_neighborStab column in core object
a_tc.arrMIcore{4, 1}.find_k_value

% check if opt_k has the best_neighborStab column
a_tc.arrMIcore{4, 1}.opt_k

% test if the script is working; outputs 
assert(size(a_tc.arrMIcore{4, 1}.opt_k, 2) == 5, 'not working')

% check that the correct value is being ouputted 
a_tc.arrMIcore{4, 1}.opt_k{1,5}

%%% test mi_ksg_core.m >> get_opk_mi

a_tc.arrMIcore{4, 1}.get_opk_mi

% check that the correct k value is being outputted 
a_tc.arrMIcore{4, 3}

%% a_tc.arrMIcore{3, 1}

clear all

load('20200427_bl21lb21_05192022unitDunitG_analysis_t2c1.mat')

%%% test mi_ksg_core.m >> find_k_value neighbor eval 

% get best_neighborStab column in core object
a_tc.arrMIcore{3, 1}.find_k_value

% check if opt_k has the best_neighborStab column
a_tc.arrMIcore{3, 1}.opt_k

% test if the script is working; outputs 
assert(size(a_tc.arrMIcore{3, 1}.opt_k, 2) == 5, 'not working')

% check that the correct value is being ouputted 
a_tc.arrMIcore{3, 1}.opt_k{1,5}

%%% test mi_ksg_core.m >> get_opk_mi

a_tc.arrMIcore{3, 1}.get_opk_mi

% check that the correct k value is being outputted 
a_tc.arrMIcore{3, 3}

%% a_tt.arrMIcore{16, 1}

clear all

load('20200427_bl21lb21_05192022unitDunitG_analysis_t1t2.mat')

%%% test mi_ksg_core.m >> find_k_value neighbor eval 

% get best_neighborStab column in core object
a_tt.arrMIcore{16, 1}.find_k_value

% check if opt_k has the best_neighborStab column
a_tt.arrMIcore{16, 1}.opt_k

% test if the script is working; outputs 
assert(size(a_tt.arrMIcore{16, 1}.opt_k, 2) == 5, 'not working')

% check that the correct value is being ouputted 
a_tt.arrMIcore{16, 1}.opt_k{1,5}

%%% test mi_ksg_core.m >> get_opk_mi

a_tt.arrMIcore{16, 1}.get_opk_mi

% check that the correct k value is being outputted 
a_tt.arrMIcore{16, 3}

%% %% a_tt.arrMIcore{35, 1}

clear all

load('20200427_bl21lb21_05192022unitDunitG_analysis_t1t2.mat')

%%% test mi_ksg_core.m >> find_k_value neighbor eval 

% get best_neighborStab column in core object
a_tt.arrMIcore{35, 1}.find_k_value

% check if opt_k has the best_neighborStab column
a_tt.arrMIcore{35, 1}.opt_k

% test if the script is working; outputs 
assert(size(a_tt.arrMIcore{35, 1}.opt_k, 2) == 5, 'not working')

% check that the correct value is being ouputted 
a_tt.arrMIcore{35, 1}.opt_k{1,5}

%%% test mi_ksg_core.m >> get_opk_mi

a_tt.arrMIcore{35, 1}.get_opk_mi

% check that the correct k value is being outputted 
a_tt.arrMIcore{35, 3}