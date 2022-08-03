clear all

load 20200427_bl21lb21_04052022unitDunitG_analysis_c1c2.mat 

core_object = a_cc.arrMIcore{1, 1};

% a = 7
% 
% x = core_object.data_frac_stab_mat;
% 
% z = core_object.k_val_stab_mat; 
% 
% x_matrix = x(:,:,a);
% 
% x_vector = reshape(x_matrix(:,:,1),1,[]); % change the matrix into an array
% 
% edges = 0:0.1:round(max(x_vector)); 
% 
% figure(a)
% 
% hist = histogram(x_vector, edges)
% 
% hist_counts = hist.BinCounts;
% 
% title('Count of Values in Data Fraction Stability Matrix')
% xlabel('Matrix Value')
% ylabel('Count')
% 
% xlim([0, 3])
% 
% xticks(0:0.1:round(max(x_vector)))
% 
% yticks(0:5:max(hist_counts))
% 
% figure_name = ['a_cc_', 'data_frac_stab_mat_', num2str(a), '.png']
% 
% saveas(figure(a), figure_name)


a = 1

z = core_object.k_val_stab_mat; 

z_matrix = z(:,:,a);

z_vector = reshape(z_matrix(:,:,1),1,[]); % change the matrix into an array

edges = 0:0.1:round(max(z_vector));

figure(a)

hist = histogram(z_vector, edges)

hist_counts = hist.BinCounts;

title('Count of Values in K Stability Matrix')
xlabel('Matrix Value')
ylabel('Count')

xlim([0, 3])

xticks(0:0.1:round(max(z_vector)))

yticks(0:5:max(hist_counts))

figure_name = ['a_cc_', 'k_stab_mat_', num2str(a), '.png']

saveas(figure(a), figure_name)


