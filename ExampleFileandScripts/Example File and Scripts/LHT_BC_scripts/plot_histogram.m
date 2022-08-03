clear all

load 20200427_bl21lb21_04052022unitDunitG_analysis_c1c2.mat 

core_object = a_cc.arrMIcore{1, 1};

 x = core_object.data_frac_stab_mat;

 for i = 1:size(x,3)

     x_matrix = x(:,:,i);

     x_vector = reshape(x_matrix(:,:,1),1,[]); % change the matrix into an array
    
     edges = 0:0.5:round(max(x_vector));
    
     figure(i)

     hist = histogram(x_vector, edges)
    
     hist_counts = hist.BinCounts;

     title('Count of Values in Data Fraction Stability Matrix')
     xlabel('Matrix Value')
     ylabel('Count')

     xlim([0, round(max(x_vector))])
    
     xticks(0:1:round(max(x_vector)))
    
     yticks(0:5:max(hist_counts))

     figure_name = ['a_cc', 'data_frac_stab_mat_', num2str(i), '.png']

     saveas(figure(i), figure_name)

 end
 

 