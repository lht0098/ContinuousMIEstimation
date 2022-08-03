
spk_ct_tempD2 = NaN(1, size(temp_D2, 2));

spk_ct_tempD3 = NaN(1, size(temp_D3, 2));

for i = 1:size(spk_ct_tempD2, 2)

    spk_ct_tempD2(1, i) = sum((temp_D2 > temp_cyc(i,1)) & (temp_D2 < temp_cyc(i,2))); 

    spk_ct_tempD3(1, i) = sum((temp_D3 > temp_cyc(i,1)) & (temp_D3 < temp_cyc(i,2))); 

end 