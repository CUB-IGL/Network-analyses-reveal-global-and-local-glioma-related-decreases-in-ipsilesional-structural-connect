fwe_no_thr = csvread('tfnbs_fwe_1mpvalue.csv',1,0);
connectome_mask = fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask = connectome_mask > 0.95;
fwe_no_thr = fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_fwe_1mpvalue.edge',fwe_no_thr,'\t');


pre_fwe_no_thr = csvread('tfnbs_precentral_fwe_1mpvalue.csv',1,0);
connectome_mask =pre_fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
pre_fwe_no_thr= pre_fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_precentral_fwe_1mpvalue.edge',pre_fwe_no_thr,'\t');


post_fwe_no_thr = csvread('tfnbs_postcentral_fwe_1mpvalue.csv',1,0);
connectome_mask = post_fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
post_fwe_no_thr= post_fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_postcentral_fwe_1mpvalue.edge',post_fwe_no_thr,'\t');


insula_fwe_no_thr = csvread('tfnbs_insular_fwe_1mpvalue.csv',1,0);
connectome_mask = insula_fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask = connectome_mask > 0.95;
insula_fwe_no_thr= insula_fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_insula_fwe_1mpvalue.edge',insula_fwe_no_thr,'\t');


fwe_thr = csvread('tfnbs_thr_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_thr= fwe_thr .* connectome_mask;
dlmwrite('tfnbs_fwe_thr_1mpvalue.edge',fwe_thr,'\t');

fwe_pre_thr = csvread('tfnbs_precentral_thr_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_pre_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_pre_thr= fwe_pre_thr .* connectome_mask;
dlmwrite('tfnbs_precentral_thr_fwe_1mpvalue.edge',fwe_pre_thr,'\t');

fwe_no_thr = csvread('tfnbs_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_no_thr= fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_fwe_1mpvalue.edge',fwe_no_thr,'\t');

fwe_no_thr = csvread('tfnbs_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_no_thr= fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_fwe_1mpvalue.edge',fwe_no_thr,'\t');
