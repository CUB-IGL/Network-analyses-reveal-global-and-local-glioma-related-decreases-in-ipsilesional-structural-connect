clear;clc;
%read connectome matrices and seperate them into healthy and pathologic
%hemipsheres' connectomes
List = dir('S*.csv')
left_patho = [81 82 139 164 284 293 299 302 354 364 374 381 390 391 395 407]
for i = 1:length(List)
    filename = List(i).name;
    Index = sscanf(filename, 'S%d');
    control_name = ['S',num2str(Index),'_c.csv'];
    pathologic_name = ['S',num2str(Index),'_p.csv'];
    connectome = csvread(filename);
    
    %assign brainstem node (node79) to both hemisphere matrices
    connectome_left = connectome([1:39,79],[1:39,79]);  
    connectome_right = connectome(40:79,40:79);
    
    if any(left_patho==Index)
        csvwrite(pathologic_name,connectome_left);
        csvwrite(control_name,connectome_right);
    else
        csvwrite(pathologic_name,connectome_right);
        csvwrite(control_name,connectome_left);
    end
end
%% motor-sensory
% List = dir('S*.csv')
% left_patho = [15 81 82 139 164 284 293 299 302 354 364 374 381 390 391 395 407]
% for i = 1:length(List)
%     filename = List(i).name;
%     Index = sscanf(filename, 'S%d');
%     control_name = ['S',num2str(Index),'_c.csv'];
%     pathologic_name = ['S',num2str(Index),'_p.csv'];
%     connectome = csvread(filename);
%     
%     connectome_left = connectome([1:7,15],[1:7,15]);
%     connectome_right = connectome(8:15,8:15);
%     
%     if any(left_patho==Index)
%         csvwrite(pathologic_name,connectome_left);
%         csvwrite(control_name,connectome_right);
%     else
%         csvwrite(pathologic_name,connectome_right);
%         csvwrite(control_name,connectome_left);
%     end
% end

%%for tfnbs all
fwe_no_thr = csvread('tfnbs_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_no_thr= fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_fwe_1mpvalue.edge',fwe_no_thr,'\t');

%%for tfnbs insular group
fwe_no_thr = csvread('tfnbs_insular_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_no_thr= fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_insular_fwe_1mpvalue.edge',fwe_no_thr,'\t');

%%for tfnbs frontal group
fwe_no_thr = csvread('tfnbs_frontal_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_no_thr= fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_frontal_fwe_1mpvalue.edge',fwe_no_thr,'\t');

%%for tfnbs precentral group
fwe_no_thr = csvread('tfnbs_precentral_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_no_thr= fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_precentral_fwe_1mpvalue.edge',fwe_no_thr,'\t');

%%for tfnbs postcentral group
fwe_no_thr = csvread('tfnbs_postcentral_fwe_1mpvalue.csv',1,0);
connectome_mask =fwe_no_thr;
%connectome_mask_no_thr(connectome_mask_no_thr == 0) = 1;
connectome_mask=connectome_mask > 0.95;
fwe_no_thr= fwe_no_thr .* connectome_mask;
dlmwrite('tfnbs_postcentral_fwe_1mpvalue.edge',fwe_no_thr,'\t');

%%for tfnbs minimally connected components
fwe_mc = csvread('tfnbs_mcfwe_1mpvalue.csv');
connectome_mask_mc =fwe_mc;
%connectome_mask_mc(connectome_mask_mc == 0) = 1;
connectome_mask_mc=connectome_mask_mc >0.95;
fwe_mc= fwe_mc .* connectome_mask_mc;
dlmwrite('fwe_mc.edge',fwe_mc,'\t')



tvalue = csvread('tfnbs_nothrtvalue.csv');
tvalue_nozero= tvalue .* connectome_mask;
dlmwrite('tvalue_nozero.edge',tvalue_nozero,'\t');

%calculate the minimal connectivity

%input subject number into list_num
list_num = [28,31,32,81,82,110,112,130,139,149,155,164,180,208,212,284,285,286,287,290,291,292,293,299,302,354,356,364,374,381,382,390,391,395,400,401,407];
Index_all=[]
for i = 1:length(list_num)
    %set input and output names of csv files 
    Index_num = list_num(i);
    con_file_name = ['S',num2str(Index_num),'_c.csv'];
    patho_file_name = ['S',num2str(Index_num),'_p.csv'];
    Index_all = [Index_all,Index_num];
    save_name_con = ['S',num2str(Index_num),'_c_thr.csv'];
    save_name_patho = ['S',num2str(Index_num),'_p_thr.csv'];
   
    %read csv and call computing function
    connectome_c = csvread(con_file_name);
    connectome_p = csvread(patho_file_name);
    connecomte_c_thr = compute_min_conn_matrix(connectome_c);
    connecomte_p_thr = compute_min_conn_matrix(connectome_p);
   
    %save files
    dlmwrite(save_name_con,full(connecomte_c_thr));   
    dlmwrite(save_name_patho,full(connecomte_p_thr));
   
end
    
%% BCT analysis
%% Original connectome data and node degree
list_num = [28,81,82,110,130,149,155,164,208,284,285,286,287,291,293,374,382,390,391,400];
Index_all=[]
connectome_data={}
degrees={}
for i = 1:length(list_num)
    Index_num = list_num(i);
    con_file_name = ['S',num2str(Index_num),'_c.csv'];
    patho_file_name = ['S',num2str(Index_num),'_p.csv'];
    thr_con_file_name = ['S',num2str(Index_num),'_c_thr.csv'];
    thr_patho_file_name = ['S',num2str(Index_num),'_p_thr.csv'];
    
    con_name = ['S',num2str(Index_num),'_c'];
    patho_name = ['S',num2str(Index_num),'_p'];
    thr_con_name = ['S',num2str(Index_num),'_c_thr'];
    thr_patho_name = ['S',num2str(Index_num),'_p_thr'];
    
    %list_name_con={list_name_con, con_name};
    
    connectome_data{i,1} = con_name;
    connectome_data{i,2} = csvread(con_file_name);
    connectome_data{i,3} = patho_name;
    connectome_data{i,4} = csvread(patho_file_name);
    connectome_data{i,5} = thr_con_name;
    connectome_data{i,6} = csvread(thr_con_file_name);
    connectome_data{i,7} = thr_patho_name;
    connectome_data{i,8} = csvread(thr_patho_file_name);
    
    %%%save connectome data in the 'connectome_data'
    %degrees
    degrees{i,1} = con_name;
    degrees{i,2} = degrees_und(csvread(con_file_name));
    degrees{i,3} = patho_name;
    degrees{i,4} = degrees_und(csvread(patho_file_name));
    degrees{i,5} = thr_con_name;
    degrees{i,6} = degrees_und(csvread(thr_con_file_name));
    degrees{i,7} = thr_patho_name;
    degrees{i,8} = degrees_und(csvread(thr_patho_file_name));
end
%% node strength
strength = {}
for i = 1:size(connectome_data,1)
    for e = 1:size(connectome_data,2)
        if mod(e,2)==0
            strength{i,e} =strengths_und(connectome_data{i,e})
        else    
            strength{i,e} = connectome_data{i,e}
        
        end 
    end
end
%% CPL
characteristic_path_length = {}
for i = 1:size(connectome_data,1)
    for e = 1:size(connectome_data,2)
        if mod(e,2)==0
            characteristic_path_length{i,e} =charpath(connectome_data{i,e})
        else    
            characteristic_path_length{i,e} = connectome_data{i,e}
        
        end 
    end
end

%% rich club

rich_club = {}
for i = 1:size(connectome_data,1)
    for e = 1:size(connectome_data,2)
        if mod(e,2)==0
            rich_club{i,e} =rich_club_wu(connectome_data{i,e})
        else    
            rich_club{i,e} = connectome_data{i,e}
        
        end 
    end
end

%% binarize the weighted matrices
connectome_binarise = {}
for i = 1:size(connectome_data,1)
    for e = 1:size(connectome_data,2)
        if mod(e,2)==0
            connectome_binarise{i,e} = weight_conversion(connectome_data{i,e},'binarize')
        else    
            connectome_binarise{i,e} = connectome_data{i,e}
        
        end 
    end
end

%% get components based on binary matrices
component = {}
for i = 1:size(connectome_binarise,1)
    for e = 1:size(connectome_binarise,2)
        if mod(e,2)==0
            component{i,e} = get_components(connectome_binarise{i,e})
        else    
            component{i,e} = connectome_binarise{i,e}
        
        end 
    end
end



%% characteristic path length paired ttest
CPL_con = characteristic_path_length(:,2)
CPL_patho = characteristic_path_length(:,4)
CPL_con=cell2mat(CPL_con)
CPL_patho=cell2mat(CPL_patho)
CPL_con_thr = characteristic_path_length(:,6)
CPL_patho_thr = characteristic_path_length(:,8)
CPL_con_thr=cell2mat(CPL_con_thr)
CPL_patho_thr=cell2mat(CPL_patho_thr)
[h,p,ci,stats]= ttest(CPL_con_thr,CPL_patho_thr)
[h,p,ci,stats]= ttest(CPL_con,CPL_patho)
boxplot([CPL_con,CPL_patho])

mean(CPL_con)
std(CPL_con)
mean(CPL_patho)
std(CPL_patho)
mean(CPL_con_thr)
std(CPL_con_thr)
mean(CPL_patho_thr)
std(CPL_patho_thr)


%% paired boxplot


%%% Load some sample data:

subplot(1,2,1)

measures = [CPL_con,CPL_patho];
nCats = 2;
nDatas = 20;
%%% Plot
% figure();
boxplot(measures(1:nDatas, 1:nCats),'Labels',{'Healthy hemisphere','Pathologic hemisphere'}, 'Symbol', 'k.'); hold on;
line(repmat([(1:nCats).';NaN], [nDatas,1]), ...
  reshape(measures(1:nDatas,[1:nCats, 1]).', [], 1), ...
  'Color', 0.7*[1 1 1], 'Marker', '.', 'MarkerSize', 10);


xlabel('Hemispheres')
ylabel('Average shortest path length')
title('A Characteristic path length (no threshold)')

subplot(1,2,2)
measures_thr = [CPL_con_thr,CPL_patho_thr];


%nCats = 2;
%nDatas = 20;
%%% Plot
% figure();
boxplot(measures_thr(1:nDatas, 1:nCats),'Labels',{'Healthy hemisphere','Pathologic hemisphere'}, 'Symbol', 'k.'); hold on;
line(repmat([(1:nCats).';NaN], [nDatas,1]), ...
  reshape(measures_thr(1:nDatas,[1:nCats, 1]).', [], 1), ...
  'Color', 0.7*[1 1 1], 'Marker', '.', 'MarkerSize', 10);

xlabel('Hemispheres')
ylabel('Average shortest path length')
title('B Characteristic path length (minimal connected network threshold)')





%% average strength and degrees & paired ttest
strength_mean = {}
for i = 1:size(strength,1)
    for e = 1:size(strength,2)
        if mod(e,2)==0
            strength_mean{i,e} = mean(strength{i,e})
        else    
            strength_mean{i,e} = strength{i,e}
        
        end 
    end
end

strength_mean_con = strength_mean(:,2);
strength_mean_patho = strength_mean(:,4);
strength_mean_con = cell2mat(strength_mean_con);
strength_mean_patho = cell2mat(strength_mean_patho);
strength_mean_con_thr = strength_mean(:,6);
strength_mean_patho_thr = strength_mean(:,8);
strength_mean_con_thr = cell2mat(strength_mean_con_thr);
strength_mean_patho_thr = cell2mat(strength_mean_patho_thr);
[h,p,ci,stats] = ttest(strength_mean_con_thr,strength_mean_patho_thr)
[h,p,ci,stats] = ttest(strength_mean_con,strength_mean_patho)
boxplot([strength_mean_con,strength_mean_patho])
%%%%%we found the average strength is the chracteristic path length times 37

%% paired ttest over each node of strength
strength_node_ttest = {};
strength_node = {};
for node = 1:38
    for subj = 1:20
        strength_node{node,1}(subj,1) = strength{subj,2}(node);
        strength_node{node,2}(subj,1) = strength{subj,4}(node);
        strength_node{node,3}(subj,1) = strength{subj,6}(node);
        strength_node{node,4}(subj,1) = strength{subj,8}(node);
    end
    [h,p,ci,stats] = ttest(strength_node{node,1},strength_node{node,2});
    strength_node_ttest{1,1}(node,1) = h;
    strength_node_ttest{1,2} = 'hypothesis'
    strength_node_ttest{2,1}(node,1) = p;
    strength_node_ttest{2,2} = 'pvalues';
    strength_node_ttest{3,1}{node,1} = ci;
    strength_node_ttest{3,2} = 'CI';
    strength_node_ttest{4,1}{node,1} = stats;
    strength_node_ttest{4,2} = 'stats';
    
    [h,p,ci,stats] = ttest(strength_node{node,3},strength_node{node,4});
    strength_node_ttest{1,3}(node,1) = h;
    strength_node_ttest{1,4} = 'hypothesis_thr'
    strength_node_ttest{2,3}(node,1) = p;
    strength_node_ttest{2,4} = 'pvalues_thr';
    strength_node_ttest{3,3}{node,1} = ci;
    strength_node_ttest{3,4} = 'CI_thr';
    strength_node_ttest{4,3}{node,1} = stats;
    strength_node_ttest{4,4} = 'stats_thr';
    
    
end
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(strength_node_ttest{2,1},0.05,'pdep','yes')
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(strength_node_ttest{2,3},0.05,'pdep','yes')

[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(strength_node_ttest{2,1},0.05,'pdep','yes')
strength_node_ttest{5,1}{1,1}=h
strength_node_ttest{5,1}{2,1}=crit_p;
strength_node_ttest{5,1}{3,1}=adj_ci_cvrg;
strength_node_ttest{5,1}{3,1}=adj_p;
strength_node_ttest{5,2}='adj_fdr';
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(strength_node_ttest{2,3},0.05,'pdep','yes')
strength_node_ttest{5,3}{1,1}=h
strength_node_ttest{5,3}{2,1}=crit_p;
strength_node_ttest{5,3}{3,1}=adj_ci_cvrg;
strength_node_ttest{5,3}{4,1}=adj_p;
strength_node_ttest{5,4}='daj_fdr_thr';
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(strength_node_ttest{2,1},0.05,'pdep','yes')

%% degree
degree_node_ttest = {};
degree_node = {};

for node = 1:38
    for subj = 1:20
        degree_node{node,1}(subj,1) = degrees{subj,2}(node);
        degree_node{node,2}(subj,1) = degrees{subj,4}(node);
        degree_node{node,3}(subj,1) = degrees{subj,6}(node);
        degree_node{node,4}(subj,1) = degrees{subj,8}(node);
    end
    [h,p,ci,stats] = ttest(degree_node{node,1},degree_node{node,2});
    degree_node_ttest{1,1}(node,1) = h;
    degree_node_ttest{1,2} = 'hypothesis'
    degree_node_ttest{2,1}(node,1) = p;
    degree_node_ttest{2,2} = 'pvalues';
    degree_node_ttest{3,1}{node,1} = ci;
    degree_node_ttest{3,2} = 'CI';
    degree_node_ttest{4,1}{node,1} = stats;
    degree_node_ttest{4,2} = 'stats';
    
    [h,p,ci,stats] = ttest(degree_node{node,3},degree_node{node,4});
    degree_node_ttest{1,3}(node,1) = h;
    degree_node_ttest{1,4} = 'hypothesis_thr'
    degree_node_ttest{2,3}(node,1) = p;
    degree_node_ttest{2,4} = 'pvalues_thr';
    degree_node_ttest{3,3}{node,1} = ci;
    degree_node_ttest{3,4} = 'CI_thr';
    degree_node_ttest{4,3}{node,1} = stats;
    degree_node_ttest{4,4} = 'stats_thr';
    
    
end
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(degree_node_ttest{2,1},0.05,'pdep','yes')
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(degree_node_ttest{2,3},0.05,'pdep','yes')



[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(degree_node_ttest{2,1},0.05,'pdep','yes')
degree_node_ttest{5,1}{1,1}=h
degree_node_ttest{5,1}{2,1}=crit_p;
degree_node_ttest{5,1}{3,1}=adj_ci_cvrg;
degree_node_ttest{5,1}{4,1}=adj_p;
degree_node_ttest{5,2}='adj_fdr';
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(degree_node_ttest{2,3},0.05,'pdep','yes')
degree_node_ttest{5,3}{1,1}=h
degree_node_ttest{5,3}{2,1}=crit_p;
degree_node_ttest{5,3}{3,1}=adj_ci_cvrg;
degree_node_ttest{5,3}{4,1}=adj_p;
degree_node_ttest{5,4}='daj_fdr_thr';
%% average degrees & paired ttest
degree_mean = {}
for i = 1:size(degrees,1)
    for e = 1:size(degrees,2)
        if mod(e,2)==0
            degree_mean{i,e} = mean(degrees{i,e})
        else    
            degree_mean{i,e} = degrees{i,e}
        
        end 
    end
end

degree_mean_con = degree_mean(:,2);
degree_mean_patho = degree_mean(:,4);
degree_mean_con = cell2mat(degree_mean_con);
degree_mean_patho = cell2mat(degree_mean_patho);
degree_mean_con_thr = degree_mean(:,6);
degree_mean_patho_thr = degree_mean(:,8);
degree_mean_con_thr = cell2mat(degree_mean_con_thr);
degree_mean_patho_thr = cell2mat(degree_mean_patho_thr);
[h,p,ci,stats] = ttest(degree_mean_con_thr,degree_mean_patho_thr)
[h,p,ci,stats] = ttest(degree_mean_con,degree_mean_patho)
boxplot([strength_mean_con,strength_mean_patho])



%% boxplot


%% local efficiency
local_efficiency= {};
for i = 1:size(connectome_data,1)
    for e = 1:size(connectome_data,2)
        if mod(e,2)==0
            local_efficiency{i,e} = efficiency_wei(connectome_data{i,e},2)
        else    
            local_efficiency{i,e} = connectome_data{i,e}
        
        end 
    end
end


local_efficiency_node_ttest = {};
local_efficiency_node = {};
for node = 1:38
    for subj = 1:20
        local_efficiency_node{node,1}(subj,1) = local_efficiency{subj,2}(node);
        local_efficiency_node{node,2}(subj,1) = local_efficiency{subj,4}(node);
        local_efficiency_node{node,3}(subj,1) = local_efficiency{subj,6}(node);
        local_efficiency_node{node,4}(subj,1) = local_efficiency{subj,8}(node);
    end
    [h,p,ci,stats] = ttest(local_efficiency_node{node,1},local_efficiency_node{node,2});
    local_efficiency_node_ttest{1,1}(node,1) = h;
    local_efficiency_node_ttest{1,2} = 'hypothesis'
    local_efficiency_node_ttest{2,1}(node,1) = p;
    local_efficiency_node_ttest{2,2} = 'pvalues';
    local_efficiency_node_ttest{3,1}{node,1} = ci;
    local_efficiency_node_ttest{3,2} = 'CI';
    local_efficiency_node_ttest{4,1}{node,1} = stats;
    local_efficiency_node_ttest{4,2} = 'stats';
    
    [h,p,ci,stats] = ttest(local_efficiency_node{node,3},local_efficiency_node{node,4});
    local_efficiency_node_ttest{1,3}(node,1) = h;
    local_efficiency_node_ttest{1,4} = 'hypothesis_thr'
    local_efficiency_node_ttest{2,3}(node,1) = p;
    local_efficiency_node_ttest{2,4} = 'pvalues_thr';
    local_efficiency_node_ttest{3,3}{node,1} = ci;
    local_efficiency_node_ttest{3,4} = 'CI_thr';
    local_efficiency_node_ttest{4,3}{node,1} = stats;
    local_efficiency_node_ttest{4,4} = 'stats_thr';
    
    
end
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(local_efficiency_node_ttest{2,1},0.05,'pdep','yes')
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(local_efficiency_node_ttest{2,3},0.05,'pdep','yes')

[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(local_efficiency_node_ttest{2,1},0.05,'pdep','yes')
local_efficiency_node_ttest{5,1}{1,1}=h
local_efficiency_node_ttest{5,1}{2,1}=crit_p;
local_efficiency_node_ttest{5,1}{3,1}=adj_ci_cvrg;
local_efficiency_node_ttest{5,1}{3,1}=adj_p;
local_efficiency_node_ttest{5,2}='adj_fdr';
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(local_efficiency_node_ttest{2,3},0.05,'pdep','yes')
local_efficiency_node_ttest{5,3}{1,1}=h
local_efficiency_node_ttest{5,3}{2,1}=crit_p;
local_efficiency_node_ttest{5,3}{3,1}=adj_ci_cvrg;
local_efficiency_node_ttest{5,3}{4,1}=adj_p;
local_efficiency_node_ttest{5,4}='daj_fdr_thr';
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(local_efficiency_node_ttest{2,1},0.05,'pdep','yes')


%%
local_efficiency_mean = {}
for i = 1:size(local_efficiency,1)
    for e = 1:size(local_efficiency,2)
        if mod(e,2)==0
            local_efficiency_mean{i,e} = mean(local_efficiency{i,e})
        else    
            local_efficiency_mean{i,e} = local_efficiency{i,e}
        
        end 
    end
end

local_efficiency_mean_con = local_efficiency_mean(:,2);
local_efficiency_mean_patho = local_efficiency_mean(:,4);
local_efficiency_mean_con = cell2mat(local_efficiency_mean_con);
local_efficiency_mean_patho = cell2mat(local_efficiency_mean_patho);
local_efficiency_mean_con_thr = local_efficiency_mean(:,6);
local_efficiency_mean_patho_thr = local_efficiency_mean(:,8);
local_efficiency_mean_con_thr = cell2mat(local_efficiency_mean_con_thr);
local_efficiency_mean_patho_thr = cell2mat(local_efficiency_mean_patho_thr);
[h,p,ci,stats] = ttest(local_efficiency_mean_con_thr,local_efficiency_mean_patho_thr)
[h,p,ci,stats] = ttest(local_efficiency_mean_con,local_efficiency_mean_patho)








%% global efficiency
global_efficiency= {};
for i = 1:size(connectome_data,1)
    for e = 1:size(connectome_data,2)
        if mod(e,2)==0
            global_efficiency{i,e} = efficiency_wei(connectome_data{i,e})
        else    
            global_efficiency{i,e} = connectome_data{i,e}
        
        end 
    end
end




global_efficiency_con = global_efficiency(:,2);
global_efficiency_patho = global_efficiency(:,4);
global_efficiency_con = cell2mat(global_efficiency_con);
global_efficiency_patho = cell2mat(global_efficiency_patho);
global_efficiency_con_thr = global_efficiency(:,6);
global_efficiency_patho_thr = global_efficiency(:,8);
global_efficiency_con_thr = cell2mat(global_efficiency_con_thr);
global_efficiency_patho_thr = cell2mat(global_efficiency_patho_thr);
[h,p,ci,stats] = ttest(global_efficiency_con_thr,global_efficiency_patho_thr)
[h,p,ci,stats] = ttest(global_efficiency_con,global_efficiency_patho)


%% Clustering coefficient

cluster_coefficient= {};
for i = 1:size(connectome_data,1)
    for e = 1:size(connectome_data,2)
        if mod(e,2)==0
            cluster_coefficient{i,e} = clustering_coef_wu(connectome_data{i,e})
        else    
            cluster_coefficient{i,e} = connectome_data{i,e}
        
        end 
    end
end


cluster_coefficient_node_ttest = {};
cluster_coefficient_node = {};
for node = 1:38
    for subj = 1:20
        cluster_coefficient_node{node,1}(subj,1) = cluster_coefficient{subj,2}(node);
        cluster_coefficient_node{node,2}(subj,1) = cluster_coefficient{subj,4}(node);
        cluster_coefficient_node{node,3}(subj,1) = cluster_coefficient{subj,6}(node);
        cluster_coefficient_node{node,4}(subj,1) = cluster_coefficient{subj,8}(node);
    end
    [h,p,ci,stats] = ttest(cluster_coefficient_node{node,1},cluster_coefficient_node{node,2});
    cluster_coefficient_node_ttest{1,1}(node,1) = h;
    cluster_coefficient_node_ttest{1,2} = 'hypothesis'
    cluster_coefficient_node_ttest{2,1}(node,1) = p;
    cluster_coefficient_node_ttest{2,2} = 'pvalues';
    cluster_coefficient_node_ttest{3,1}{node,1} = ci;
    cluster_coefficient_node_ttest{3,2} = 'CI';
    cluster_coefficient_node_ttest{4,1}{node,1} = stats;
    cluster_coefficient_node_ttest{4,2} = 'stats';
    
    [h,p,ci,stats] = ttest(cluster_coefficient_node{node,3},cluster_coefficient_node{node,4});
    cluster_coefficient_node_ttest{1,3}(node,1) = h;
    cluster_coefficient_node_ttest{1,4} = 'hypothesis_thr'
    cluster_coefficient_node_ttest{2,3}(node,1) = p;
    cluster_coefficient_node_ttest{2,4} = 'pvalues_thr';
    cluster_coefficient_node_ttest{3,3}{node,1} = ci;
    cluster_coefficient_node_ttest{3,4} = 'CI_thr';
    cluster_coefficient_node_ttest{4,3}{node,1} = stats;
    cluster_coefficient_node_ttest{4,4} = 'stats_thr';
    
    
end
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(cluster_coefficient_node_ttest{2,1},0.05,'pdep','yes')
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(cluster_coefficient_node_ttest{2,3},0.05,'pdep','yes')

[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(cluster_coefficient_node_ttest{2,1},0.05,'pdep','yes')
cluster_coefficient_node_ttest{5,1}{1,1}=h
cluster_coefficient_node_ttest{5,1}{2,1}=crit_p;
cluster_coefficient_node_ttest{5,1}{3,1}=adj_ci_cvrg;
cluster_coefficient_node_ttest{5,1}{4,1}=adj_p;
cluster_coefficient_node_ttest{5,2}='adj_fdr';
[h, crit_p, adj_ci_cvrg, adj_p] =fdr_bh(cluster_coefficient_node_ttest{2,3},0.05,'pdep','yes')
cluster_coefficient_node_ttest{5,3}{1,1}=h
cluster_coefficient_node_ttest{5,3}{2,1}=crit_p;
cluster_coefficient_node_ttest{5,3}{3,1}=adj_ci_cvrg;
cluster_coefficient_node_ttest{5,3}{4,1}=adj_p;
cluster_coefficient_node_ttest{5,4}='daj_fdr_thr';


%% average cluster_coefficient
cluster_coefficient_mean = {}
for i = 1:size(cluster_coefficient,1)
    for e = 1:size(cluster_coefficient,2)
        if mod(e,2)==0
            cluster_coefficient_mean{i,e} = mean(cluster_coefficient{i,e})
        else    
            cluster_coefficient_mean{i,e} = cluster_coefficient{i,e}
        
        end 
    end
end

cluster_coefficient_mean_con = cluster_coefficient_mean(:,2);
cluster_coefficient_mean_patho = cluster_coefficient_mean(:,4);
cluster_coefficient_mean_con = cell2mat(cluster_coefficient_mean_con);
cluster_coefficient_mean_patho = cell2mat(cluster_coefficient_mean_patho);
cluster_coefficient_mean_con_thr = cluster_coefficient_mean(:,6);
cluster_coefficient_mean_patho_thr = cluster_coefficient_mean(:,8);
cluster_coefficient_mean_con_thr = cell2mat(cluster_coefficient_mean_con_thr);
cluster_coefficient_mean_patho_thr = cell2mat(cluster_coefficient_mean_patho_thr);
[h,p,ci,stats] = ttest(cluster_coefficient_mean_con_thr,cluster_coefficient_mean_patho_thr)
[h,p,ci,stats] = ttest(cluster_coefficient_mean_con,cluster_coefficient_mean_patho)
%% assortativity paired ttest
assortativity = {}
for i = 1:size(connectome_data,1)
    for e = 1:size(connectome_data,2)
        if mod(e,2)==0
            assortativity{i,e} =assortativity_wei(connectome_data{i,e},0)
        else    
            assortativity{i,e} = connectome_data{i,e}
        
        end 
    end
end

assortativity_con = assortativity(:,2)
assortativity_patho = assortativity(:,4)
assortativity_con=cell2mat(assortativity_con)
assortativity_patho=cell2mat(assortativity_patho)
assortativity_con_thr = assortativity(:,6)
assortativity_patho_thr = assortativity(:,8)
assortativity_con_thr=cell2mat(assortativity_con_thr)
assortativity_patho_thr=cell2mat(assortativity_patho_thr)
[h,p,ci,stats]= ttest(assortativity_con_thr,assortativity_patho_thr)
[h,p,ci,stats]= ttest(assortativity_con,assortativity_patho)
boxplot([assortativity_con,assortativity_patho])


mean(CPL_con)
std(CPL_con)
mean(CPL_patho)
std(CPL_patho)
mean(CPL_con_thr)
std(CPL_con_thr)
mean(CPL_patho_thr)
std(CPL_patho_thr)
%% average the matrices into one
connectome_average={};
connectome_average_con=zeros(38,38);
connectome_average_patho=zeros(38,38);
connectome_average_con_thr=zeros(38,38);
connectome_average_patho_thr=zeros(38,38);

for subj = 1:20
    connectome_average_con = connectome_average_con + connectome_data{subj,2};
    connectome_average_patho = connectome_average_patho + connectome_data{subj,4};
    connectome_average_con_thr = connectome_average_con_thr + connectome_data{subj,6};
    connectome_average_patho_thr = connectome_average_patho_thr + connectome_data{subj,8};
end

connectome_average_con = connectome_average_con/20;
connectome_average_patho = connectome_average_patho/20;
connectome_average_con_thr = connectome_average_con_thr/20;
connectome_average_patho_thr = connectome_average_patho_thr/20;

%% calculate the degrees per group and test for differences
degree_average_con = degrees_und(connectome_average_con)
degree_average_patho =degrees_und(connectome_average_patho)
degree_average_con_thr = degrees_und(connectome_average_con_thr)
degree_average_patho_thr = degrees_und(connectome_average_patho_thr)
[h,p,ci,stats]=ttest(degree_average_con_thr',degree_average_patho_thr')

%% calculate the strength per group and test for differences
strength_average_con =strengths_und(connectome_average_con)
strength_average_patho =strengths_und(connectome_average_patho)
strength_average_con_thr =strengths_und(connectome_average_con_thr)
strength_average_patho_thr =strengths_und(connectome_average_patho_thr)
[h,p,ci,stats]=ttest(strength_average_con_thr',strength_average_patho_thr')
[h,p,ci,stats]=ttest(strength_average_con',strength_average_patho')

%% calculate the assortativity per group and test for differences
assortativity_average_con =assortativity_wei(connectome_average_con,0)
assortativity_average_patho =assortativity_wei(connectome_average_patho,0)
assortativity_average_con_thr =assortativity_wei(connectome_average_con_thr,0)
assortativity_average_patho_thr =assortativity_wei(connectome_average_patho_thr,0)



%% calculate the global efficiency
efficiency_average_con =efficiency_wei(connectome_average_con)
efficiency_average_patho =efficiency_wei(connectome_average_patho)
efficiency_average_con_thr =efficiency_wei(connectome_average_con_thr)
efficiency_average_patho_thr =efficiency_wei(connectome_average_patho_thr)

%% calculate the local efficiency
local_efficiency_average_con =efficiency_wei(connectome_average_con,2)
local_efficiency_average_patho =efficiency_wei(connectome_average_patho,2)
local_efficiency_average_con_thr =efficiency_wei(connectome_average_con_thr,2)
local_efficiency_average_patho_thr =efficiency_wei(connectome_average_patho_thr,2)
[h,p,ci,stats]=ttest(local_efficiency_average_con_thr',local_efficiency_average_patho_thr')
[h,p,ci,stats]=ttest(local_efficiency_average_con',local_efficiency_average_patho')


%%rich_club
rich_club_average_con =rich_club_wu(connectome_average_con)
rich_club_average_patho =rich_club_wu(connectome_average_patho)
rich_club_average_con_thr =rich_club_wu(connectome_average_con_thr)
rich_club_average_patho_thr =rich_club_wu(connectome_average_patho_thr)

%% clustering coefficient
clustering_coef_average_con =clustering_coef_wu(connectome_average_con)
clustering_coef_average_patho =clustering_coef_wu(connectome_average_patho)
clustering_coef_average_con_thr =clustering_coef_wu(connectome_average_con_thr)
clustering_coef_average_patho_thr =clustering_coef_wu(connectome_average_patho_thr)
[h,p,ci,stats]=ttest(clustering_coef_average_con_thr',clustering_coef_average_patho_thr')
[h,p,ci,stats]=ttest(clustering_coef_average_con',clustering_coef_average_patho')




high_grade = [28,81,82,110,130,149,155,164,293,374,382,400];
low_grad = [208,284,285,286,287,291,390,391];

high_grade_index = [1,2,3,4,5,6,7,8,15,16,17,20];
low_grade_index = [9,10,11,12,13,14,18,19];

high_grade_assortativity_con=[]
low_grade_assortativity_con=[]
for i = 1:20
    if any(high_grade_index==i)
        high_grade_assortativity_con = [high_grade_assortativity_con,assortativity{i,2}];
    else
        low_grade_assortativity_con = [low_grade_assortativity_con,assortativity{i,2}];
  
    end
end
[h,p,ci,stats]=ttest2(high_grade_assortativity_con',low_grade_assortativity_con')
    
high_grade_assortativity_patho=[]
low_grade_assortativity_patho=[]
for i = 1:20
    if any(high_grade_index==i)
        high_grade_assortativity_patho = [high_grade_assortativity_patho,assortativity{i,4}];
    else
        low_grade_assortativity_patho = [low_grade_assortativity_patho,assortativity{i,4}];
  
    end
end
[h,p,ci,stats]=ttest2(high_grade_assortativity_patho',low_grade_assortativity_patho')


high_grade_assortativity_con_thr=[]
low_grade_assortativity_con_thr=[]
for i = 1:20
    if any(high_grade_index==i)
        high_grade_assortativity_con_thr = [high_grade_assortativity_con_thr,assortativity{i,6}];
    else
        low_grade_assortativity_con_thr = [low_grade_assortativity_con_thr,assortativity{i,6}];
  
    end
end
[h,p,ci,stats]=ttest2(high_grade_assortativity_con_thr',low_grade_assortativity_con_thr')
    
high_grade_assortativity_patho_thr=[]
low_grade_assortativity_patho_thr=[]
for i = 1:20
    if any(high_grade_index==i)
        high_grade_assortativity_patho_thr = [high_grade_assortativity_patho_thr,assortativity{i,8}];
    else
        low_grade_assortativity_patho_thr = [low_grade_assortativity_patho_thr,assortativity{i,8}];
  
    end
end
[h,p,ci,stats]=ttest2(high_grade_assortativity_patho_thr',low_grade_assortativity_patho_thr')










high_grade_global_efficiency_con=[]
low_grade_global_efficiency_con=[]
for i = 1:20
    if any(high_grade_index==i)
        high_grade_global_efficiency_con = [high_grade_global_efficiency_con,global_efficiency{i,2}];
    else
        low_grade_global_efficiency_con = [low_grade_global_efficiency_con,global_efficiency{i,2}];
  
    end
end
[h,p,ci,stats]=ttest2(high_grade_global_efficiency_con',low_grade_global_efficiency_con')
    
high_grade_global_efficiency_patho=[]
low_grade_global_efficiency_patho=[]
for i = 1:20
    if any(high_grade_index==i)
        high_grade_global_efficiency_patho = [high_grade_global_efficiency_patho,global_efficiency{i,4}];
    else
        low_grade_global_efficiency_patho = [low_grade_global_efficiency_patho,global_efficiency{i,4}];
  
    end
end
[h,p,ci,stats]=ttest2(high_grade_global_efficiency_patho',low_grade_global_efficiency_patho')


high_grade_global_efficiency_con_thr=[]
low_grade_global_efficiency_con_thr=[]
for i = 1:20
    if any(high_grade_index==i)
        high_grade_global_efficiency_con_thr = [high_grade_global_efficiency_con_thr,global_efficiency{i,6}];
    else
        low_grade_global_efficiency_con_thr = [low_grade_global_efficiency_con_thr,global_efficiency{i,6}];
  
    end
end
[h,p,ci,stats]=ttest2(high_grade_global_efficiency_con_thr',low_grade_global_efficiency_con_thr')
    
high_grade_global_efficiency_patho_thr=[]
low_grade_global_efficiency_patho_thr=[]
for i = 1:20
    if any(high_grade_index==i)
        high_grade_global_efficiency_patho_thr = [high_grade_global_efficiency_patho_thr,global_efficiency{i,8}];
    else
        low_grade_global_efficiency_patho_thr = [low_grade_global_efficiency_patho_thr,global_efficiency{i,8}];
  
    end
end
[h,p,ci,stats]=ttest2(high_grade_global_efficiency_patho_thr',low_grade_global_efficiency_patho_thr')


%% correlation
%%%All: read p-value file
clc;clear;
p_values_all = csvread('tfnbs/tfnbs_fwe_1mpvalue.csv',1,0);
mask_all = p_values_all > 0.95;
ltria_mask = tril(mask_all);
list_num = [15,28,31,32,81,82,110,112,130,139,149,155,164,180,208,212,284,285,286,287,290,291,292,293,299,302,354,356,364,374,381,382,390,391,395,400,401,407];
values_all = {};
for i = 1:length(list_num)
    Index_num = list_num(i);
    con_file_name = ['S',num2str(Index_num),'_c.csv'];
    patho_file_name = ['S',num2str(Index_num),'_p.csv'];
   
    %read csv and call computing function
    connectome_c = csvread(con_file_name);
    connectome_p = csvread(patho_file_name);
    connectome_c = connectome_c .* ltria_mask;
    connectome_p = connectome_p .* ltria_mask;
    connectome_c = connectome_c(:);
    connectome_p = connectome_p(:);
    
    
    
    values_all{i,1} = connectome_c;
    values_all{i,2} = con_file_name;
    values_all{i,3} = connectome_p;
    values_all{i,4} = patho_file_name;   
end

motor_status = [3
5
4
5
4
5
4
5
4
4
2
4
5
5
4
5
5
5
5
4
5
5
5
5
5
5
5
5
5
3
4
5
4
5
4
5
5
4
]

values_all_vectors = {};

for i = 1:1600
    value_all_p = [];
    value_all_c = [];
    value_all_difference = [];
    for e = 1:38
        value_all_p=[value_all_p;values_all{e,3}(i)];
        value_all_c=[value_all_c;values_all{e,1}(i)];
        value_all_difference=[value_all_difference;values_all{e,1}(i)-values_all{e,3}(i)];
    end
    values_all_vectors{i,1}=value_all_p;
    values_all_vectors{i,2}=value_all_c;
    values_all_vectors{i,3}=value_all_difference;
end

cor_coefficient_p = [];
cor_pvalues_p = [];
for i = 1:1600
    edge_vector = values_all_vectors{i,1};
    [RHO,PVAL] = corr(edge_vector,motor_status,'Type','Spearman');
    cor_coefficient_p = [cor_coefficient_p;RHO];
    cor_pvalues_p = [cor_pvalues_p;PVAL];
    
end

cor_coefficient_matrix_p = reshape(cor_coefficient_p,[40 40]);
cor_pvalues_matrix_p = reshape(cor_pvalues_p,[40,40]);




cor_coefficient_c = [];
cor_pvalues_c = [];
for i = 1:1600
    edge_vector = values_all_vectors{i,2};
    [RHO,PVAL] = corr(edge_vector,motor_status,'Type','Spearman');
    cor_coefficient_c = [cor_coefficient_c;RHO];
    cor_pvalues_c = [cor_pvalues_c;PVAL];
    
end

cor_coefficient_matrix_c = reshape(cor_coefficient_c,[40 40]);
cor_pvalues_matrix_c = reshape(cor_pvalues_c,[40,40]);

cor_coefficient_difference = [];
cor_pvalues_difference = [];
for i = 1:1600
    edge_vector = values_all_vectors{i,3};
    [RHO,PVAL] = corr(edge_vector,motor_status,'Type','Spearman');
    cor_coefficient_difference = [cor_coefficient_difference;RHO];
    cor_pvalues_difference = [cor_pvalues_difference;PVAL];
    
end



%%%%difference 40-15, difference is at the third column
h = figure ('Color', [1 1 1]);
s1 = plot(motor_status, values_all_vectors{600,3},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 1);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('paracentral – brainstem', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength difference', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(2.5,400,['rs(38) = -.3341; p = .0404'],'FontSize',10)


%%%%difference 25-2
h = figure ('Color', [1 1 1]);
s1 = plot(motor_status, values_all_vectors{65,3},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 1);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('caudal middle frontal – rostral middle frontal', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength difference', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(2.5,2000,['rs(38) = -3326; p = .0483'],'FontSize',10)


%%%%difference 34-12
h = figure ('Color', [1 1 1]);
s1 = plot(motor_status, values_all_vectors{474,3},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 1);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('medial orbitofrontal – Thalamus proper', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength difference', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(2.5,1200,['rs(38) = .3370; p = .0386'],'FontSize',10)



%%%%pathologic 18-15 pathologic is in the first column
h = figure ('Color', [1 1 1]);
s1 = plot(motor_status, values_all_vectors{578,1},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 1);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('paracentral – pars triangularis', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')

text(2.5,15,['rs(38) = .4467; p = .0049'],'FontSize',10)

%% correlation of precentral patients
%%%All: read p-value file
p_values_pre = csvread('tfnbs/tfnbs_precentral_fwe_1mpvalue.csv',1,0);
mask_pre = p_values_pre > 0.95;
ltria_pre = tril(mask_pre);
list_num_pre = [15,81,82,110,112,139,208,286,292,293,299,354,382,391,395,400,407];
values_pre = {};
for i = 1:length(list_num_pre)
    Index_num = list_num_pre(i);
    con_file_name = ['S',num2str(Index_num),'_c.csv'];
    patho_file_name = ['S',num2str(Index_num),'_p.csv'];
   
    %read csv and call computing function
    connectome_c = csvread(con_file_name);
    connectome_p = csvread(patho_file_name);
    connectome_c = connectome_c .* ltria_mask;
    connectome_p = connectome_p .* ltria_mask;
    connectome_c = connectome_c(:);
    connectome_p = connectome_p(:);
    
    
    
    values_pre{i,1} = connectome_c;
    values_pre{i,2} = con_file_name;
    values_pre{i,3} = connectome_p;
    values_pre{i,4} = patho_file_name;   
end
motor_status_pre=[3;4;5;4;5;4;4;5;5;5;5;5;5;5;4;5;4]

values_pre_vectors = {};

for i = 1:1600
    value_pre_p = [];
    value_pre_c = [];
    value_pre_difference = [];
    for e = 1:17
        value_pre_p=[value_pre_p;values_pre{e,3}(i)];
        value_pre_c=[value_pre_c;values_pre{e,1}(i)];
        value_pre_difference=[value_pre_difference;values_pre{e,1}(i)-values_pre{e,3}(i)];
    end
    values_pre_vectors{i,1}=value_pre_p;
    values_pre_vectors{i,2}=value_pre_c;
    values_pre_vectors{i,3}=value_pre_difference;
end

cor_pre_coefficient_p = [];
cor_pre_pvalues_p = [];
for i = 1:1600
    edge_vector = values_pre_vectors{i,1};
    [RHO,PVAL] = corr(edge_vector,motor_status_pre,'Type','Spearman');
    cor_pre_coefficient_p = [cor_pre_coefficient_p;RHO];
    cor_pre_pvalues_p = [cor_pre_pvalues_p;PVAL];
    
end

cor_pre_coefficient_matrix_p = reshape(cor_pre_coefficient_p,[40 40]);
cor_pre_pvalues_matrix_p = reshape(cor_pre_pvalues_p,[40,40]);




cor_pre_coefficient_c = [];
cor_pre_pvalues_c = [];
for i = 1:1600
    [RHO,PVAL] = corr(values_pre_vectors{i,2},motor_status_pre,'Type','Spearman');
    cor_pre_coefficient_c = [cor_pre_coefficient_c;RHO];
    cor_pre_pvalues_c = [cor_pre_pvalues_c;PVAL];
    
end

cor_pre_coefficient_matrix_c = reshape(cor_pre_coefficient_c,[40 40]);
cor_pre_pvalues_matrix_c = reshape(cor_pre_pvalues_c,[40,40]);

cor_pre_coefficient_difference = [];
cor_pre_pvalues_difference = [];
for i = 1:1600
    edge_vector = values_pre_vectors{i,3};
    [RHO,PVAL] = corr(edge_vector,motor_status_pre,'Type','Spearman');
    cor_pre_coefficient_difference = [cor_pre_coefficient_difference;RHO];
    cor_pre_pvalues_difference = [cor_pre_pvalues_difference;PVAL];
    
end

cor_pre_coefficient_matrix_difference = reshape(cor_pre_coefficient_difference,[40 40]);
cor_pre_pvalues_matrix_difference = reshape(cor_pre_pvalues_difference,[40,40]);

%%%correlation plots
%%%%Healty 29-2
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{69,2},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 1);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('supramarginal – caudal middle frontal', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,500,['rs(17) = .6432; p = .0053'],'FontSize',10)

%%%%Healty 26-6
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{226,2},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 1);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('superior frontal – inferior parietal', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,550,['rs(17) = .5980; p = .0112'],'FontSize',10)

%%%%Healty 26-13
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{506,2},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 1);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('superior frontal – middle temporal', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,300,['rs(17) = .5095; p = .0367'],'FontSize',10)




%%%%Healty 39-15
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{599,2},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('superior frontal – middle temporal', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,90,['rs(17) = -.5719; p = .0164'],'FontSize',10)

%%%%Healty 29-18
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{709,2},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('supramarginal – pars triangularis', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,80,['rs(17) = .5078; p = .0375'],'FontSize',10)



%%%%Healty 34-20
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{794,2},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('putamen – postcentral', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,2500,['rs(17) = .5529; p = .0213'],'FontSize',10)

%%%%Healty 35-22
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{875,2},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('pallidum – precentral', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,3200,['rs(17) = .5529; p = .0213'],'FontSize',10)

%%%%%%%%%pathologic 26-6 
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{226,1},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('superior frontal – inferior parietal', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,200,['rs(17) = .5081; p = .0373'],'FontSize',10)

%%%%%%%%%pathologic 29-18
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{709,1},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('supramarginal – pars triangularis', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,100,['rs(17) = .4869; p = .0475'],'FontSize',10)

%%%%difference 25-2
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_pre, values_pre_vectors{65,3},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('caudal middle frontal – rostral middle frontal', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength difference', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,5500,['rs(17) = .5007; p = .0406'],'FontSize',10)


%% correlation of insular patients
%%%All: read p-value file
p_values_insu = csvread('tfnbs/tfnbs_insular_fwe_1mpvalue.csv',1,0);
mask_insu = p_values_insu > 0.95;
ltria_insu = tril(mask_insu);
list_num_insu = [28,31,32,130,180,290,302,381];
values_insu = {};
for i = 1:length(list_num_insu)
    Index_num = list_num_insu(i);
    con_file_name = ['S',num2str(Index_num),'_c.csv'];
    patho_file_name = ['S',num2str(Index_num),'_p.csv'];
   
    %read csv and call computing function
    connectome_c = csvread(con_file_name);
    connectome_p = csvread(patho_file_name);
    connectome_c = connectome_c .* ltria_mask;
    connectome_p = connectome_p .* ltria_mask;
    connectome_c = connectome_c(:);
    connectome_p = connectome_p(:);
    
    
    
    values_insu{i,1} = connectome_c;
    values_insu{i,2} = con_file_name;
    values_insu{i,3} = connectome_p;
    values_insu{i,4} = patho_file_name;   
end
motor_status_insu=[5;4;5;4;5;5;5;4]

values_insu_vectors = {};

for i = 1:1600
    value_insu_p = [];
    value_insu_c = [];
    value_insu_difference = [];
    for e = 1:8
        value_insu_p=[value_insu_p;values_insu{e,3}(i)];
        value_insu_c=[value_insu_c;values_insu{e,1}(i)];
        value_insu_difference=[value_insu_difference;values_insu{e,1}(i)-values_insu{e,3}(i)];
    end
    values_insu_vectors{i,1}=value_insu_p;
    values_insu_vectors{i,2}=value_insu_c;
    values_insu_vectors{i,3}=value_insu_difference;
end

cor_insu_coefficient_p = [];
cor_insu_pvalues_p = [];
for i = 1:1600
    edge_vector = values_insu_vectors{i,1};
    [RHO,PVAL] = corr(edge_vector,motor_status_insu,'Type','Spearman');
    cor_insu_coefficient_p = [cor_insu_coefficient_p;RHO];
    cor_insu_pvalues_p = [cor_insu_pvalues_p;PVAL];
    
end

cor_insu_coefficient_matrix_p = reshape(cor_insu_coefficient_p,[40 40]);
cor_insu_pvalues_matrix_p = reshape(cor_insu_pvalues_p,[40,40]);


%%%difference
cor_insu_coefficient_c = [];
cor_insu_pvalues_c = [];
for i = 1:1600
    edge_vector = values_insu_vectors{i,2};
    [RHO,PVAL] = corr(edge_vector,motor_status_insu,'Type','Spearman');
    cor_insu_coefficient_c = [cor_insu_coefficient_c;RHO];
    cor_insu_pvalues_c = [cor_insu_pvalues_c;PVAL];
    
end

cor_insu_coefficient_matrix_c = reshape(cor_insu_coefficient_c,[40 40]);
cor_insu_pvalues_matrix_c = reshape(cor_insu_pvalues_c,[40,40]);

cor_insu_coefficient_difference = [];
cor_insu_pvalues_difference = [];
for i = 1:1600
    edge_vector = values_insu_vectors{i,3};
    [RHO,PVAL] = corr(edge_vector,motor_status_insu,'Type','Spearman');
    cor_insu_coefficient_difference = [cor_insu_coefficient_difference;RHO];
    cor_insu_pvalues_difference = [cor_insu_pvalues_difference;PVAL];
    
end

cor_insu_coefficient_matrix_difference = reshape(cor_insu_coefficient_difference,[40 40]);
cor_insu_pvalues_matrix_difference = reshape(cor_insu_pvalues_difference,[40,40]);


%%%%pathologic 26-16
h = figure ('Color', [1 1 1],'Position',[1000,890,700,600]);
s1 = plot(motor_status_insu, values_insu_vectors{626,1},'bo');
set(s1, 'MarkerSize', 4, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
%set the parameters of least square line 
l.Color = 'k';
l.LineWidth = 2
%%% axis display 
title('caudal middle frontal – rostral middle frontal', 'FontSize', 10)
xlabel('Motor status', 'FontSize', 10)
ylabel('Edge strength difference', 'FontSize', 10)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')
text(3.5,100,['rs(17) =-.8452; p = .0357'],'FontSize',10)


%% assortitivity
cd '/media/bwgig/5b18abd4-0119-4492-b321-05825dc35687/along_tract_study/connectome_local_analysis/connectome_bs/matrices_excl'
clear
list_num = [15,28,31,32,81,82,110,112,130,139,149,155,164,180,208,212,284,285,286,287,290,291,292,293,299,302,354,356,364,374,381,382,390,391,395,400,401,407];
Index_all=[];
assortativity_c =[];
assortativity_p =[];
for i = 1:length(list_num)
    %set input and output names of csv files 
    Index_num = list_num(i);
    con_file_name = ['S',num2str(Index_num),'_c.csv'];
    patho_file_name = ['S',num2str(Index_num),'_p.csv'];
    Index_all = [Index_all,Index_num];
    %save_name_con = ['S',num2str(Index_num),'_c_thr.csv'];
    %save_name_patho = ['S',num2str(Index_num),'_p_thr.csv'];
   
    %read csv and call computing function
    connectome_c = csvread(con_file_name);
    connectome_p = csvread(patho_file_name);
    %connecomte_c_thr = csvread(connectome_c);
    %connecomte_p_thr = csvread(connectome_p);
   
    %calculate the measure of assotativity coefficinet
    
    measure_c = assortativity_wei(connectome_c,0);
    measure_p = assortativity_wei(connectome_p,0);
    
    assortativity_c = [assortativity_c;measure_c];
    assortativity_p = [assortativity_p;measure_p];
    
    
end



%% calculate the FA csv, seperate the helthy and pathological hemispheres; interpolate data via mean value
clear;clc;
%read connectome matrices and seperate them into healthy and pathologic
%hemipsheres' connectomes
List = dir('S*.csv')
left_patho = [81 82 139 164 284 293 299 302 354 364 374 381 390 391 395 407]
for i = 1:length(List)
    filename = List(i).name;
    Index = sscanf(filename, 'S%d');
    control_name = ['S',num2str(Index),'_c.csv'];
    pathologic_name = ['S',num2str(Index),'_p.csv'];
    connectome = csvread(filename);
    
    connectome_left = connectome([1:39,79],[1:39,79]);
    connectome_right = connectome(40:79,40:79);
    
    [row_left,col_left] = find(isnan(connectome_left));
    [row_right,col_right] = find(isnan(connectome_right));
    
    if ~isempty(row_left)
        connectome_left_mean = nanmean(connectome_left,1);
        for nan_index = 1:length(row_left);
            connectome_left(row_left(nan_index),col_left(nan_index)) = connectome_left_mean(col_left(nan_index));
        end
    end
    
    if ~isempty(row_right)
        connectome_right_mean = nanmean(connectome_right,1);
        for nan_index = 1:length(row_right);
            connectome_right(row_right(nan_index),col_right(nan_index)) = connectome_right_mean(col_right(nan_index));
        end
    end  
    
    %%% exclude upper triangular
    connectome_left = tril(connectome_left);
    connectome_right = tril(connectome_right);
    %%% copy from lower triangular to upper
    connectome_left=triu(connectome_left.',1) + tril(connectome_left);
    connectome_right=triu(connectome_right.',1) + tril(connectome_right);
    
    
    if any(left_patho==Index)
        csvwrite(pathologic_name,connectome_left);
        csvwrite(control_name,connectome_right);
    else
        csvwrite(pathologic_name,connectome_right);
        csvwrite(control_name,connectome_left);
    end
end