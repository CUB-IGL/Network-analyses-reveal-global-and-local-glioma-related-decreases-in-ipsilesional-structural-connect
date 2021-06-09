#############1.
library(readxl)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(labelled)
library(rlang)
library(GGally)
library(xlsx)


#set working directory
setwd("path of the folder to save your figures")
#define path fo connectome matrices (in this case sd_stream results)
matrices_path = 'path to sd_stream connectoms matrices'
#read the nodes names and functional data
patient_info = read_excel('excel file contains the information of patints regarding clinical measures')
name_nodes = read_excel('excel file contains the name of nodes')
##Define th prefix if you have prefix before the number of patients. If not, please define matrices_prefix = ''.
matrices_prefix = 'S'

#load the subgroup
patient_groups= c('all','precentral','postcentral','insular','frontal')

#set the alpha level for tests
alpha = 0.05
matrices_suffix = '.csv'

#patient_groups = c('all','precentral','postcentral','insular','frontal')
#we do it separately
patient_groups = 'all'

#load clinical measures
motor_status = patient_info$`Motor status`
rmt_ratio = patient_info$`RMT ratio`
data_matrix=cbind(motor_status,rmt_ratio)
nihss = patient_info$`NIHSS(0-42)`
data_matrix = cbind(data_matrix,nihss)
histology = patient_info$Histology
data_matrix = cbind(data_matrix,histology)


  ###########################################----------------------------------- load the matrices and calculate the difference-----------------------------------------------------------#######################
for (subj in patient_info$Subj) {
  #use paste0 to paste strings together
  healthy_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_c',matrices_suffix)
  healthy_matrices =read.csv(healthy_matrices_path,header = F)
  pathological_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_p',matrices_suffix)
  pathological_matrices =read.csv(pathological_matrices_path,header = F)
  diff_name = paste0('S',subj,'_diff')
  #assign variable name for calculated difference matrices
  assign(diff_name,healthy_matrices-pathological_matrices)
}
#####################################-------------------------------------------------------calculation of difference--------------------------------------------------------------------#######################
###################################################################################################################################################################################################################
for (subgroup in patient_groups){
  ######################all
  if (subgroup =='all'){
    data_loop<- patient_info %>% filter(all_patients==1)
    tfnbs_p =read.csv(paste0(matrices_path,'tfnbs/tfnbs_fwe_1mpvalue.csv'),header = F,skip=1)
  }
  ######################precentral
  else if (subgroup =='precentral'){
    data_loop<- patient_info %>% filter(precentral==1)
    tfnbs_p =read.csv(paste0(matrices_path,'tfnbs/tfnbs_precentral_fwe_1mpvalue.csv'),header = F,skip=1)
  }
  ######################postcentral
  else if (subgroup =='postcentral'){
    data_loop<- patient_info %>% filter(postcentral==1)
    tfnbs_p =read.csv(paste0(matrices_path,'tfnbs/tfnbs_postcentral_fwe_1mpvalue.csv'),header = F,skip=1)
  }
  ######################insular
  else if (subgroup =='insular'){
    data_loop<- patient_info %>% filter(insular==1)
    tfnbs_p =read.csv(paste0(matrices_path,'tfnbs/tfnbs_insular_fwe_1mpvalue.csv'),header = F,skip=1)
  }
  ######################frontal
  else if (subgroup =='frontal'){
    data_loop<- patient_info %>% filter(frontal==1)
    tfnbs_p =read.csv(paste0(matrices_path,'tfnbs/tfnbs_frontal_fwe_1mpvalue.csv'),header = F,skip=1)
  }
  else {print('There is no this subgroup exist.');break}
  #################################################################
  tfnbs_p_mask=tfnbs_p>(1-alpha)
  tfnbs_p_mask[upper.tri(tfnbs_p_mask)]<-FALSE
  n_correlation = sum(tfnbs_p_mask)
  #################################################################
  diff_edge=c()
  diff_motor_edge_p = c()
  diff_who_edge_p = c()
  diff_nihss_edge_p =c()
  diff_rmt_edge_p = c()
  for (i in 1:40){
    for (e in 1:40){
      if (tfnbs_p_mask[i,e]){
        for (subj in data_loop $Subj) {
          diff_name = paste0('S',subj,'_diff')
          #assign variable name for calculated difference matrices
          diff_edge=c(diff_edge,get(diff_name)[i,e])
        }
        
        assign(paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e]),diff_edge)
        data_matrix=cbind(data_matrix,get(paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e])))
        colnames(data_matrix)[ncol(data_matrix)] <- paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e])
        
        
        diff_edge=c()
      }
    }
  }
  ggpairs(as.data.frame(data_matrix),
          title="correlogram edge difference",
          lower = list(continuous = wrap("points",
                                         #alpha = 0.3,
                                         size= 0.3)),
          upper= list(continuous=wrap("cor",size = 2)))+
    theme_minimal(base_size = 6) 
  ggsave("correlogram_diff.pdf")
  # ###########################------------------------------------------------calculation of pathological side---------------------------------------------------------------------#######################################
  patho_edge=c()
  for (i in 1:40){
    for (e in 1:40){
      if (tfnbs_p_mask[i,e]){
        for (subj in data_loop$Subj) {
          pathological_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_p',matrices_suffix)
          pathological_matrices =read.csv(pathological_matrices_path,header = F)
          patho_edge=c(patho_edge,pathological_matrices[i,e])
        }
        
        assign(paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e]),patho_edge)
        data_matrix=cbind(data_matrix,get(paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e])))
        colnames(data_matrix)[ncol(data_matrix)] <- paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e])
        
        
        patho_edge=c()
      }
    }
  }
  ggpairs(as.data.frame(data_matrix),
          title="correlogram edge pathological",
          lower = list(continuous = wrap("points",
                                         #alpha = 0.3,
                                         size= 0.3)),
          upper= list(continuous=wrap("cor",size = 2)))+
    theme_minimal(base_size = 6) 
  ggsave("correlogram_patho.pdf")
  # ############---------------------------------------------------------calculation of healthy side -------------------------------------------------------------------------------------------#########################
  healthy_edge=c()
  for (i in 1:40){
    for (e in 1:40){
      if (tfnbs_p_mask[i,e]){
        for (subj in data_loop $Subj) {
          healthy_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_c',matrices_suffix)
          healthy_matrices =read.csv(healthy_matrices_path,header = F)
          healthy_edge=c(healthy_edge,healthy_matrices[i,e])
        }
        
        assign(paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e]),healthy_edge)
        data_matrix=cbind(data_matrix,get(paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e])))
        colnames(data_matrix)[ncol(data_matrix)] <- paste0(name_nodes$abbreviation[i],"-",name_nodes$abbreviation[e])
        healthy_edge=c()
      }
    }
  }
  ggpairs(as.data.frame(data_matrix),
          title="correlogram edge healthy",
          lower = list(continuous = wrap("points",
                                         #alpha = 0.3,
                                         size= 0.3)),
          upper= list(continuous=wrap("cor",size = 2)))+
    theme_minimal(base_size = 6) 
  ggsave("correlogram_diff.pdf")
}
###################################correlation between clinical measures and network metrics########################
#####we have non-significant
# setwd("path to your excel file which contains the netowrk metrics")
# sd_global_correlation = read.xlsx2('sd_stream excel file',sheetIndex = 1)
# sd_global_correlation <- mutate_all(sd_global_correlation, function(x) as.numeric(as.character(x)))
# ggpairs(sd_global_correlation,cardinality_threshold=NULL)+
#         #title="correlogram global coefficients")+
#   theme_minimal(base_size = 6) 
# 
# 
# ifod_global_correlation = read.xlsx2('ifod excel file',sheetIndex = 1)
# ifod_global_correlation <- mutate_all(ifod_global_correlation, function(x) as.numeric(as.character(x)))
# ggpairs(ifod_global_correlation,cardinality_threshold=NULL,upper = list(corMethod="spearman"))+
#   theme_minimal(base_size = 6) 
# 
# 
#########################Compare the healthy and pathological hemispheres based on network metrics##################
####load file and convert from string to number
sd_global_correlation = read.xlsx2('sd_stream excel file',sheetIndex = 1)
sd_global_correlation <- mutate_all(sd_global_correlation, function(x) as.numeric(as.character(x)))
ifod_global_correlation = read.xlsx2('ifod excel file',sheetIndex = 1)
ifod_global_correlation <- mutate_all(ifod_global_correlation, function(x) as.numeric(as.character(x)))
####paired t test
t.test(sd_global_correlation$Assortativity_C,sd_global_correlation$Assortativity_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$Hierarchy_C,sd_global_correlation$Hierarchy_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$Eglob_C,sd_global_correlation$Eglob_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$Eloc_C,sd_global_correlation$Eloc_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$SmWo_Sigma_C,sd_global_correlation$SmWo_Sigma_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$Synchron_C,sd_global_correlation$Synchron_P,paired = TRUE, alternative = "two.sided")

t.test(ifod_global_correlation$Assortativity_C,ifod_global_correlation$Assortativity_P,paired = TRUE, alternative = "two.sided")
t.test(ifod_global_correlation$Hierarchy_C,ifod_global_correlation$Hierarchy_P,paired = TRUE, alternative = "two.sided")
t.test(ifod_global_correlation$Eglob_C,ifod_global_correlation$Eglob_P,paired = TRUE, alternative = "two.sided")
t.test(ifod_global_correlation$Eloc_C,ifod_global_correlation$Eloc_P,paired = TRUE, alternative = "two.sided")
p <- ggplot(ifod_global_correlation, x = "Eglob_C")
####Control hemisphere has higher global and local coefficiency
t.test(ifod_global_correlation$SmWo_Sigma_C,ifod_global_correlation$SmWo_Sigma_P,paired = TRUE, alternative = "two.sided")
t.test(ifod_global_correlation$Synchron_C,ifod_global_correlation$Synchron_P,paired = TRUE, alternative = "two.sided")
## add a number for each subject
subject = c(1:'number of sample size')
ifod_global_correlation$subject = subject
sd_global_correlation$subject = subject
####convert wide data to long data 
library(tidyr)
ifod_global_correlation_long <- gather(ifod_global_correlation,metrics,values,Assortativity_C:Synchron_P,factor_key=TRUE)
sd_global_correlation_long <- gather(sd_global_correlation,metrics,values,Assortativity_C:Synchron_P,factor_key=TRUE)
###Eglob box plots
sd_global_correlation_long_Eglob <- sd_global_correlation_long %>% filter(metrics %in% c("Eglob_P","Eglob_C"))
ggpaired(sd_global_correlation_long_Eglob, x = "metrics", y = "values",
         color = "metrics", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(method = "t.test",paired = TRUE,label.x=1.3)+xlab("Groups")

ifod_global_correlation_long_Eglob <- ifod_global_correlation_long %>% filter(metrics %in% c("Eglob_P","Eglob_C"))
ggpaired(ifod_global_correlation_long_Eglob, x = "metrics", y = "values",
         color = "metrics", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(method = "t.test",paired = TRUE,label.x=1.3)+xlab("Groups")


###Eloc box plots
sd_global_correlation_long_Eloc <- sd_global_correlation_long %>% filter(metrics %in% c("Eloc_P","Eloc_C"))
ggpaired(sd_global_correlation_long_Eloc, x = "metrics", y = "values",
         color = "metrics", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(method = "t.test",paired = TRUE,label.x=1.3)+xlab("Groups")


ifod_global_correlation_long_Eloc <- ifod_global_correlation_long %>% filter(metrics %in% c("Eloc_P","Eloc_C"))
ggpaired(ifod_global_correlation_long_Eloc, x = "metrics", y = "values",
         color = "metrics", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(method = "t.test",paired = TRUE,label.x=1.3)+xlab("Groups")

###Synchronization box plots
sd_global_correlation_long_Synchron <- sd_global_correlation_long %>% filter(metrics %in% c("Synchron_P","Synchron_C"))
ggpaired(sd_global_correlation_long_Synchron, x = "metrics", y = "values",
         color = "metrics", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(method = "t.test",paired = TRUE,label.x=1.3)+xlab("Groups")


ifod_global_correlation_long_Synchron <- ifod_global_correlation_long %>% filter(metrics %in% c("Synchron_P","Synchron_C"))
ggpaired(ifod_global_correlation_long_Synchron, x = "metrics", y = "values",
         color = "metrics", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(method = "t.test",paired = TRUE,label.x=1.3)+xlab("Groups")





