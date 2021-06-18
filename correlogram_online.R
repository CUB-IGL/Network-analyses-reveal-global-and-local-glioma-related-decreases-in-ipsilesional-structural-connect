#############1.
library(readxl)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(labelled)
library(rlang)
library(GGally)
library(xlsx)

##ifod data
setwd("working_directory")

###################sd data
matrices_path = 'folder saving connectome matrices/folder saving connectome matrices/'

###################
matrices_prefix = 'S'
matrices_suffix = '.csv'
alpha = 0.05

patient_info = read_excel('connectome_table_online.xlsx')
name_nodes = read_excel('name_nodes.xlsx')


#set the group
patient_groups = 'all'

#load clinical measures
motor_status = patient_info$`Motor status`
RMT_ratio = patient_info$`RMT ratio`
data_matrix=cbind(motor_status,RMT_ratio)
NIHSS = patient_info$`NIHSS(0-42)`
data_matrix = cbind(data_matrix,NIHSS)
WHO = patient_info$Histology
data_matrix = cbind(data_matrix,WHO)
tumor_volume = patient_info$tumor_volume
data_matrix = cbind(data_matrix,tumor_volume)
data_matrix_org=data_matrix

###########################################----------------------------------- load the matrices and calculate the difference-----------------------------------------------------------#######################
###3
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
#####################################-------------------------------------------------------calculation for difference--------------------------------------------------------------------#######################
data_loop<- patient_info %>% filter(all_patients==1)
tfnbs_p =read.csv(paste0(matrices_path,'tfnbs_fwe_1mpvalue.csv'),header = F,skip=1)
tfnbs_p_mask=tfnbs_p>(1-alpha)
tfnbs_p_mask[upper.tri(tfnbs_p_mask)]<-FALSE
n_correlation = sum(tfnbs_p_mask)

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
pdf(file = "Correlogram_sd_diff_a4_new.pdf",   # The directory you want to save the file in
    width = 8.27, # The width of the plot in inches
    height = 11.69) 
ggcorr(as.data.frame(data_matrix),
       method=c("pairwise","spearman"),
       label = TRUE,
       label_size=1.7,
       size = 2.5,
       hjust = 1,
       layout.exp=5)
dev.off()
############################------------------------------------------------calculation for pathological side---------------------------------------------------------------------#######################################
data_matrix=data_matrix_org
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
pdf(file = "Correlogram_sd_patho_a4_new.pdf",   # The directory you want to save the file in
    width = 8.27, # The width of the plot in inches
    height = 11.69) 
ggcorr(as.data.frame(data_matrix),
       method=c("pairwise","spearman"),
       label = TRUE,
       label_size=1.7,
       size = 2.5,
       hjust = 1,
       layout.exp=5)
dev.off()
# ############---------------------------------------------------------calculation for healthy side -------------------------------------------------------------------------------------------#########################
data_matrix=data_matrix_org
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
pdf(file = "Correlogram_sd_health_a4_new.pdf",   # The directory you want to save the file in
    width = 8.27, # The width of the plot in inches
    height = 11.69) 
ggcorr(as.data.frame(data_matrix),
       method=c("pairwise","spearman"),
       label = TRUE,
       label_size=1.7,
       size = 2.5,
       hjust = 1,
       layout.exp=5)
dev.off()


#########################Compare the healthy and pathological hemispheres based on network metrics##################
####paired t test
t.test(sd_global_correlation$Assortativity_C,sd_global_correlation$Assortativity_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$Hierarchy_C,sd_global_correlation$Hierarchy_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$Eglob_C,sd_global_correlation$Eglob_P,paired = TRUE, alternative = "two.sided")
##significant
t.test(sd_global_correlation$Eloc_C,sd_global_correlation$Eloc_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$SmWo_Sigma_C,sd_global_correlation$SmWo_Sigma_P,paired = TRUE, alternative = "two.sided")
t.test(sd_global_correlation$Synchron_C,sd_global_correlation$Synchron_P,paired = TRUE, alternative = "two.sided")
mean(sd_global_correlation$Eglob_C)
sd(sd_global_correlation$Eglob_C)
summary(sd_global_correlation$Eglob_P)
sd(sd_global_correlation$Eglob_P)
mean(sd_global_correlation$Eloc_C)
sd(sd_global_correlation$Eloc_C)
mean(sd_global_correlation$Eloc_P)
sd(sd_global_correlation$Eloc_P)

t.test(ifod_global_correlation$Assortativity_C,ifod_global_correlation$Assortativity_P,paired = TRUE, alternative = "two.sided")
t.test(ifod_global_correlation$Hierarchy_C,ifod_global_correlation$Hierarchy_P,paired = TRUE, alternative = "two.sided")
t.test(ifod_global_correlation$Eglob_C,ifod_global_correlation$Eglob_P,paired = TRUE, alternative = "two.sided")
t.test(ifod_global_correlation$Eloc_C,ifod_global_correlation$Eloc_P,paired = TRUE, alternative = "two.sided")
mean(ifod_global_correlation$Eglob_C)
sd(ifod_global_correlation$Eglob_C)
mean(ifod_global_correlation$Eglob_P)
sd(ifod_global_correlation$Eglob_P)
mean(ifod_global_correlation$Eloc_C)
sd(ifod_global_correlation$Eloc_C)
mean(ifod_global_correlation$Eloc_P)
sd(ifod_global_correlation$Eloc_P)
p <- ggplot(ifod_global_correlation, x = "Eglob_C")
####Healthyhemisphere has higher global and local coefficiency
t.test(ifod_global_correlation$SmWo_Sigma_C,ifod_global_correlation$SmWo_Sigma_P,paired = TRUE, alternative = "two.sided")
t.test(ifod_global_correlation$Synchron_C,ifod_global_correlation$Synchron_P,paired = TRUE, alternative = "two.sided")

mean(ifod_global_correlation$Synchron_C)
sd(ifod_global_correlation$Synchron_C)
mean(ifod_global_correlation$Synchron_P)
sd(ifod_global_correlation$Synchron_P)
## assign a subject number to each subject
subject = c(1:37)
ifod_global_correlation$subject = subject
sd_global_correlation$subject = subject
####convert wide data to long data 
library(tidyr)
ifod_global_correlation_long <- gather(ifod_global_correlation,metrics,values,Assortativity_C:Synchron_P,factor_key=TRUE)
sd_global_correlation_long <- gather(sd_global_correlation,metrics,values,Assortativity_C:Synchron_P,factor_key=TRUE)
###Eglob
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


###Eloc
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

###Synchronization
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





