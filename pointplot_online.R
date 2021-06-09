#load libraries
library(readxl)
library(ggpubr)
library(ggplot2)
library(dplyr)
#set working directory
setwd("path of the folder to save your figures")
#define path fo connectome matrices (in this case sd_stream results)
matrices_path = 'path to sd_stream connectoms matrices'
##Define th prefix if you have prefix before the number of patients. If not, please define matrices_prefix = ''.
matrices_prefix = 'S'

#read the nodes names and functional data
patient_info = read_excel('excel file contains the information of patints regarding clinical measures')
name_nodes = read_excel('excel file contains the name of nodes')

#load the subgroup
patient_groups= c('all','precentral','postcentral','insular','frontal')

#set the alpha level for tests
alpha = 0.05

###########################################------------ load the matrices and calculate the difference between pathological and healthy hemispheres--------------------------------#######################
for (subj in patient_info$Subj) {
  healthy_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_c',matrices_suffix)
  healthy_matrices =read.csv(healthy_matrices_path,header = F)
  pathological_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_p',matrices_suffix)
  pathological_matrices =read.csv(pathological_matrices_path,header = F)
  #compute edge differences and build the difference matrices
  diff_name = paste0(matrices_prefix,subj,'_diff')
  #assign variable name for calculated difference matrices
  assign(diff_name,healthy_matrices-pathological_matrices)
}
#####################################-------------------------------------------------------calculation of difference--------------------------------------------------------------------#######################
#####these loops can be ran individually
#####group-wise calculation 
for (subgroup in patient_groups){
  ######################all
  if (subgroup =='all'){
    data_loop<- patient_info %>% filter(all_patients==1)
    tfnbs_p =read.csv('group all tfnbss relut',header = F,skip=1)
  }
  ######################precentral
  else if (subgroup =='precentral'){
    data_loop<- patient_info %>% filter(precentral==1)
    tfnbs_p =read.csv('group prencentral tfnbss result',header = F,skip=1)
  }
  ######################postcentral
  else if (subgroup =='postcentral'){
    data_loop<- patient_info %>% filter(postcentral==1)
    tfnbs_p =read.csv('group postcentral tfnbss result',header = F,skip=1)
  }
  ######################insular
  else if (subgroup =='insular'){
    data_loop<- patient_info %>% filter(insular==1)
    tfnbs_p =read.csv('group insular tfnbss result',header = F,skip=1)
  }
  ######################frontal
  else if (subgroup =='frontal'){
    data_loop<- patient_info %>% filter(frontal==1)
    tfnbs_p =read.csv('group frontal tfnbss result',header = F,skip=1)
  }
  else {print('This subgroup doesnt exist.');break}
  
  #######################Only take the significant resulsts
  tfnbs_p_mask=tfnbs_p>(1-alpha)
  tfnbs_p_mask[upper.tri(tfnbs_p_mask)]<-FALSE
  n_correlation = sum(tfnbs_p_mask)
  #################################################################
  ###load empty data cell
  diff_edge=c()
  diff_motor_edge_p = c()
  diff_who_edge_p = c()
  diff_nihss_edge_p =c()
  diff_rmt_edge_p = c()
  #go through all rows of nodes
  for (i in 1:40){
    #and then go through all columns of nodes
    for (e in 1:40){
      #get into the tfnbs significant dege
      if (tfnbs_p_mask[i,e]){
        #then run this loop to this edge across the subjects
        for (subj in data_loop $Subj) {
          diff_name = paste0('S',subj,'_diff')
          #assign variable name for calculated difference matrices
          diff_edge=c(diff_edge,get(diff_name)[i,e])
        }
        
        #test for correlations for each edge and subject in the loop
        cor_motorstatus = cor.test(diff_edge,data_loop$`Motor status`,method = "spearman",exact = F)
        cor_rmt = cor.test(diff_edge,data_loop$`RMT ratio`,method = "spearman",exact = F)
        cor_who = cor.test(diff_edge,data_loop$Histology,method = "spearman",exact = F)
        cor_nihss = cor.test(diff_edge,data_loop$`NIHSS(0-42)`,method = "spearman",exact = F)
        
        #calculate the p values
        diff_motor_edge_p  = c(diff_motor_edge_p,cor_motorstatus$p.value)
        diff_who_edge_p  = c(diff_who_edge_p,cor_who$p.value)
        diff_nihss_edge_p  = c(diff_nihss_edge_p,cor_nihss$p.value)
        diff_rmt_edge_p  = c(diff_rmt_edge_p,cor_rmt$p.value)
        
        #then empty the diff_edge file
        diff_edge=c()
      }
    }
  }
  
  #after calculating the correlations, adjust with fdr 
  diff_nihss_edge_p_adjust= p.adjust(diff_nihss_edge_p,method = "fdr")
  diff_who_edge_p_adjust= p.adjust(diff_who_edge_p,method = "fdr")
  diff_motor_edge_p_adjust= p.adjust(diff_motor_edge_p,method = "fdr")
  diff_rmt_edge_p_adjust= p.adjust(diff_rmt_edge_p,method = "fdr")
  
  ##set adjusted p values index
  #this shows only the adjusted p values in the plot
  p_adj_index = 0
  for (i in 1:40){
    for (e in 1:40){
      if (tfnbs_p_mask[i,e]){
        for (subj in data_loop $Subj) {
          diff_name = paste0('S',subj,'_diff')
          #assign variable name for calculated difference matrices
          diff_edge=c(diff_edge,get(diff_name)[i,e])
        }
        #run the correlation again because we only want to plot the edges with significant correlations
        cor_motorstatus = cor.test(diff_edge,data_loop$`Motor status`,method = "spearman",exact = F)
        cor_rmt = cor.test(diff_edge,data_loop$`RMT ratio`,method = "spearman",exact = F)
        cor_who = cor.test(diff_edge,data_loop$Histology,method = "spearman",exact = F)
        cor_nihss = cor.test(diff_edge,data_loop$`NIHSS(0-42)`,method = "spearman",exact = F)
        
        
        data_loop$edge_temp=diff_edge
        
        p_adj_index = p_adj_index + 1
        #set the x coordinate of p value
        x_position = 0.8*min(data_loop$edge_temp)+0.2*max(data_loop$edge_temp)
        if (cor_motorstatus$p.value<0.05) {
          plot_motor <- ggplot(data_loop,aes(x = edge_temp,y = `Motor status`))+
            xlab('differences of n streamlines')+ylab('motor status')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_motorstatus$estimate), "p = ",sprintf("%.4f", diff_motor_edge_p_adjust[p_adj_index])))
          ##save images
          ggsave(paste0(thr,subgroup,'_motor_diff_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_motor,width = 9,height=9)
        }
        if (cor_who$p.value<0.05) {
          plot_who <- ggplot(data_loop,aes(x = edge_temp,y = Histology))+
            xlab('differences of n streamlines')+ylab('WHO°')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_who$estimate), "p = ",sprintf("%.4f", diff_who_edge_p_adjust[p_adj_index])))
          
          ggsave(paste0(thr,subgroup,'_who_diff_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_who,width = 9,height=9)
          
        }
        if (cor_nihss$p.value<0.05) {
          plot_nihss<- ggplot(data_loop,aes(x = edge_temp,y = `NIHSS(0-42)`,))+
            xlab('differences of n streamlines')+ylab('NIHSS')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_nihss$estimate), "p = ",sprintf("%.4f", diff_nihss_edge_p_adjust[p_adj_index])))
          ggsave(paste0(thr,subgroup,'_nihss_diff_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_nihss,width = 9,height=9)
          
        }
        if (cor_rmt$p.value<0.05) {
          plot_rmt <- ggplot(data_loop,aes(x = edge_temp,y = `RMT ratio`))+
            xlab('differences of n streamlines')+ylab('RMT ratio')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_rmt$estimate), "p = ",sprintf("%.4f", diff_rmt_edge_p_adjust[p_adj_index])))
          ggsave(paste0(thr,subgroup,'_rmt_diff_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_rmt,width = 9,height=9)
        }
        
        diff_edge=c()
      }
    }
  }

  ###########################------------------------------------------------calculation of pathological side---------------------------------------------------------------------#######################################
  patho_edge=c()
  patho_motor_edge_p = c()
  patho_who_edge_p = c()
  patho_nihss_edge_p =c()
  patho_rmt_edge_p = c()
  for (i in 1:40){
    for (e in 1:40){
      if (tfnbs_p_mask[i,e]){
        for (subj in data_loop$Subj) {
          pathological_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_p',matrices_suffix)
          pathological_matrices =read.csv(pathological_matrices_path,header = F)
          patho_edge=c(patho_edge,pathological_matrices[i,e])
        }
        cor_motorstatus = cor.test(patho_edge,data_loop$`Motor status`,method = "spearman",exact = F)
        cor_rmt = cor.test(patho_edge,data_loop$`RMT ratio`,method = "spearman",exact = F)
        cor_who = cor.test(patho_edge,data_loop$Histology,method = "spearman",exact = F)
        cor_nihss = cor.test(patho_edge,data_loop$`NIHSS(0-42)`,method = "spearman",exact = F)
        patho_motor_edge_p = c(patho_motor_edge_p,cor_motorstatus$p.value)
        patho_who_edge_p = c(patho_who_edge_p,cor_who$p.value)
        patho_nihss_edge_p = c(patho_nihss_edge_p,cor_nihss$p.value)
        patho_rmt_edge_p = c(patho_rmt_edge_p,cor_rmt$p.value)
        patho_edge=c()
      }
    }
  }
  patho_nihss_edge_p_adjust= p.adjust(patho_nihss_edge_p,method = "fdr")
  patho_who_edge_p_adjust= p.adjust(patho_who_edge_p,method = "fdr")
  patho_motor_edge_p_adjust= p.adjust(patho_motor_edge_p,method = "fdr")
  patho_rmt_edge_p_adjust= p.adjust(patho_rmt_edge_p,method = "fdr")
  
  p_adj_index = 0
  
  for (i in 1:40){
    for (e in 1:40){
      if (tfnbs_p_mask[i,e]){
        for (subj in data_loop$Subj) {
          pathological_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_p',matrices_suffix)
          pathological_matrices =read.csv(pathological_matrices_path,header = F)
          patho_edge=c(patho_edge,pathological_matrices[i,e])
        }      
        cor_motorstatus = cor.test(patho_edge,data_loop$`Motor status`,method = "spearman",exact = F)
        cor_rmt = cor.test(patho_edge,data_loop$`RMT ratio`,method = "spearman",exact = F)
        cor_who = cor.test(patho_edge,data_loop$Histology,method = "spearman",exact = F)
        cor_nihss = cor.test(patho_edge,data_loop$`NIHSS(0-42)`,method = "spearman",exact = F)
        
        
        
        data_loop$edge_temp=patho_edge
        p_adj_index = p_adj_index + 1
        
        #set the x coordinate of p value
        x_position = 0.8*min(data_loop$edge_temp)+0.2*max(data_loop$edge_temp)
        
        if (cor_motorstatus$p.value<0.05) {
          plot_motor <- ggplot(data_loop,aes(x = edge_temp,y = data_loop$`Motor status`))+
            xlab('n streamlines')+ylab('motor status')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_motorstatus$estimate), "p = ",sprintf("%.4f", patho_motor_edge_p_adjust[p_adj_index])))
          
          ggsave(paste0(thr,subgroup,'_motor_patho_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_motor,width = 9,height=9)
        }
        
        if (cor_who$p.value<0.05) {
          plot_who <- ggplot(data_loop,aes(x = edge_temp,y = data_loop$Histology))+
            xlab('n streamlines')+ylab('WHO°')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_who$estimate), "p = ",sprintf("%.4f", patho_who_edge_p_adjust[p_adj_index])))
          ggsave(paste0(thr,subgroup,'_who_patho_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_who,width = 9,height=9)
        }
        if (cor_nihss$p.value<0.05) {
          plot_nihss<- ggplot(data_loop,aes(x = edge_temp,y = data_loop$`NIHSS(0-42)`,))+
            xlab('n streamlines')+ylab('NIHSS')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_nihss$estimate), "p = ",sprintf("%.4f", patho_nihss_edge_p_adjust[p_adj_index])))
          
          ggsave(paste0(thr,subgroup,'_nihss_patho_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_nihss,width = 9,height=9)
        }
        if (cor_rmt$p.value<0.05) {
          plot_rmt <- ggplot(data_loop,aes(x = edge_temp,y = data_loop$`RMT ratio`))+
            xlab('n streamlines')+ylab('RMT ratio')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_rmt$estimate), "p = ",sprintf("%.4f", patho_rmt_edge_p_adjust[p_adj_index])))
          
          ggsave(paste0(thr,subgroup,'_rmt_patho_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_rmt,width = 9,height=9)
        }
        
        patho_edge=c()
      }
    }
  }

  ############---------------------------------------------------------calculation of healthy side -------------------------------------------------------------------------------------------#########################
  
  #install.packages("labelled")
  healthy_edge=c()
  healthy_motor_edge_p = c()
  healthy_who_edge_p = c()
  healthy_nihss_edge_p =c()
  healthy_rmt_edge_p = c()
  for (i in 1:40){
    for (e in 1:40){
      if (tfnbs_p_mask[i,e]){
        for (subj in data_loop $Subj) {
          healthy_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_c',matrices_suffix)
          healthy_matrices =read.csv(healthy_matrices_path,header = F)
          healthy_edge=c(healthy_edge,healthy_matrices[i,e])
        }
        cor_motorstatus = cor.test(healthy_edge,data_loop$`Motor status`,method = "spearman",exact = F)
        cor_rmt = cor.test(healthy_edge,data_loop$`RMT ratio`,method = "spearman",exact = F)
        cor_who = cor.test(healthy_edge,data_loop$Histology,method = "spearman",exact = F)
        cor_nihss = cor.test(healthy_edge,data_loop$`NIHSS(0-42)`,method = "spearman",exact = F)
        healthy_motor_edge_p = c(healthy_motor_edge_p,cor_motorstatus$p.value)
        healthy_who_edge_p = c(healthy_who_edge_p,cor_who$p.value)
        healthy_nihss_edge_p = c(healthy_nihss_edge_p,cor_nihss$p.value)
        healthy_rmt_edge_p = c(healthy_rmt_edge_p,cor_rmt$p.value)
        healthy_edge=c()
      }
    }
  }
  
  healthy_nihss_edge_p_adjust= p.adjust(healthy_nihss_edge_p,method = "fdr")
  healthy_who_edge_p_adjust= p.adjust(healthy_who_edge_p,method = "fdr")
  healthy_motor_edge_p_adjust= p.adjust(healthy_motor_edge_p,method = "fdr")
  healthy_rmt_edge_p_adjust= p.adjust(healthy_rmt_edge_p,method = "fdr")
  
  p_adj_index = 0      
  
  for (i in 1:40){
    for (e in 1:40){
      if (tfnbs_p_mask[i,e]){
        for (subj in data_loop $Subj) {
          healthy_matrices_path = paste0(matrices_path,matrices_prefix,subj,'_c',matrices_suffix)
          healthy_matrices =read.csv(healthy_matrices_path,header = F)
          healthy_edge=c(healthy_edge,healthy_matrices[i,e])
        }
        
        cor_motorstatus = cor.test(healthy_edge,data_loop$`Motor status`,method = "spearman",exact = F)
        cor_rmt = cor.test(healthy_edge,data_loop$`RMT ratio`,method = "spearman",exact = F)
        cor_who = cor.test(healthy_edge,data_loop$Histology,method = "spearman",exact = F)
        cor_nihss = cor.test(healthy_edge,data_loop$`NIHSS(0-42)`,method = "spearman",exact = F)
        
        data_loop$edge_temp=healthy_edge
        p_adj_index = p_adj_index + 1
        #set the x coordinate of p value
        x_position = 0.8*min(data_loop$edge_temp)+0.2*max(data_loop$edge_temp)
        
        if (cor_motorstatus$p.value<0.05) {
          plot_motor <- ggplot(data_loop,aes(x = edge_temp,y = data_loop$`Motor status`))+
            xlab('n streamlines')+ylab('motor status')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_motorstatus$estimate), "p = ",sprintf("%.4f", healthy_motor_edge_p_adjust[p_adj_index])))
          
          ggsave(paste0(thr,subgroup,'_motor_healthy_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_motor,width = 9,height=9)
        }
        
        if (cor_who$p.value<0.05) {
          plot_who <- ggplot(data_loop,aes(x = edge_temp,y = data_loop$Histology))+
            xlab('n streamlines')+ylab('WHO°')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_who$estimate), "p = ",sprintf("%.4f", healthy_who_edge_p_adjust[p_adj_index])))
          
          ggsave(paste0(thr,subgroup,'_who_healthy_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_who,width = 9,height=9)
        }
        if (cor_nihss$p.value<0.05) {
          plot_nihss<- ggplot(data_loop,aes(x = edge_temp,y = data_loop$`NIHSS(0-42)`,))+
            xlab('n streamlines')+ylab('NIHSS')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_nihss$estimate), "p = ",sprintf("%.4f", healthy_nihss_edge_p_adjust[p_adj_index])))
          
          ggsave(paste0(thr,subgroup,'_nihss_healthy_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_nihss,width = 9,height=9)
        }
        if (cor_rmt$p.value<0.05) {
          plot_rmt <- ggplot(data_loop,aes(x = edge_temp,y = data_loop$`RMT ratio`))+
            xlab('of n streamlines')+ylab('RMT ratio')+
            geom_point(size=3) +geom_smooth(method='glm',alpha=I(0.2))+labs(title=paste0(name_nodes$name[i]," - ",name_nodes$name[e]))+
            theme_classic(base_size = 15)+
            annotate("text",x=x_position,y=5.5,size=6,label=paste("rs = ",sprintf("%.2f",cor_rmt$estimate), "p = ",sprintf("%.4f", healthy_rmt_edge_p_adjust[p_adj_index])))
          
          ggsave(paste0(thr,subgroup,'_rmt_healthy_',name_nodes$name[i]," - ",name_nodes$name[e],".png"),plot = plot_rmt,width = 9,height=9)
        }
        healthy_edge=c()
      }
    }
  }
}
