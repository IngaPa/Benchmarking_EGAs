# get density plots for: distance to nearest TSS and width statistics

.libPaths("/home/ipatarc/R/x86_64-redhat-linux-gnu-library/3.5")
library(pheatmap)            
library(GenomicInteractions)
library(rtracklayer)
library(genomation)
library(ggplot2)
library(reshape)
library(stringr)
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/generalPlot_functions.R")


#############################################
#---------------------------------------------
# input 
#############################################




#############################################
# Size distribution PLOT
############################################



#---------------------------------------------
# input


names <- c("inhouseModels","flexible","stringent","jeme","GeneHancer","EnahncerAtlas","FOCS")

path.GInteractions = c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/JEME/downloaded/fANTOM_encoderoadmap_lasso_EN_pooled_processed_GR_forceBYname_19_06_19.rds",
                       "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GeneHancer/Fishilevich_2017_GeneHancer_hg19_GInteractions_190606.rds",
                       "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/EnhancerAtlas/Processed/GInteractions_EnhancerAtlas_ENSG_addedd_20918_forceByName180930.rds",
                       "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/Hait_2017_FOCS/processed/19_06_17_FOCS_pooled_processed_interactions.rds")


names(path.GInteractions) <- names



EperGStatistics <- function(data){
  

  genes <- as.character(anchorTwo(data))
  
  EperG <-  (tapply(anchorOne(data),genes,function(x){length(unique(x))}))
  
  return(EperG)
  
}

GperEStatistics <- function(data){
  
  
  enhancer <- as.character(anchorOne(data))
  
  EperG <-  (tapply(anchorTwo(data),enhancer,function(x){length(unique(x))}))
  
  return(EperG)
  
}


##############################################
#=------------------------------------------
# get statistics
# get statistics N genes per enhancer, and vice-versa

# GperEStatistics_results <- lapply(names(path.GInteractions), function(x){
#   
#   data <- readRDS(path.GInteractions[x])
#   
#   return(GperEStatistics(data))
#   
#   })
# 
# saveRDS(GperEStatistics_results,"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/SupplementaryFiles//GperEStatistics_results.rds")
# 
# 
# EperGStatistics_results <- lapply(names(path.GInteractions), function(x){
#   
#   data <- readRDS(path.GInteractions[x])
#   
# 
#   
#   return(EperGStatistics(data))
#   
# })
# 
# saveRDS(EperGStatistics_results,"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SupplementaryFiles/EperGStatistics_results.rds")
# 

    EperGStatistics_results <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SupplementaryFiles/EperGStatistics_results.rds")
    GperEStatistics_results <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SupplementaryFiles/GperEStatistics_results.rds")
    
    
    names(GperEStatistics_results) <- names(EperGStatistics_results) <- names


    
    
    
    
    
################################################
#--------------------------------------------
# 1c.enahncer per gene
    
    EperGStatistics_results <- EperGStatistics_results[names(EperGStatistics_results)%in%c("FOCS","GeneHancer","JEME","EnhancerAtlas")]
    
    plotData <- as.data.frame(unlist(EperGStatistics_results))
    
    # get appropriate names
    tmp <- sapply(EperGStatistics_results, length)
    names_Publ <- rep(names(tmp),tmp)
    
    # get appropriate order
    ggplotData <- cbind(plotData,names_Publ)
    names(ggplotData) <- c("value","Publication")
    
    
    # colors
    colorScheme <- c("black","darkorange","blue","darkcyan" )
    names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")
    
    
    # plot
    
    png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure2_publicationEGA/2C_Enh_per_Gene_distribution_Plot.png"
        ,width = 900,height = 900)
    
    ggplot(ggplotData, aes(value, 
                           colour = Publication,
                           fill=Publication)) +
      xlab("\n N[enhancers per gene]")+
      ylab("Density \n")+
      xlim(1,100)+
      geom_density(adjust=8, size=4)+
      scale_color_manual(values =  alpha(colorScheme[levels(ggplotData$Publication)], .7))+
      scale_fill_manual(values =  alpha(colorScheme[levels(ggplotData$Publication)], .1))+
      #theme_minimal()+
      theme(legend.position="none",
        panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title = element_text(size=50,face="bold"),
            axis.text =element_text(size=40,face="bold"),
            text = element_text(size=40,face="bold", vjust = 2),
            legend.text = element_text(size=20,face="bold"))

    
    
    
    dev.off()    
    
    
    
    



################################################
#--------------------------------------------
# 1c.genes per enhancer

    GperEStatistics_results <- GperEStatistics_results[names(GperEStatistics_results)%in%c("FOCS","GeneHancer","JEME","EnhancerAtlas")]
    
    plotData <- as.data.frame(unlist(GperEStatistics_results))
    
    # get appropriate names
    tmp <- sapply(GperEStatistics_results, length)
    names_Publ <- rep(names(tmp),tmp)
    
    # get appropriate order
    ggplotData <- cbind(plotData,names_Publ)
    names(ggplotData) <- c("value","Publication")


# colors
colorScheme <- c("black","darkorange","blue","darkcyan" )
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")


# plot

    png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure2_publicationEGA/2B_Genes_per_Enhancer_distribution_Plot.png"
        ,width = 900,height = 900)
    
    ggplot(ggplotData, aes(value, 
                           colour = Publication,
                           fill=Publication)) +
      xlim(1,10)+
      xlab("\n N[genes per enhancer]")+
      ylab("Density \n")+
      geom_density(adjust=8, size=4)+
      scale_color_manual(values =  alpha(colorScheme[levels(ggplotData$Publication)], .7))+
      scale_fill_manual(values =  alpha(colorScheme[levels(ggplotData$Publication)], .1))+
      #theme_minimal()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title = element_text(size=50,face="bold"),
            axis.text =element_text(size=40,face="bold"),
            text = element_text(size=40,face="bold", vjust = 2),
            legend.text = element_text(size=40,face="bold"))

    
    
    dev.off()

















