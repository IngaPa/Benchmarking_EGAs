# plots mean of mean of enhancer activity across cells
# step 1 calculation mean(mean(enhancerActivity for one enhancer) across enhancers)
# step 2 plot
# NOTE: chr1 based analysis

library(GenomicInteractions)
library(pheatmap)


meanExpression <- function(path="/data/local/ipatarc/AAkalin_CatalogRI/Data/regActivityImputed/EnhancerAtlas_chr1/"){
  
  # for each enhancer calculated mean expression across cell types
  EA <- list.files(path,
                   "mean.rds", full.names = T)
  
  
      All_methods <- lapply(EA,readRDS)
      
      names(All_methods) <- str_extract(EA,"DNase|H3K27ac|H3K4me1|Methylation")
      
      perMethod <- lapply(All_methods,function(x){
        
        mm <- as.data.frame(mcols(x))
        
        rowMeans(mm)
        
        
      })
      
  
      return(perMethod)
  
  
}


# step 1a calculation mean(enhancerActivity for one enhancer) 

# calculate mean expression across cell types for each enahncer
  EAexp <- meanExpression(path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Publication_Enhancers_and_EGA_Analysis/EnhancerAtlas_regActivity_chr1/")
  FOCSexp <- meanExpression("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Publication_Enhancers_and_EGA_Analysis/FOCS_chr1/")
  GHexp <- meanExpression("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Publication_Enhancers_and_EGA_Analysis/GeneHancer_chr1/")
  jemeExp <- meanExpression("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Publication_Enhancers_and_EGA_Analysis/JEME_chr1/")
  

  
 
  
# step 1b calculation mean(mean(enhancerActivity for one enhancer) across enhancers)
  
  EA <- as.data.frame(unlist(EAexp))
  EA$Publication <- "EnhancerAtlas"
  EA$Method <- rep(names(EAexp),each=length(EAexp[[1]]))
  colnames(EA)[1] <- "Value"
  
  sapply(FOCSexp,length)
  
  FOCS <- as.data.frame(unlist(FOCSexp))
  FOCS$Publication <- "FOCS"
  FOCS$Method <- rep(names(FOCSexp),each=length(FOCSexp[[1]]))
  colnames(FOCS)[1] <- "Value"
  
  GH <- as.data.frame(unlist(GHexp))
  GH$Publication <- "GeneHancer"
  GH$Method <- rep(names(GHexp),each=length(GHexp[[1]]))
  colnames(GH)[1] <- "Value"
  
  
  JEME <- as.data.frame(unlist(jemeExp))
  JEME$Publication <- "JEME"
  JEME$Method <- rep(names(jemeExp),each=length(jemeExp[[1]]))
  colnames(JEME)[1] <- "Value"
  
  
  All.together <- rbind(EA,GH,FOCS,JEME)
  All.together$Publication <- as.factor(All.together$Publication)
  
  colorScheme <- c("black","darkorange","blue","darkcyan" )
  names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")
  colors <- colorScheme[levels(All.together$Publication)]
  colors <- colors[complete.cases(colors)]

  
  
  
  
  png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/fig1F_enhancerActivity_barplot_ggplot_.png",
      width = 1000,height = 500)
  
  
  ggplot(All.together,
         aes(order =Publication , 
             y = log10(Value),
             x=Method),fill=Publication) +
     theme_minimal()+
    #geom_violin(aes(fill = factor(Publication),position = "dodge",width=0.6)) + 
   scale_fill_manual(values =  alpha(rep(colors,4), .7))+
  

    geom_boxplot(aes(fill = Publication),
                     position = "dodge",
                     width=0.6, outlier.fill = NULL, outlier.size = 0.1)+
    theme(text = element_text(size=30,face="bold"),
            axis.text =element_text(size=25,face="bold"),
            legend.text = element_text(size=30,face="bold"))+
    #geom_point(aes(size = 0.0001), colour = "blue", position = "jitter",alpha=0.2) +
    xlab("\n Method") +  ylab ("log10 (enhancer activity)\n")

  
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  signal <- do.call("cbind.data.frame",lapply(list(EAexp,
                                                   FOCSexp,
                                                   GHexp,
                                                   jemeExp),function(x){
    
      # mean expression across cell types
                                                     
    do.call("rbind.data.frame",lapply(x,mean,na.rm=T))
    
    
  }))
  



    colnames(signal) <- c("EnhancerAtlas","FOCS","GeneHancer","JEME")
    rownames(signal) <- names(EAexp)



# step 2 plot 
    
    png("/data/akalin/Projects/AAkalin_Catalog_RI/AAkalin_CatalogRI/Results/Plots/PAPER/EnhancerActivity_pheatMap.png", width = 750, height = 600)


pheatmap(signal,
         display_numbers = T,
         number_format = "%.2f",  
         fontsize = 29,
         cluster_rows =  F,
         cluster_cols = F )

dev.off()






#########################
# z-score based

chr1_mean_exp <- sapply(chr1_chopped,mean,na.rm=T)
chr1_sd_exp <- sapply(chr1_chopped,sd,na.rm=T)

z_score <- apply(signal,1,function(x,chr1_mean_exp,chr1_sd_exp){
  
  (x-chr1_mean_exp)/chr1_sd_exp
  
},chr1_mean_exp=chr1_mean_exp,chr1_sd_exp=chr1_sd_exp)




png("/data/akalin/Projects/AAkalin_Catalog_RI/AAkalin_CatalogRI/Results/Plots/PAPER/EnhancerActivity_pheatMap_ascore.png", width = 750, height = 600)


pheatmap(z_score,
         display_numbers = T,
         number_format = "%.2f",  
         fontsize = 29,
         cluster_rows =  F,
         cluster_cols = F ,lwd=3)

dev.off()





