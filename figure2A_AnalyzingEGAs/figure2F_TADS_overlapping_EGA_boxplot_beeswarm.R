


#########################################
# STEP 1: Calculating overlap with TADs
.libPaths("/home/ipatarc/R/x86_64-redhat-linux-gnu-library/3.5")
library(GenomicInteractions)
library(stringr)
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/projectHelpFunctions.R")

pathsData <-  c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/FOCS.rds",
                "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/GeneHancer.rds",
                "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/JEME.rds",
                "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/EnhancerAtlas.rds",
                "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/inhouse_models.rds",
                "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/flexible.rds",
                "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/stringent.rds")





DIXON=list.files("/data/akalin/Base/AccessoryData/Dixon_2012_TADs/","txt", full.names = T)




withinTADsOverlap <- function(pathTADsPerCell,
                              pathData){
  
  
  TADsPerCell <- bed_to_granges(pathTADsPerCell)
  Data <- readRDS(pathData)
  
  enhTADsOverlap <- as.data.frame(findOverlaps(TADsPerCell,anchorOne(Data)))    
  geneTADsOverlap <- as.data.frame(findOverlaps(TADsPerCell,anchorTwo(Data)))    
  
  withinTAD <- sum(duplicated(rbind(enhTADsOverlap,geneTADsOverlap)))
  
  return(withinTAD/length(Data))
  
}


perDataSetStatisitcs <- lapply(pathsData, 
                               function(x,DIXON){
                                 
                                 print(x)
                                 
                                 perDataSet <- sapply(DIXON, withinTADsOverlap,pathData=x)
                                 
                                 
                                 names(perDataSet) <- basename(names(perDataSet))
                                 
                                 return(perDataSet)
                                 
                                 
                               },DIXON=DIXON)


names(perDataSetStatisitcs) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas","InhouseM","FlexibleC","StringentC")

saveRDS(perDataSetStatisitcs,"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SupplementaryFiles/19_08_05_Dixon_2012_TADS_DperDataset_Overlap_Statistics.rds")








#########################################
# STEP 2: Plotting overlap with TADs

library(reshape)
library('beeswarm')



perDataSetStatisitcs <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SupplementaryFiles/19_08_05_Dixon_2012_TADS_DperDataset_Overlap_Statistics.rds")

#names(colorScheme3) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas","InhouseM","FlexibleC","StringentC")


tst <- t(do.call("rbind.data.frame",perDataSetStatisitcs))
colnames(tst) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas","InhouseM","FlexibleC","StringentC")
rownames(tst) <- 1:nrow(tst)

#tst$X2 <- as.character(tst$X2)


colnames(tst) <- names(perDataSetStatisitcs)
rownames(tst) <- 1:nrow(tst)
tst <- tst[,1:4]
tst <- melt(tst)
tst$X2 = as.character(tst$X2)


#colorScheme3 <- c("black","darkorange","blue","darkcyan","orange","red","darkred")
#names(colorScheme3) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas","InhouseM","FlexibleC","StringentC")

colorScheme3 <- c("black","darkorange","blue","darkcyan")
names(colorScheme3) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")




##############################
# beeswarm + boxplot


png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure4_benchmarking/Figure5D_boxplot_perMethodsTADS.png",
    width = 600,height = 550 )


      par(mar=c(20,8,3.4,2.5))
      
      par(font.axis = 2)
      par(font.lab=2)
      
      
      boxplot(value ~ X2, data = tst, 
              outline = FALSE,     ## avoid double-plotting outliers, if any
              main = '',
              lwd=2,
              las=2,
              cex.lab=2.4,
              cex.axis=3,
              xaxt = "n",
              yaxt="n",
              ylab="% of EGA within TADs \n",
              xlab="",srt = 60)
      
      text(x = c(0.5,1.8,2.5,3.8),
           y = c(0.43,0.49,0.45,0.49),
           labels = c("EnhancerAtlas","FOCS","GeneHancer","JEME"),
           xpd = NA,
           ## Rotate the labels by 35 degrees.
           srt = 45,
           cex = 3, font.lab=2 )

        beeswarm(value ~ X2, data = tst,   
               col = colorScheme3[order(unique(tst$X2))], pch = 16, add = TRUE)
      
     # axis(1,cex.axis=3)
       axis(2, cex.axis=3)
      
      dev.off()









      
      
      
##############################
# supplementary table 
      
      
      
      TADS <- lapply(perDataSetStatisitcs, summary)
      
      tads_stat <- do.call("rbind.data.frame",TADS)
      
      rownames(tads_stat) <- names(TADS)
      colnames(tads_stat) <- names(TADS[[1]])
      tads_stat <- round(tads_stat,3)
      
      
      png.path="/data/akalin/Projects/AAkalin_Catalog_RI/AAkalin_CatalogRI/Results/Plots/PAPER/"
          png(paste0(png.path,"SupplTable6_TADsCoverage.png") ,
              height=220, width=415)
          library(gridExtra)
          p<-tableGrob(tads_stat)
          grid.arrange(p)
          dev.off()


