
.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.5/")
library(rtracklayer)
library(ggplot2)
# library(reg2gene)
library(stringr)
library(parallel)
library(GenomicRanges)
library(GenomicInteractions)
library(InteractionSet)
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/BenchmarkGInteractions.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/projectHelpFunctions.R")
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/generalPlot_functions.R")
options("scipen"=100, "digits"=4)

# 4* higher overlap with HIC but ranking stayed the same
Fullgenome <- 3094490490


paths <- c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/JEME/downloaded/fANTOM_encoderoadmap_lasso_EN_pooled_processed_GR_forceBYname_19_06_19.rds",
           "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GeneHancer/Fishilevich_2017_GeneHancer_hg19_GInteractions_190606.rds",
           "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/EnhancerAtlas/Processed/GInteractions_EnhancerAtlas_ENSG_addedd_20918_forceByName180930.rds",
           "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/Hait_2017_FOCS/processed/19_06_17_FOCS_pooled_processed_interactions.rds")


names(paths) <- c("JEME","GeneHancer","EnhancerAtlas","FOCS")



bench_path <- "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GTEx/Processed/GTX.df.txt_formatted_uniqueEPpairs_forcebyName_190702.rds"



getRepressedRegions <- function(path.mnem="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Mnemomics/Enhancers/"){
  
  
  mne <- list.files(path.mnem,full.names = T)[c(1:3,5)]
  
  
  OriginalData.path <- c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/FOCS.rds",
                    "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/JEME.rds",
                    "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/GeneHancer.rds",
                    "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/EnhancerAtlas.rds")
  
  
  names(mne) <- names(OriginalData.path) <-  c("FOCS","JEME","GeneHancer","EnhancerAtlas")
  
  
  nofRepressed <- sapply(names(mne),function(x){
    
          #x=names(mne)[1]
          Mnem_FOR_Enh <- readRDS(mne[x])
          
          # step 1. Identify repressed enhancers
          # creating Granges object
            EnhRegions <- GRanges(str_replace(str_replace(rownames(Mnem_FOR_Enh),"\\.",":"),"\\.","-"))
              mcols(EnhRegions) <- DataFrame(Mnem_FOR_Enh)
          
          
          # binarize mnemonics marks, count and report in how many cell types this region is repressed 
          perEnhancerRepressedBinary <- binarizeFunction(EnhRegions,
                                                         categories= c("12_EnhBiv","11_BivFlnk","10_TssBiv","8_ZNF/Rpts","7_Enh","6_EnhG","5_TxWk","4_Tx","2_TssAFlnk","3_TxFlnk","1_TssA" ))
          
          RepressedEnhRegions <- EnhRegions[which(perEnhancerRepressedBinary==127)]
          
          
          
          # step 2.identify EGA with repressed enhancers
          OriginalData <- readRDS(OriginalData.path[x])
          
          tmp <- DataFrame(findOverlaps(RepressedEnhRegions,anchorOne(OriginalData)))
          
          # length repressed EGA
          repressedEGA <- length(unique(tmp$subjectHits))
          
          return(repressedEGA)
  
  })
  
  
  names(nofRepressed) <- names(OriginalData.path)
  
  return(nofRepressed)
}

getRankingOrginal <- function(path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/PCHiCBenchVal2BenchVal.rds",
                              names=c("JEME","GeneHancer","EnhancerAtlas","FOCS" )){
  
  # get benchmarking statistics
  HIC.benchmarking <- readRDS(path)
  
  metadata <- as.data.frame(mcols(HIC.benchmarking))[,names]
  
  statistics <- colSums(metadata)
  
  
  # get size of benchmark datasets
  paths <- c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/JEME/downloaded/fANTOM_encoderoadmap_lasso_EN_pooled_processed_GR_forceBYname_19_06_19.rds",
             "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GeneHancer/Fishilevich_2017_GeneHancer_hg19_GInteractions_190606.rds",
             "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/EnhancerAtlas/Processed/GInteractions_EnhancerAtlas_ENSG_addedd_20918_forceByName180930.rds",
             "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/Hait_2017_FOCS/processed/19_06_17_FOCS_pooled_processed_interactions.rds")
  
  
  data <- lapply(paths,readRDS)
  names(data) <- c("JEME","GeneHancer","EnhancerAtlas","FOCS")
  size.vector <- (sapply(data, function(x){length((x))}))
  names(size.vector) <- names(data)
  size.v <- 1:length(paths)
  names(size.v) <-  names(sort(size.vector))
  
  
  df <- as.data.frame(cbind(Noverlaps=sort(statistics),
                            FullSize=size.vector[names(sort(statistics))]#,
                            #rank1=1:length(paths),
                            #rank2=size.v[names(sort(statistics,decreasing = T))]
  ))
  
  
  df$Ratio <- round(df$Noverlaps*100/df$FullSize,2)
  
  
  
  return(df)
  
}





###############################################
#--------------Analysis----------------------
###############################################

# Get TPs

Gasperini_TP <- getRankingOrginal(path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/GasperiniTPBenchVal2BenchVal.rds")


# ranking is the same, regardless of the Benchmark used

# Get FPs
nofRepressed <- getRepressedRegions(path.mnem="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Mnemomics/Enhancers/")




Gasperini_TP$FP <- nofRepressed[row.names(Gasperini_TP)]
Gasperini_TP$Publication <- row.names(Gasperini_TP)
Gasperini_TP$TP  <- Gasperini_TP$Noverlaps 

#########################################################
# ----------------------------------------------------




colorScheme <- c("black","darkorange","blue","darkcyan")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")




# plot 1 normal plot 
path="/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/fig3_TP_FP_normal_plot.png"

png(path,width = 900, height = 600)
ggplot(data=Gasperini_TP, aes(x=FP, 
                       y=TP,
                       fill=Publication,
                       size=FullSize)) +
  geom_point(aes(colour=Publication)) +
  scale_color_manual(values=colorScheme[Gasperini_TP$Publication])+
  scale_size_continuous(range=c(20,30)) +
  theme(legend.key = element_rect(colour = NA, fill = NA),
    plot.margin = unit(c(6,0,0,1), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    text = element_text(size=40,face="bold"),
    axis.text.x = element_text(size=40,angle = 45, hjust = 1))+
  ylim(0,320)+xlim(0,37000)+
  ylab(" N[true positives]\n") +
  xlab("\n N[false positives]") +
  guides( size=FALSE,
         colour=guide_legend(override.aes = list(size=25)))

dev.off()








# plot 2 zoom in plot
#path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure5_TP_FP/fig5B_TP_FP_zoomin_plot.png"

png(path,width = 400, height = 400)
ggplot(data=Gasperini_TP, aes(x=FP, 
                              y=TP,
                              fill=Publication,
                              size=FullSize)) +
  geom_point(aes(colour=Publication)) +
  scale_color_manual(values=(colorScheme[Gasperini_TP$Publication]))+
  scale_size_continuous(range=c(5,10)) +
  theme(legend.position =  "none",
    legend.key = element_rect(colour = NA, fill = NA),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=40,angle = 45, hjust = 1))+
  ylim(0,320)+xlim(0,2000)+
  ylab(" N[TP]\n") +
  xlab("\n N[FP]") +
  guides( size=FALSE,
          colour=guide_legend(override.aes = list(size=25)))

dev.off()