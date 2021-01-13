# IDEA: using Gasperini dataset define TP EGAs in the genome as well as FN - 
# this is done by identifying T enhancers, and testing which interactions were predicted or not
# in the 2nd step we test number of enhancer that overlap enhaners which activity was not identified well

# thus, I have 2 directions of analysis:
# 1) overlap succesful enhancers: and get TP EGAs, FN EGAs, and TN EGAs, FP enhaneers
 

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
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/BenchmarkGInteractions.R")


#Filtering out enhancers that could not be benchmarked at all
filterEnhancer <- function(EGA,GaspEnhancers){
  
  Enhancers <- anchorOne(EGA)
  
  Data <- as.data.frame(findOverlaps(Enhancers,
                                     GaspEnhancers))
  
  return(EGA[unique(Data$queryHits)])
}

# function that reports number of overlaps between Gasperini true enhancers and false enhancers
Enhancer_overlap <- function(EGAs,
                             GaspEnhancers){
  
  # separating positive and negative enhnacers
  PositiveEn <- GaspEnhancers[GaspEnhancers$Activity=="Y"]
  NegativeEn <- GaspEnhancers[GaspEnhancers$Activity!="Y"]
  
  #detecting overlap with positive enhancers
  Pos_overlap <- data.frame(findOverlaps(EGAs,PositiveEn))
  N_P <- length(unique(Pos_overlap$queryHits))
  
  # detetin overlap with negative ehancers  
  Neg_overlap <- data.frame(findOverlaps(EGAs,NegativeEn))
  N_N <- length(unique(Neg_overlap$queryHits))
  
  
  return(c(N_P,N_N))
  
}



# function that reports number of overlaps between Gasperini enhancer-gene assoiations and EGAs
EGAs_overlap <- function(EGAs,
                         TestEGAs=Gasperini_enh_gene,
                         report="overlappingPairs"){
  
  
  # prior runnning overlap extend +/-1000bp
  tmp <- mcols(EGAs)
  EGAs <- GInteractions(anchorOne(EGAs),
                        promoters(anchorTwo(EGAs),1000,1000))
  mcols(EGAs) <- tmp
  
  # running benchmarking!
  tmp <- benchmarkInteractions(EGAs,Gasperini_enh_gene)
  if (report=="overlappingPairs"){res <- unique(tmp[tmp$Bench!=0])}
  if (report=="non_overlappingPairs"){res <- unique(tmp[tmp$Bench==0])}
  
  return(res)
}



# IMPORT
# test enhancer dataset from Gasperini et al. 2019. Two sets of enhancers were analyzed: positive and negative interactions
# Tested enhaner regions in Gasperini
Datasets <- list.files("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/",full.names = T,pattern = "rds")
GaspEnhancers = readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Data/Gasperini/Gasperini_S2A_tested_enhancers_GRanges.rds")
# assessed positive interactions
Gasperini_enh_gene <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/Gasperini2019_Cell_CRISPR_enhGene/Gasperini2019_Cell_CRISPR_enhGene.rds")
# results of the inhouse models
EA <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/EnhancerAtlas.rds")
GH <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/GeneHancer.rds")
JEME <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/JEME.rds")
FOCS <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/FOCS.rds")

All_Tested_Models <- list(EA,FOCS,GH,JEME)
# STEP1 Filtering out enhancers that could not be benchmarked at all

All_Tested_Models.filEnh <- lapply(All_Tested_Models,
                                   filterEnhancer,
                                   GaspEnhancers)
names(All_Tested_Models.filEnh) <- c("EA","FOCS","GH","JEME")

sapply(All_Tested_Models.filEnh,length)

# there are 67151  1806 14546 13257 interactions that can be benchmarked using Gasperini Interactions 
# the next line is per method - how much it can be benchmarked in the first place
#apply(as.data.frame(mcols(All_Tested_Models.filEnh)[-c(1:2)]),2,table)



# names of the results of the modelling procedure
Modeling_results <- colnames(mcols(All_Tested_Models.filEnh))[-c(1:2)]

# get me a function of EGAsand ehancers

###############
# assessing overlap with Gaperinin TP

TP_FP <- lapply(names(All_Tested_Models.filEnh),function(x){
  
  EGA <- All_Tested_Models.filEnh[[x]]

      TP <- EGAs_overlap(EGA,
                        Gasperini_enh_gene,"overlappingPairs")
      
      FP <- EGAs_overlap(EGA,
                   Gasperini_enh_gene,
                   "non_overlappingPairs")
  return(c(length(TP),
           length(FP)))
  
})

Gasperini_TP <- do.call("rbind.data.frame",TP_FP)
colnames(Gasperini_TP) <- c("TP","FP")

#Gasperini_TP <- melt(Gasperini_TP)
Gasperini_TP$Publication <-  c("EnhancerAtlas","FOCS","GeneHancer","JEME")


# plot 1 normal plot 
path="/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/fig3_TP_FP_normal_plot.png"

colorScheme <- c("black","darkorange","blue","darkcyan")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")

png(path,width = 900, height = 600)
ggplot(data=Gasperini_TP, aes(x=TP,
                              y=FP,
                              color=Publication)) +
  geom_point(size=16) +
  scale_color_manual(values=colorScheme[Gasperini_TP$Publication])+
  scale_size_continuous(range=c(20,30)) +
  theme(legend.key = element_rect(colour = NA, fill = NA),
        plot.margin = unit(c(2,0,0,1), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=40,angle = 45, hjust = 1))+
  ylab(" N[FP]\n") +
  xlab("\n N[TP]") +
  guides( size=FALSE,
          colour=guide_legend(override.aes = list(size=25)))

dev.off()




