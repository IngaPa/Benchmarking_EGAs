# script that creates supplementary table for PhD with following statistics: 
# supplement table
# 1. size enh 
# 2. size EGA
# 3. size Genes
# 4. N enhancers per gene
# 5. N genes per enahncer 
# 6. width of enhancers
# 7. N of repressed genes 
# ?8. overlap with TP
# ?9. vista enhancers?

#---------------------------------------------
# input


names <- c("JEME","GeneHancer","EnahncerAtlas","FOCS")

path.GInteractions = c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/JEME/downloaded/fANTOM_encoderoadmap_lasso_EN_pooled_processed_GR_forceBYname_19_06_19.rds",
                       "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GeneHancer/Fishilevich_2017_GeneHancer_hg19_GInteractions_190606.rds",
                       "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/EnhancerAtlas/Processed/GInteractions_EnhancerAtlas_ENSG_addedd_20918_forceByName180930.rds",
                       "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/Hait_2017_FOCS/processed/19_06_17_FOCS_pooled_processed_interactions.rds")


# corresponding mnemonics objects
path.Mnemonics <- c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Mnemomics/Enhancers/FullfANTOM_encoderoadmap_lasso_EN_pooled_processed_GR_forceBYname_19_06_19RegionsMnemonics_19_01_23.rds",
                    "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Mnemomics/Enhancers/FullFishilevich_2017_GeneHancer_hg19_GInteractions_190606RegionsMnemonics_19_01_23.rds",
                    "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Mnemomics/Enhancers/FullGInteractions_EnhancerAtlas_ENSG_addedd_20918_forceByName180930RegionsMnemonics_19_01_23.rds",
                    "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Mnemomics/Enhancers/Full19_06_17_FOCS_pooled_processed_interactionsRegionsMnemonics_19_01_23.rds")

names(path.GInteractions) <- names(path.Mnemonics) <- names





#---------------------------------------------
# functions  

.libPaths("/home/ipatarc/R/x86_64-redhat-linux-gnu-library/3.5")
# get basic statistics
basicStatistics <- function(data){
  
  # returns 
  # 1. size enh 
  # 2. size EGA  
  # 3. size Genes  
  # 4. e per gene  
  # 5. gene per e 
  # 6. N of repressed genes 
  # 7. overlap with TP
  # vista enhancers?
  
  genes <- as.character(anchorTwo(data))
  
  EperG <-  (tapply(anchorOne(data),genes,function(x){length(unique(x))}))
  
  stat.EperG <- paste0(c("range: ",min(EperG),"-",max(EperG),
                         " [med: ",median(EperG), "; mean: ",round(mean(EperG),1),"]"),
                       collapse = "")
  
  # genes per enhancer
  
  
  enhancer <- as.character(anchorOne(data))
  
  GperE <-  (tapply(genes,enhancer,function(x){length(unique(x))}))
  
  stat.GperE <- paste0(c("range: ",min(GperE),"-",max(GperE),
                         " [med: ",median(GperE), "; mean: ",round(mean(GperE),1),"]"),
                       collapse = "")
  
  enhancers <-length(unique(enhancer))
  
  genes <- length(unique(genes))
  
  EGA <- length(data)
  
  
  
  # enh width
  
  enhancer.length <- width(anchorOne(data))
  
  stat.enhancer.length <- paste0(c("range: ",min(enhancer.length),"-",max(enhancer.length),
                                   " [med: ",round(median(enhancer.length),0), "; mean: ",
                                   round(mean(enhancer.length),1),"]"),
                                 collapse = "")
  
  statistics <- c(enhancers,genes,EGA,stat.EperG,stat.GperE,stat.enhancer.length)
  names(statistics)<- c("N[enhancers]",
                        "N[genes]",
                        "N[EGA]",
                        "N[enh/gene]",
                        "N[gene/enh]",
                        "enh.length")
  return(statistics)
  
}
# get statistics about repressed regions
find_repressed_regions <- function(Mnem_FOR_Enh){
  
  require(GenomicInteractions)
  require(stringr)
  require(pheatmap)
  source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/BenchmarkGInteractions.R")
  source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/projectHelpFunctions.R")
  
  EnhRegions <- GRanges(str_replace(str_replace(rownames(Mnem_FOR_Enh),"\\.",":"),"\\.","-"))
  mcols(EnhRegions) <- DataFrame(Mnem_FOR_Enh)
  
  
  # binarize mnemonics marks, count and report in how many cell types this region is repressed 
  perEnhancerRepressedBinary <- binarizeFunction(EnhRegions,
                                                 categories= c("12_EnhBiv","11_BivFlnk","10_TssBiv","8_ZNF/Rpts","7_Enh","6_EnhG","5_TxWk","4_Tx","2_TssAFlnk","3_TxFlnk","1_TssA" ))
  
  
  
  RepressedEnhRegions <- EnhRegions[which(perEnhancerRepressedBinary==127)]
  
  
  return(length(RepressedEnhRegions))
  
  
}

# basic + repressed statistics
full_Dataset_Statistics <- function(path.GInteractions,
                                    path.Mnemonics){
  
  lapply(names(path.GInteractions), function(x){
    
    
    data <- readRDS(path.GInteractions[x])
    Mnem_FOR_Enh <- readRDS(path.Mnemonics[x])
    
    basic.stat <- basicStatistics(data)
    
    N_repressed_regions <- find_repressed_regions(Mnem_FOR_Enh)
    
    
    return(c(basic.stat,N_repressed_regions=N_repressed_regions))
    
  })
  
}


#---------------------------------------------
# analyis  

# Final_PhD_DF <- full_Dataset_Statistics(path.GInteractions,
#                                         path.Mnemonics)
# 
#saveRDS(Final_PhD_DF,
"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SupplementaryFiles/SupplementaryTable1_General_Input_Statistics.rds")
#



Final_PhD_DF <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SupplementaryFiles/SupplementaryTable1_General_Input_Statistics.rds")

#1. reorganiaztion of data
tmp <- rownames(Final_PhD_DF)
names <- c("JEME","GeneHancer","EnhancerAtlas","FOCS")
Final_PhD_DF <- Final_PhD_DF[names(Final_PhD_DF)%in%names]

Final_PhD_DF <- do.call("cbind.data.frame",Final_PhD_DF)

#2. add expected numbers (numbers reported in publications)


Expected_Enhancers <- c("92,603","489,581 Roadmap","284,834",rep(NA,4))
Expected_EGA <- c("167,988","177,532 FANTOM","1,019,746","2,488,394",rep(NA,4))


names(Expected_Enhancers) <- names(Expected_EGA) <- c("FOCS","JEME","GeneHancer","EnhancerAtlas")


# cp for ega
Expected_Enhancers <-  Expected_Enhancers[colnames(Final_PhD_DF)]
names(Expected_Enhancers)[1:3] <- colnames(Final_PhD_DF)[1:3]
t_df <- t(Final_PhD_DF)
colnames(t_df) <- tmp
finaldf <-  cbind(t_df,as.data.frame(Expected_Enhancers))


# 3. change names   
colnames(finaldf)[colnames(finaldf)=="N_repressed_regions"] <- "N[repressed regions]"
colnames(finaldf)[colnames(finaldf)=="enh.length"] <- "Enhancer length"
colnames(finaldf)[colnames(finaldf)=="Expected_Enhancers"] <- "N[enhancer, expected]"

rownames(finaldf)[rownames(finaldf)=="inhouseModels"] <- "InhouseM"
rownames(finaldf)[rownames(finaldf)=="flexible"] <- "FlexibleC"
rownames(finaldf)[rownames(finaldf)=="stringent"] <- "StringentC"


# 4.plot



png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/SupplementaryTable1_General_Input_Statistics_1.png",
    height=100, width=700)
library(gridExtra)
p<-tableGrob(finaldf[,1:5])
grid.arrange(p)
dev.off()

png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/SupplementaryTable1_General_Input_Statistics_2.png",
    height=100, width=650)
library(gridExtra)
p<-tableGrob(finaldf[,6:8])
grid.arrange(p)
dev.off()
