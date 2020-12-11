
.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.5/")
library(stringr)
library(pheatmap)
library(gwascat)
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/generalPlot_functions.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/projectHelpFunctions.R")
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/GWASCatalog/19_06_12_helpFunction_for_GWAS.R")

GWASCatalog="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_20_12_05_non_codingSNPs.rds"
GWAS_PATHS <- "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/GWAS_Catalog/annotatedGWAS_EGA/"
DISGENET <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/DISGENET/DGN_all_gene_disease_associations_noGWASC.rds") 


EFO.Terms <- c("EFO:0000305",
               "EFO_0001360",
               "EFO_0002690",
               "EFO_0003767",
               "EFO_0000384",
               "EFO_0004339",
               "EFO_0004340",
               "EFO:0005842",
               "EFO:0000692",
               "EFO:0001645")


Disease.Terms <- c("BC",
               "T2D",
               "SLE",
               "IBD",
               "CD",
               "Height",
               "BMI",
               "CRC",
               "SCZ",
               "CAD")

EFO2diseases <- cbind(EFO.Terms,
             Disease.Terms)

# Breast carcinoma: EFO:0000305 C0678222
# T2D EFO_0001360 CN244395
# SLE EFO_0002690 
# IBD EFO_0003767 C0021390 
# CD EFO_0000384 CN043071
# height EFO_0004339
# BMI EFO_0004340

GWASCatalog="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_20_12_05_non_codingSNPs.rds"
SNPannotations <- GWAS_PATHS


extractGWAS <- function(EFO2disease,
                        GWASCatalog,
                        SNPannotations){
  
  #EFO2disease <- EFO2diseases[x,]

  GWASCatalogGeneSets <- getGWASCatalogGeneSets(terms=c(EFO2disease[1]),
                                              GWASCatalog=GWASCatalog,
                                              SNPannotations=SNPannotations)

  #filtering out additional analyses
  GWASCatalogGeneSets <- GWASCatalogGeneSets[names(GWASCatalogGeneSets)%in%c("GeneHancer",
                                                                           "FOCS",
                                                                           "EnhancerAtlas",
                                                                           "JEME")]


  GWASCatalogGeneSets <- lapply(GWASCatalogGeneSets,
                              function(x){return((unique(x)))})

        nGenes <- sapply(GWASCatalogGeneSets,length)
        
        Overlap <- as.data.frame(cbind(Publication=names(nGenes),
                                              value=nGenes))
        
        Overlap$Disease <- EFO2disease[2]
        
        return(Overlap)
      
}  
    
    
    
    
    # lapply(1:nrow(EFO2disease), 
disease.stat <- lapply(1:nrow(EFO2diseases), 
           function(x){
             print(x)
             return(extractGWAS(EFO2diseases[x,],
                         GWASCatalog=GWASCatalog,
                         SNPannotations=SNPannotations))})
           
           
    
 stats <- do.call("rbind.data.frame",
                  disease.stat)   
  
 stats$Disease[stats$Disease=="height"] <- "Heigth"  
    
    



png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/Bench_N_GWAS_Disease_genes_10diseases.png",
    width = 1150,height = 500)


colorScheme <- c("black","darkorange","blue","darkcyan")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")


ggplot(stats,aes(x=Disease,
                 y=as.integer(as.character(value)),
                 fill=Publication,
                 order=Disease), color=colorScheme) +  
  stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
  scale_fill_manual(values =  alpha(colorScheme[sort(names(colorScheme))], .9))+
  ylab("N [genes] \n") +
  xlab("\n Diseases") + 
  #ylim(0,85) +
  theme(text = element_text(size=40,face="bold"),
        axis.text =element_text(size=40,face="bold"),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=40,angle = 45, hjust = 1),
        legend.text = element_text(size=35,face="bold"))
geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()

    
    
    
