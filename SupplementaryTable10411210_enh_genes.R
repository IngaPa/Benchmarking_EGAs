# PLOT pheatmap for per gene per method statistics

.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.5/")
library(GenomicInteractions)
library(stringr)

source("/data/akalin/Projects/AAkalin_SNP2gene/19_12_04_SNP2gene_Import_functions.R")
source("/data/akalin/Projects/AAkalin_SNP2gene/19_12_04_SNP2gene_Overlap_functions.R")
source("/data/akalin/Projects/AAkalin_SNP2gene/19_12_03_TFBS_finder.R")
source("/data/akalin/Projects/AAkalin_SNP2gene/TFBS_finder.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/plotF.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/BenchmarkGInteractions.R")


splitTable <- 
function(dataFrame,
         everyN=30,
         ...){
  
  
  dataFrame <- apply(dataFrame,1,as.character)
  
  addNewlines <- function(string){
    
    d <- 1:nchar(string)
    # split string into chunks of 100 pieceses and ad newline
    Splitted <- split(d, ceiling(seq_along(d)/everyN))
    String <- unlist(str_split(string,""))
    
    newline_list <- sapply(Splitted, function(x){
      #print(x)
      return(paste(paste0(String[x],collapse = ""),"\n"))})
    
    return(paste0(newline_list,collapse = ""))
    
    
    
  }
  
  ST1.a <- dataFrame  
  
  for (i in 1:ncol(dataFrame)){
    for (j in 1:nrow(dataFrame)){
      print(i)
      print(j)
      if (!is.na(dataFrame[j,i])){ST1.a[j,i] <- addNewlines(dataFrame[j,i]) }
      
    }
  }
  
  return(ST1.a)
  
}


## Data import 


gwascatalog <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_rtracklayer_19_02_05_hg19_NearestGeneAdded.rds")
# annotated GWAS Catalog; annotated by different EGA models
annotatedGWASC <- importAnnotatedGWAS(path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/GWAS_Catalog/annotatedGWAS_EGA/")

GWAS_FOCS <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/FOCSBenchVal2BenchVal.rds")

# enahncer DB
``
SNP <- "rs10411210"


## Analysis 1

###Which genes are annotated for this SNP?

#```{r}
# 
s2g <- snp2gene(SNP,
                annotatedGWASCatalogs = annotatedGWASC,
                gwascatalog)
genes <- sapply(s2g,paste0,collapse=", ")

#```

### Which enhancers are annotated for this SNP?

#```{r}
s2e <- snp2enhancer(SNP,annotatedGWASC)
enhancers <- sapply(s2e,paste0,collapse=", ")


enhancers <- enhancers[names(enhancers)%in%c("EnhancerAtlas","JEME","GeneHancer","FOCS")]
genes <- genes[names(genes)%in%c("EnhancerAtlas","JEME","GeneHancer","FOCS")]


Df <- rbind(genes,enhancers)
DF <- splitTable(Df,30)
rownames(DF) <- colnames(Df)
png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/SupplementaryTable10411210_enh_genes.png",
    height=450, width=520)
library(gridExtra)
p<-tableGrob(DF)
grid.arrange(p)
dev.off()


#```


