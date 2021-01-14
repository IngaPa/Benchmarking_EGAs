#  
# 4) check how analysis changes if SNPs in LD included
# IDEA: Analyze rs10411210 using different sources of enhancer-gene assocations


## Functions

#```{r setup, include=FALSE}
.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.5/")
library(GenomicInteractions)
library(stringr)

source("/data/akalin/Projects/AAkalin_SNP2gene/19_12_04_SNP2gene_Import_functions.R")
source("/data/akalin/Projects/AAkalin_SNP2gene/19_12_04_SNP2gene_Overlap_functions.R")
source("/data/akalin/Projects/AAkalin_SNP2gene/19_12_03_TFBS_finder.R")
source("/data/akalin/Projects/AAkalin_SNP2gene/TFBS_finder.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/plotF.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/BenchmarkGInteractions.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/downstreamFunctions.R")
#```

## Data import 

#```{r }

DISGENET <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/DISGENET/DGN_all_gene_disease_associations_noGWASC.rds")
gwascatalog <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_rtracklayer_19_02_05_hg19_NearestGeneAdded.rds")
# annotated GWAS Catalog; annotated by different EGA models
annotatedGWASC <- importAnnotatedGWAS(path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/GWAS_Catalog/annotatedGWAS_EGA/")
# repository of EGA interactions
EGA.GI <- importEGA(path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/")

# enahncer DB
endb_enhancer <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/ENdb/eadb_enhancers_PromotersGENCODE_added_GI.rds")

SNP <- "rs10411210"
#```


## Analysis 1: get SNPs in LD

SNPinLD <- snpsInLD(rsID="rs10411210",
                    pop="CEU",
                    R.squared=0.8)

SNPsinLD <- unique(SNPinLD$ID)


###  Which EGA are revealed for this SNP?

#```{r}

# 
s2ega <- lapply(SNPsinLD,snp2EGA,EGA.GI)

adjustingF <- function(s2ega,
                       SNPsinLD){

      names(s2ega) <- SNPsinLD
      
      s2ega <- s2ega[sapply(s2ega,length)!=0]
      
      s2ega.adj <- lapply(names(s2ega),function(x){
        #x="rs11880141"
        
        if (typeof(s2ega[[x]])=="list"){ rmp <- do.call("unlist",s2ega[[x]])}
        if (typeof(s2ega[[x]])!="list"){ rmp <- s2ega[[x]]}
        
          rmp$rs <- x
      
        return(rmp)
        })
      
      s2ega.df <- do.call("c",s2ega.adj)

      }  


df <- adjustingF(s2ega,SNPsinLD)
g2t_table <- ega.df(df)

Enhancer <- paste0(g2t_table[,1],":",g2t_table[,2],"-",g2t_table[,3])
EnhancerWidth <- g2t_table[,4]

DF <- data.frame(cbind(Enhancer,EnhancerWidth,
                       Gene=as.character(g2t_table[,8]),
                       Method=as.character(g2t_table[,9]),
                       SNP=as.character(g2t_table[,10])))


additionalResults <- tapply(DF$Gene,DF$SNP,unique)

DF_rs_genes <- data.frame(rbind(names(additionalResults),
                     Genes=sapply(additionalResults,paste0,collapse=", ")))

DF <- splitTable(DF_rs_genes,60)
colnames(DF) <- c("SNPs","Genes")




# SNPs and corresponding genes
png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SuppementaryTables/SupplementaryTable10411210_LDgenes.png",
    height=850, width=520)
library(gridExtra)
p<-tableGrob(DF)
grid.arrange(p)
dev.off()




####################
# What are genes that differ?

# get a list of original genes for rs10411210 and SNPs in 0.8 LD
genes <- df$X
all_Genes <- as.character(unique(genes[genes!=0]))

# get a list of original genes for rs10411210
s2g <- snp2EGA(snp_ID = "rs10411210",EGA.GI)
gens <- unique(s2g$X)


# how many new genes is detected

all_Genes[!all_Genes%in%gens]









