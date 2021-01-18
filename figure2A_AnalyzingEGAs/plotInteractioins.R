####################################
# -------------------------------
# plotinteractions for rs10411210




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

png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/fiigure7_rs10411210/plotInteractions_acrsP_rs10411210.png",height = 550, width = 425)

plotInteractions(interactions=EGA.GI,
                 rangeGenes=NULL,
                 selectGene=NULL,
                 selectRegulatoryRegion=getsnplocation("rs10411210"),
                 filters="hgnc_symbol",
                 benchInteractions=NULL,
                 statistics=NULL,
                 coloring=NULL,
                 interactionsNaming="Enh~Prom Interactions",
                 benchmarkNaming="Benchmark",
                 sizes=0.3,
                 cex.title=1, 
                 plotAllAnchors1=FALSE)

dev.off()

#```

plotinteractions for rs10411210 + eaDB nothing special

#```{r}

png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/fiigure7_rs10411210/plotInteractions_acrsP_rs10411210_endb.png",height = 550, width = 425)

plotInteractions(interactions=EGA.GI,
                 rangeGenes=NULL,
                 selectGene=NULL,
                 selectRegulatoryRegion=getsnplocation("rs10411210"),
                 filters="hgnc_symbol",
                 benchInteractions=endb_enhancer,
                 statistics=NULL,
                 coloring=NULL,
                 #interactionsNaming="Enh~Prom Interactions",
                 benchmarkNaming="eaDB",
                 sizes=0.3,
                 cex.title=1, 
                 plotAllAnchors1=FALSE)
dev.off()
#```