# Supplementary table: CRC genes associated with each method

# ---
#   title: "Annotating GWAS Catalog with different EGA methods"
# output:
#   html_document:
#   df_print: paged
# ---

#Showcase CRC:

# IDEA
Check how different sets of genes that were associated with colorectal cancer based on the overlap 
with enhancers from different EGA publications overlap 3 gene sets: 
  
  
  
  2) overlap based on genes annotated based on overlap of CRC-associated SNPS with promoter regions of identified SNPs 
3) DISGENET reported genes in disease of interest, and corresponding disease ancestors and descendents




#colorectal carcinoma: http://www.ebi.ac.uk/efo/EFO_0005842
# colorectal carcinoma: UMLS:C1527249


##```{r,echo=FALSE,message=FALSE}
.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.5/")
library(stringr)
library(pheatmap)
library(gwascat)
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/generalPlot_functions.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/projectHelpFunctions.R")
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/GWASCatalog/19_06_12_helpFunction_for_GWAS.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/downstreamFunctions.R")


GWASCatalog="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_rtracklayer_19_02_05_hg19_NearestGeneAdded.rds"

GWAS_PATHS <- "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/GWAS_Catalog/annotatedGWAS_EGA/"

##```

##0.Get GWAS Catalog colorectal cancer publications
##```{r}
GC <- readRDS(GWASCatalog)


GC_colorectal <- GC[str_detect(GC$"anchor1.MAPPED_TRAIT_URI","EFO_0005842")]


paste0(c("The number of SNPs associated with CRC is:",
         length(unique(GC_colorectal$"anchor1.SNPS"))),
       collapse = "")



paste0(c("The number of SNPs associated with CRC and PMID: 29228715 is:",
         length(GC_colorectal$"anchor1.SNPS"[which(GC_colorectal$"anchor1.PUBMEDID"=="29228715")]
         )), collapse = "")

####################################################

paste0(c("The number of genes associated with CRC is:",
         length(unique(GC_colorectal$"anchor1.REPORTED GENE(S)"))),
       collapse = "")



paste0(c("The number of genes associated with CRC and PMID: 29228715 is:",
         length(GC_colorectal$"anchor1.REPORTED GENE(S)"[which(GC_colorectal$"anchor1.PUBMEDID"=="29228715")]
         )), collapse = "")




tapply(GC_colorectal$"anchor1.INITIAL SAMPLE SIZE",
       GC_colorectal$"anchor1.PUBMEDID",
       unique)

table(GC_colorectal$"anchor1.PUBMEDID")

GC_colorectal

# 29228715 
# 26151821 - 18,299 European ancestry cases, 19,656 European ancestry controls;19 genes

##```


##1.Get GWAS Catalog colorectal cancer genes
##```{r}

# extract GWAS genes for CRC
GWASCatalogGeneSets <- getGWASCatalogGeneSets(terms=c("EFO:0005842"),
                                              GWASCatalog=GWASCatalog,
                                              SNPannotations=GWAS_PATHS)

GWASCatalogGeneSets <- lapply(GWASCatalogGeneSets,
                              function(x){return((unique(x)))})




sapply(GWASCatalogGeneSets,length)
GWASCatalogGeneSets <- GWASCatalogGeneSets[names(GWASCatalogGeneSets)%in%c("EnhancerAtlas","FOCS","JEME","GeneHancer")]

# manual check
GC_colorectal$"anchor1.DISEASE/TRAIT"[which(GC_colorectal$"anchor1.PUBMEDID"=="29228715")] #if everything equals colorectal cancer then OK - OK      

# PLOT: identifying genes that are predicted as CRC with at least at least 4 method


#```{r}

GeneSets <- GWASCatalogGeneSets
# select top10 genes that are in the
TOPGenes=sort(table(unlist(GWASCatalogGeneSets)))
TOPGenes <- names(TOPGenes)[TOPGenes>=3]

GeneSetsDF <- cbind(Method=names(unlist(GeneSets)),
                    Genes=unlist(GeneSets))

GeneSetsDF[,"Method"] <- str_replace_all(GeneSetsDF[,"Method"],"[0-9]*","")

GeneSetsDF <- table(GeneSetsDF[,"Method"],GeneSetsDF[,"Genes"])

# Select genes that appear in 3 or more methods
GeneSetsDF <- GeneSetsDF[,colnames(GeneSetsDF)%in%TOPGenes]

dim(GeneSets)
library(reshape)

melted_cormat <- melt(GeneSetsDF)
colnames(melted_cormat) <- c("Method","Gene","value")

melted_cormat$Method <- as.character(melted_cormat$Method)

library(ggplot2)

path="/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/fig7CRC_pergene_perPublicationStat_red_white_stat.png"

png(path,height =900, width = 700)

ggplot(data = melted_cormat, 
       aes(Method,Gene ,  fill = value))+
  
  geom_tile(color = "white")+
  
  scale_fill_gradient2(low = "white", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,1), space = "Lab", 
                       name="Gene identified Y/N") +
  ylab("")+
  xlab("")+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 40, hjust = 1,face="bold"),
        axis.title = element_text(size=50,face="bold"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y =element_text(size=40,face="bold"),
        legend.text = element_text(size=28,face="bold"),
        legend.title =  element_text(size=30,face="bold"))+
  #geom_text(aes(Gene , Method, label = value), 
  #         color = "black", size =10) +
  coord_fixed()


dev.off()



######################################################
# plot supplementary table of CRC genes

genes.df <- (as.data.frame(sapply(GWASCatalogGeneSets, 
                                  paste,
                                  collapse=", ")))
colnames(genes.df) <- "Genes"
genes.df <- apply(genes.df,2,as.character)

rownames(genes.df) <- names(GWASCatalogGeneSets)

PerGeneTOPcells <- splitTable(genes.df,100)

png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/SupplT_CRC_genes_per_EGA.png",
    height=1200, width=800)
library(gridExtra)
p<-tableGrob(PerGeneTOPcells)
grid.arrange(p)
dev.off()


###############################################
# plot pairwise comparison



#```{r}


pairwiseComp <- function(geneList){
  
  #list=GWASCatalogGeneSets  
  df <- matrix(NA,
               length(geneList),
               length(geneList))
  
  for (i in (1:length(geneList))){
    for (j in (1:length(geneList))){
      
      df[i,j] <- sum(geneList[[i]]%in%geneList[[j]])
      
      colnames(df) <- rownames(df) <- names(geneList)
    }
  }
  return(df)
}

PairwiseComparison <- pairwiseComp(GWASCatalogGeneSets)
#PairwiseComparison_pmid <- pairwiseComp(GWASCatalogGeneSets_1PMID)

t <- round(100*PairwiseComparison/c(PairwiseComparison[1,1],
                     PairwiseComparison[2,2],
                     PairwiseComparison[3,3],
                     PairwiseComparison[4,4]))



pheatmap(PairwiseComparison,number_format = "%.0f", display_numbers = T)
pheatmap(t,number_format = "%.0f", display_numbers = T)
#pheatmap(PairwiseComparison_pmid,number_format = "%.0f", display_numbers = T)



library(reshape2)
melted_cormat <- melt(t)
colnames(melted_cormat) <- c("Method","Method2","value")



# Heatmap
library(ggplot2)

path="/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/figCRC_heatmap_across_methods.png"

png(path,height =800, width = 1100)

ggplot(data = melted_cormat, aes(Method , Method2, fill = value))+
  
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 50, limit = c(0,100), space = "Lab", 
                       name="Gene overlap (%)") +
  ylab("")+
  xlab("")+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 40, hjust = 1,face="bold"),
        axis.title = element_text(size=50,face="bold"),
        axis.text.y =element_text(size=40,face="bold"),
        legend.text = element_text(size=28,face="bold"),
        legend.title =  element_text(size=30,face="bold"))+
  geom_text(aes(Method , Method2, label = value), 
            color = "black", size =10) +
  coord_fixed()


dev.off()

