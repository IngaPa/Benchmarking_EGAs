gwasCatalog <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_rtracklayer_19_02_05_hg19_NearestGeneAdded.rds")




# IDEA: remove all SNPs that are in coding regions + get statistics
#import Annotations
gene.parts <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GENCODE/geneAnnotations_grlist.rds")


  coding.region <- c(reduce(gene.parts$exons),
                     reduce(gene.parts$promoters))


  coding.region.gwas <- data.frame(findOverlaps(anchorOne(gwasCatalog),
                                                coding.region))

  
  #unique(anchorOne(gwasCatalog.reduced))
  # there are 60,132 non-coding SNPs in the GWAS Catalog
  
  gwasCatalog.reduced <- gwasCatalog[-(unique(coding.region.gwas$queryHits))]

  
  # saveRDS(gwasCatalog.reduced,
  #         "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_20_12_05_non_codingSNPs.rds")
  # 
  
  # manually selected top10 diseases with highest number of associated SNPs

  head(sort(table(gwasCatalog.reduced$"anchor1.DISEASE/TRAIT"),decreasing = T),50)
  #tail(sort(table(gwasCatalog$"anchor1.DISEASE/TRAIT"),decreasing = T),50)
  
TOP10diseases <- c("Adolescent idiopathic scoliosis",
                   "Schizophrenia",
                   "Breast cancer",
                   "Type 2 diabetes",
                   "Neuroticism",
                   "Coronary artery disease",
                   "Systemic lupus erythematosus",
                   "Crohn's disease",
                   "Inflammatory bowel disease",
                   "Prostate cancer")



df <- sort(table(gwasCatalog.reduced$"anchor1.DISEASE/TRAIT"[gwasCatalog.reduced$"anchor1.DISEASE/TRAIT"%in%TOP10diseases]),decreasing = T)
names(df) <- c("AIS",
                  "SCZ",
                  "BC",
                  "T2D",
                  "Neuroticism",
                  "CAD",
                  "SLE",
                  "CD",
                  "IBD",
                  "PC")

df <- data.frame(df)



#########################
#---- top10 diseases plot
#########################

path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure6_downstream/top10_diseases.png"
png(path,width = 700, height = 500)



ggplot(df, aes( y=Freq, x=Var1,fill=Var1)) + 
  geom_bar(position="dodge", stat="identity",width = 0.7,color="black")+
  xlab("Diseases") +
  ylab("N[SNPs]\n") +
  theme(legend.position = "none",
        #plot.margin = unit(c(1,1,1.5,1.5), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=30,angle = 45, hjust = 1),
        text = element_text(size=40,face="bold"))+
scale_fill_brewer(palette="Blues")

dev.off()




##################################################
#---- analysis of genes - how much nearest, unknown, etc.---
##################################################


gwas.genes <- str_split(gwasCatalog.reduced$"anchor1.REPORTED GENE(S)","/,")
     


#focusing only on reported genees

gene2nearest <- lapply(1:length(gwasCatalog.reduced),
                       function(i,gwasCatalog.reduced){
  

  tmp <- str_detect(gwasCatalog.reduced$"anchor1.REPORTED GENE(S)"[i],
             gwasCatalog.reduced$nearestGene.name2[i])
  # tmp2 <- str_detect(gwasCatalog.reduced$"anchor1.MAPPED_GENE"[i],
  #            gwasCatalog.reduced$nearestGene.name2[i])
  
  # return(c(tmp|tmp2))
  return(tmp)
  
  },gwasCatalog.reduced=gwasCatalog.reduced)


N.nearest <- sum(unlist(gene2nearest),na.rm = T)
N.other <- sum(!unlist(gene2nearest),na.rm = T)

# 25% genes could be tracked to the nearest gene
N.nearest/(N.nearest+N.other)


unreported <- sum(gwasCatalog.reduced$"anchor1.REPORTED GENE(S)"%in%c("NR",
                                                                      "intergenic",
                                                                      "Intergenic",
                                                                      "genic"))

reported <- sum(!gwasCatalog.reduced$"anchor1.REPORTED GENE(S)"%in%c("NR",
                                                                     "intergenic",
                                                                     "Intergenic",
                                                                     "genic"))


NR <- unreported/(N.nearest+N.other)
nearest <- N.nearest/(N.nearest+N.other)
different <- 1-((N.nearest+unreported)/(N.nearest+N.other))


df.f <- data.frame(c(different,NR,nearest))
colnames(df.f) <- "value"
df.f$Categories <- as.factor(c("Different","NR","Nearest"))



path="/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/GWAScat_genes_cat_NR_nearest_diff.png"
png(path,width = 900, height = 800)


color <- c("white","lightred","red3")



ggplot(df.f, aes(x="", y=value, fill=Categories,color=Categories))+
  geom_bar(width = 0.3, stat = "identity")+
  xlab("\n  GWAS  genes") +
  ylab("% gene overlap \n") +
  theme(#legend.position = "none",
    plot.margin = unit(c(1,1,2.5,2.5), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size=45,face="bold"),
    axis.line = element_line(colour = "black"),
    text = element_text(size=50,face="bold"))+
  scale_fill_manual("Gene Cat.", values = c("Different" = "ivory2", 
                                            "NR" = "tomato", 
                                            "Nearest"="red3"))+
   scale_color_manual("Gene Cat.", values = c("Different" = "ivory2", 
                                             "NR" = "tomato", 
                                             "Nearest" = "red3"))


dev.off()














##################
# analysis for mapped genes
gene2nearest <- lapply(1:length(gwasCatalog.reduced),
                       function(i,gwasCatalog.reduced){
                         
                         
                         # tmp <- str_detect(gwasCatalog.reduced$"anchor1.REPORTED GENE(S)"[i],
                         #                   gwasCatalog.reduced$nearestGene.name2[i])
                         tmp2 <- str_detect(gwasCatalog.reduced$"anchor1.MAPPED_GENE"[i],
                                     gwasCatalog.reduced$nearestGene.name2[i])
                         
                         # return(c(tmp|tmp2))
                         return(tmp2)
                         
                       },gwasCatalog.reduced=gwasCatalog.reduced)


N.nearest <- sum(unlist(gene2nearest),na.rm = T)
N.other <- sum(!unlist(gene2nearest),na.rm = T)

# 36% genes could be tracked to the nearest gene
N.nearest/(N.nearest+N.other)



1-(unreported/(unreported+reported))