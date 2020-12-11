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

gwascatalog <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_rtracklayer_19_02_05_hg19_NearestGeneAdded.rds")
# annotated GWAS Catalog; annotated by different EGA models


annotatedGWASC <- list.files("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/GWAS_Catalog/annotatedGWAS/",
                             full.names = T,"rds")[-c(2,5,7)]


files <- lapply(annotatedGWASC,readRDS)

names(files) <- str_replace(str_replace(basename(annotatedGWASC),"GWASC_annotatedBy_",""),".rds","")
# repository of EGA interactions


annotations <- do.call("cbind.data.frame",
                       lapply(files,function(x){
  #x=files[[1]]
   return(((table(x$annotatedAs))))

  }))

annotations <- annotations[,str_detect(colnames(annotations),"Freq")]
rownames(annotations) <- c("enhancer","nearestGene","promoter")
colnames(annotations) <- str_replace(colnames(annotations),".Freq","")




# df

# png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SuppementaryTables/SupplementaryTable_FullGWAS_Annotation_Genetic_Context.png",
#     height=110, width=590)
# library(gridExtra)
# p<-tableGrob(annotations)
# grid.arrange(p)
# dev.off()



tmp <- t(round(t(annotations)/colSums(annotations)[names(annotations)],3)*100)

# png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/SuppementaryTables/SupplementaryTable_FullGWAS_Annotation_Genetic_Context_percentage.png",
#     height=110, width=590)
# library(gridExtra)
# p<-tableGrob(tmp)
# grid.arrange(p)
# dev.off()


# plot

library(reshape)
annotations2 <- melt(annotations)
annotations2$annotation <- rownames(annotations)
colnames(annotations2) <- c("Publication","value","Annotation")


annotations2$Percentage <- round(annotations2$value/
          rep(tapply(annotations2$value,annotations2$Publication,sum),each=3),2)


annotations2$Publication <- as.character(annotations2$Publication)

# change names for this
annotations2$v2 <- round(annotations2$value/1000,1)

annotations2$Annotation[annotations2$Annotation=="enhancer"] <- "2. Enhancer"
annotations2$Annotation[annotations2$Annotation=="promoter"] <- "1. Promoter"
annotations2$Annotation[annotations2$Annotation=="nearestGene"] <- "3. NearestG"

#per Experiment statistics

 png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/GWAS_Annotations_Publ_v2.png",
     width = 870, height = 700)

        par(font.axis = 2)
        par(font.lab=2)
        par(font=2)
        
        ggplot(annotations2,aes(x=Publication,
                                y=Percentage,
                                fill=Annotation), color=Annotation) +  
          stat_summary(fun.y=mean,
                       position="stack",
                       geom="bar",
                       width=0.55) +  
            # coord_flip() +  
          ylab("% of annotation") +
          xlab("Publication") + 
          scale_fill_brewer(palette = "RdBu") + 
          theme(text = element_text(size=50,colour = "black",face="bold"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.title.x=element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text.x = element_text(size=40,angle = 45, hjust = 1),
                legend.text = element_text(size=35,face="bold"),
                axis.text.y =element_text(colour = "black",face="bold"))#+ 
          # geom_text(aes(label = v2), size = 13,
          #           angle = 90, 
          #           hjust = 0.55, 
          #           vjust = 0.1, position = "stack")


 dev.off()
 
 
 
 