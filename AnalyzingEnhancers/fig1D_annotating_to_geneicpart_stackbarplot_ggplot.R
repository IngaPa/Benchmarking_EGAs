

# IDEA:
annotate enhancer regions from different methods together

.libPaths("/home/ipatarc/R/x86_64-redhat-linux-gnu-library/3.5")
library(GenomicInteractions)
library(rtracklayer)
library(genomation)
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/projectHelpFunctions.R")




#import Annotations
gene.parts <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GENCODE/geneAnnotations_grlist.rds")

setwd("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/")

path <- list.files(full.names = T,pattern="rds")







perMethodStat <- lapply(path, function(x,
                                       geneAnnotations){
  print(x)
  # x="./EnhancerAtlas.rds"
  regions=readRDS(x)
  
  # process (remove chrY)
  peaks=anchorOne(regions)
  peaks <- peaks[!seqnames(peaks)=="chrY"]
  peaks <- peaks[!seqnames(peaks)=="chrM"]
  
   return(percentageannotateWithGeneAnnotations(regions=peaks,
                                                geneAnnotations))

},geneAnnotations=gene.parts)


names(perMethodStat) <- str_replace(str_replace(path,"./",""),".rds","")



genePartPercentage <- do.call("cbind.data.frame",perMethodStat)





library(pheatmap)            

pheatmap(genePartPercentage,cluster_rows = F,cluster_cols = F,display_numbers = T, number_format = "%.2f")

library(ggplot2)
library(reshape)
library(stringr)

TMP <- melt(genePartPercentage)
TMP$geneAnnot <- rownames(genePartPercentage)
TMP$value <- round(TMP$value,1)
TMP$Publication <- TMP$variable 

TMP.SS <- TMP[TMP$variable%in%c("GeneHancer","EnhancerAtlas","JEME","FOCS"),]
TMP.SS$Annotation <- str_replace(TMP.SS$geneAnnot,"Annotated","")

    colorScheme <- c("black","darkorange","blue","darkcyan" )
    names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")
TMP.SS$color <- colorScheme[as.character(TMP.SS$variable)]
  
#per Experiment statistics

TMP.SS <- transform(TMP.SS, 
                    Annotation  = factor(
                      Annotation ,
                      levels=c( 'introns','intergenic','exons', 'promoters'),
                      ordered =TRUE))

#############################
# version 1 plot
png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/fig1D_genecontext_barplot_ggplot.png",
    width = 950,height = 600)

ggplot(TMP.SS,aes(x=Publication,y=value,
                  fill=Annotation,order=Publication), color=color) +  
  stat_summary(fun.y=mean,position="dodge",geom="bar") +  
  ylab("%[enhancers in genetic context]\n") +
  xlab("\nPublication") + 
  scale_fill_brewer(palette = "Oranges") + 
  theme_bw()+
  theme(text = element_text(size=30,face="bold"),
        axis.text =element_text(size=30,face="bold"),
        legend.text = element_text(size=20,face="bold"))
  #  geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()



#############################
# version 2 plot
png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/fig1D_genecontext_barplot_ggplot_v2.png",
    width = 1000,height = 500)

colors <- colorScheme[levels(TMP.SS$Publication)]
colors <- colors[complete.cases(colors)]


ggplot(TMP.SS,aes(x=Annotation,y=value,
                  fill=Publication,order=Publication)) +  
  stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
  scale_fill_manual(values =  alpha(colors, .7))+
 
  ylab("% annotated enhancers \n") +
  xlab("\nGenetic context") + 
  
  theme_bw()+
  theme(text = element_text(size=30,face="bold"),
        axis.text =element_text(size=30,face="bold"),
        legend.text = element_text(size=30,face="bold"))
#  geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()


#############################
# version 2 plot - legend

png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/fig1D_genecontext_barplot_ggplot_nolegend.png",
    width = 950,height = 500)

ggplot(TMP.SS,aes(x=Annotation,y=value,
                  fill=Publication,order=Publication)) +  
  stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
  scale_fill_manual(values =  alpha(colors, .7))+
  
  ylab("% annotated enhancers \n") +
  xlab("\n Genetic context") + 
  
  theme_bw()+
  theme(text = element_text(size=30,face="bold"),
        axis.text =element_text(size=30,face="bold"),
        legend.position="none")
#  geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()




