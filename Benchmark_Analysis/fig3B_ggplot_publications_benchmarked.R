# script that plots ratio between counted benchmarked interactions vs full N of interactions
# N(Benchmarked EGA)/N(EGA) for eacj benchmark dataset
# xaxis is benchmark dataset 



###########################
# -------------------------
# Input
###########################

paths=c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/EnhancerAtlasBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/FOCSBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/GeneHancerBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/JEMEBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/inhousemBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/flexiblemBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/stringentmBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/CCSIBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/PCHiCBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/GTExBenchVal2BenchVal.rds",
        "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/WestraeQTLsBenchVal2BenchVal.rds")



names(paths) <- c("EnhancerAtlas","FOCS","GeneHancer","JEME","InhouseM","FlexibleC","StringentC","CCSI","PC.HiC","GTEx","Westra.eQTLs")





###########################
# -------------------------
# libraries and scripts
###########################
library(ggplot2)
# library(reg2gene)
library(stringr)
library(parallel)
library(GenomicRanges)
library(GenomicInteractions)
library(InteractionSet)

getRankingOrginalPublication <- function(path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking/PCHiCBenchVal2BenchVal.rds",
                                         names=c("JEME","GeneHancer","EnhancerAtlas","FOCS",
                                                 "in.house.m.","flexible.m.","stringent.m.",
                                                 "GTEx","Westra.eQTLs","CCSI","PC.HiC","Gasperini.TP")){
  
  # get benchmarking statistics
  HIC.benchmarking <- readRDS(path)
  
  metadata <- as.data.frame(mcols(HIC.benchmarking))[,names]
  
  statistics <- colSums(metadata)
  
  names(statistics) <- str_replace(str_replace(str_replace(names(statistics), "stringent.m.","StringentC"),
                                               "flexible.m.","FlexibleC"),
                                   "in.house.m.","InhouseM")
  
  return(statistics)
  
}





###########################
# -------------------------
# analysis
###########################

bench_results <- lapply(paths,getRankingOrginalPublication,names=c("JEME","GeneHancer","EnhancerAtlas","FOCS",
                                                                   "in.house.m.","flexible.m.","stringent.m.",
                                                                   "GTEx","Westra.eQTLs","CCSI","PC.HiC","Gasperini.TP"))
b2b <- do.call("cbind.data.frame",bench_results)




filter_results <- lapply(paths,
                         getRankingOrginalPublication,
                         names=c("Filter_JEME","Filter_GeneHancer","Filter_EnhancerAtlas","Filter_FOCS","Filter_in.house.m.","Filter_flexible.m.","Filter_stringent.m." , 
                                 "Filter_GTEx" , "Filter_Westra.eQTLs" , "Filter_CCSI","Filter_PC.HiC","Filter_Gasperini.TP"))

filter_res_df <- do.call("cbind.data.frame",filter_results)

colnames(filter_res_df)==colnames(b2b)
rownames(filter_res_df)==rownames(b2b)




df <-t( 100*b2b/filter_res_df)
colnames(df) <-str_replace(colnames(df),"Westra.eQTLs","W.eQTLs")
colnames(df) <-str_replace(colnames(df),"PC.HiC","PCHiC")


df<- df[,colnames(df)%in%c("GTEx","W.eQTLs","PCHiC","CCSI")]
df <- df[!(rownames(df)%in%c("GTEx","Westra.eQTLs","PC.HiC","CCSI")),]
library(reshape)
to.plot <- melt(df)
colnames(to.plot) <- c("Publication","Benchmark","value")



#############################
# version 1 plot

colorScheme <- c("black","darkorange","blue","darkcyan","grey","red","darkred")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas",
                        "InhouseM","FlexibleC","StringentC")


png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure4_benchmarking/fig4C_publicationsBench_barplot_ggplot.png",
    width = 950,height = 600)



ggplot(to.plot,aes(x=Benchmark,y=value,
                   fill=Publication,order=Benchmark), color=colorScheme[levels(to.plot$Publication)])  +  
  stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
  scale_fill_manual(values =  alpha(colorScheme[levels(to.plot$Publication)], .9))+
  ylab("%[Benchmarked EGA]\n") +
  xlab("\n Benchmark dataset") + 
  theme(text = element_text(size=30,face="bold"),
        axis.text =element_text(size=30,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=25,face="bold"))
#  geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()






#############################
# benvh2benvh heatmap
df <-t( 100*b2b/filter_res_df)


df <- df[,colnames(df)%in%c( "GTEx","Westra.eQTLs" , "CCSI","PC.HiC" )]
df <- df[rownames(df)%in%c( "GTEx","Westra.eQTLs" , "CCSI","PC.HiC" ),]


df <- round(df,1)


to.plot <- melt(df)
colnames(to.plot) <- c("Benchmark1","Benchmark2","value")


path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure4_benchmarking/bench2benc.png"

png(path,height =800, width = 1000)

ggplot(data = to.plot, aes(Benchmark1 , Benchmark2, fill = value))+
  
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "orange", 
                       midpoint = 50, limit = c(0,100), space = "Lab", 
                       name="Coverage of \n benchmark data(%)") +
  ylab("")+
  xlab("")+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 40, hjust = 1,face="bold"),
        axis.title = element_text(size=50,face="bold"),
        axis.text.y =element_text(size=40,face="bold"),
        legend.text = element_text(size=28,face="bold"),
        legend.title =  element_text(size=30,face="bold"))+
  geom_text(aes(Benchmark1, Benchmark2, label = value), 
            color = "black", size =10) +
  coord_fixed()

dev.off()


#################################
###############################
# benchnmark coverage by EGAs


df <-t( 100*b2b/filter_res_df)
colnames(df) <-str_replace(colnames(df),"Westra.eQTLs","W.eQTLs")
colnames(df) <-str_replace(colnames(df),"PC.HiC","PCHiC")


df<- df[,!colnames(df)%in%c("GTEx","W.eQTLs","PCHiC","CCSI","Gasperini.TP")]
df <- df[(rownames(df)%in%c("GTEx","Westra.eQTLs","PC.HiC","CCSI")),]
library(reshape)
to.plot <- melt(df)
colnames(to.plot) <- c("Benchmark","Publication","value")



#############################
# version 2 plot

colorScheme <- c("black","darkorange","blue","darkcyan","grey","red","darkred")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas",
                        "InhouseM","FlexibleC","StringentC")


png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure4_benchmarking/fig4D_barplot_benchBYegaS_ggplot.png",
    width = 1000,height = 600)



ggplot(to.plot,aes(x=Benchmark,y=value,
                   fill=Publication,order=Benchmark), color=colorScheme[levels(to.plot$Publication)])  +  
  stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
  scale_fill_manual(values =  alpha(colorScheme[levels(to.plot$Publication)], .9))+
  ylab("%[Benchmarked eQTLs or CI]\n") +
  xlab("\n Benchmarked dataset") + 
  theme(text = element_text(size=30,face="bold"),
        axis.text =element_text(size=30,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=25,face="bold"))
#  geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()

