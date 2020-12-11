# get density plots for: distance to nearest TSS and width statistics


library(pheatmap)            
library(GenomicInteractions)
library(rtracklayer)
library(genomation)
library(ggplot2)
library(reshape)
library(stringr)
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/generalPlot_functions.R")


setwd("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/")


list.dirs("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/")

path <- list.files(pattern="rds",full.names = T)[c(1,3,4,6)]

#############################################
# Size distribution PLOT
############################################


enh_regions_distanceEG <- lapply(path, function(x){
  
  EGA <- readRDS(x)
  return(distance(anchorOne(EGA),anchorTwo(EGA)))})


names(enh_regions_distanceEG) <- str_replace(str_replace(path,".rds",""),"./","")

# plotDensity4publications(path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/SizeDistributionEnhancers_Publications.png",
#                          listTOplot <- lapply(enh_regions_width[c("FOCS","GeneHancer","JEME","EnhancerAtlas")],log10),
#                          xlab=" log10(enhancer length) ",
#                          main="",
#                          #main= " Size distribution",
#                          legend.position="topright")



# get lengths
#listTOplot <- lapply(enh_regions_width[c("FOCS","GeneHancer","JEME","EnhancerAtlas")],log10)
plotData <- as.data.frame(unlist(enh_regions_distanceEG))
# get appropriate names
tmp <- sapply(enh_regions_distanceEG, length)
names_Publ <- rep(names(tmp),tmp)

# get appropriate order
ggplotData <- cbind(plotData,names_Publ)
names(ggplotData) <- c("value","Publication")


# colors
colorScheme <- c("black","darkorange","blue","darkcyan" )
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")

# v1
 png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/DistanceDistributionEnhancers_TSS_Publications.png"
     ,width = 900,height = 900)

ggplot(ggplotData, aes(value, 
                       colour = Publication,
                       fill=Publication)) +
  xlim(1,300000)+
  xlab("\n E-TSS distance")+
  ylab("Density \n")+
  geom_density(adjust=8, size=4)+
  scale_color_manual(values =  alpha(colorScheme[levels(ggplotData$Publication)], .7))+
  scale_fill_manual(values =  alpha(colorScheme[levels(ggplotData$Publication)], .1))+
  #theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=50,face="bold"),
        axis.text =element_text(size=40,face="bold"),
        text = element_text(size=40,face="bold", vjust = 2),
        legend.text = element_text(size=40,face="bold"))
  #ylim(0,1.3)



dev.off()


# v1 - legend
png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/DistanceDistributionEnhancers_TSS_Publications_nolegend.png"
    ,width = 900,height = 900)

ggplot(ggplotData, aes(value, 
                       colour = Publication,
                       fill=Publication)) +
  xlim(1,300000)+
  xlab("\n E-TSS distance")+
  ylab("Density \n")+
  geom_density(adjust=8, size=4)+
  scale_color_manual(values =  alpha(colorScheme[levels(ggplotData$Publication)], .7))+
  scale_fill_manual(values =  alpha(colorScheme[levels(ggplotData$Publication)], .1))+
  #theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=50,face="bold"),
        axis.text =element_text(size=40,face="bold"),
        text = element_text(size=40,face="bold", vjust = 2),
        legend.text = element_text(size=20,face="bold"),
        legend.position="none")



dev.off()





# plot data frame
df <- do.call("rbind.data.frame",
              lapply(enh_regions_width[c("FOCS","GeneHancer","JEME","EnhancerAtlas")],summary))

rownames(df) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")
colnames(df) <- c("Min."," 1st Qu.","Median","Mean","3rd Qu.","Max. ")



# plot table
png("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/Suppl_table_Fig1b_densityPlot.png" ,
    height=220, width=500)
library(gridExtra)
p<-tableGrob(df)
grid.arrange(p)
dev.off()











