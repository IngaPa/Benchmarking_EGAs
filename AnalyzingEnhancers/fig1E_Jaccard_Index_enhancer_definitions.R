#calculating Jaccard index

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


# AND count jaccard

Publications <- lapply(path,readRDS)
  names(Publications) <-  str_replace(str_replace(path,".rds",""),"./","")

  
  
JI <- matrix(NA, nrow=4,ncol=4, dimnames = list(names(Publications) ,names(Publications) ))


JI <-lapply(names(Publications), function(x,Publications){
  
  lapply(names(Publications), function(y,Publications){
  

    x1 <- Publications[[x]]
    unique1 <- unique(anchorOne(x1))
    x2 <- Publications[[y]]
    unique2 <- unique(anchorOne(x2))
    
    Intersection <- intersect(unique1,unique2)
    Union <- union(unique1,unique2)
    
    
    return(length(Intersection)/length(Union))
  
},Publications=Publications)
  
},Publications=Publications)

JI.df <- do.call("rbind.data.frame",JI)
colnames(JI.df) <- rownames(JI.df) <- names(Publications)
JI.df <- round(JI.df*100,1)



# get triangular plot

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(JI.df)
upper_tri


# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri)
melted_cormat$Publication <- colnames(upper_tri)
melted_cormat <- melted_cormat[complete.cases(melted_cormat),]


# Heatmap
library(ggplot2)

path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/fig1E_JaccardIndex.png"

png(path,height =680, width = 950)

ggplot(data = melted_cormat, aes(variable , Publication, fill = value))+
  
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 50, limit = c(0,100), space = "Lab", 
                       name="Jaccard index (%)") +
  ylab("")+
  xlab("")+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 40, hjust = 1,face="bold"),
        axis.title = element_text(size=50,face="bold"),
        axis.text.y =element_text(size=40,face="bold"),
        legend.text = element_text(size=28,face="bold"),
        legend.title =  element_text(size=30,face="bold"))+
  geom_text(aes(variable, Publication, label = value), color = "black", size =13) +
  coord_fixed()


dev.off()

