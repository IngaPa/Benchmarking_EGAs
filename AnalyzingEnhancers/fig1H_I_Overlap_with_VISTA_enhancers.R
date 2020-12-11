# idea, identify overlap in enhancer regions VISTA enhancers
.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.5/")
library(GenomicInteractions)
library(stringr)
library(ggplot2)
library(reshape)
library(rtracklayer)
library(reshape)


Vista <- read.table("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Data/Vista_Enhancers_201026/Human_1938_Vista_enhancers.txt",sep="\t")
Vista.gr <- GRanges(as.character(Vista$V2))

# properties of VISTA enhancers

(sum(coverage(unique(Vista.gr))))
# 3705240 bp is covered by Vista enancers
summary(width(Vista.gr)) # max 8,062 bp, min 315 enhncers


##########################################3
# publications statistics

Datasets <- c("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/EnhancerAtlas.rds",
              "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/FOCS.rds",
              "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/GeneHancer.rds",
              "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/JEME.rds")

EGAs <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/EnhancerAtlas.rds")


N_enhancer_overlap <- sapply(Datasets, function(x){
  EGAs <- readRDS(x)
  N_confirmed_EGAs_enh <- as.data.frame(findOverlaps(anchorOne(EGAs),Vista.gr))
  return(length(unique(N_confirmed_EGAs_enh$queryHits)))
})
N_enhancer_overlap



VISTA_enhancer_COVERAGE_binart <- sapply(Datasets, function(x){
  EGAs <- readRDS(x)
  N_confirmed_EGAs_enh <- as.data.frame(findOverlaps(anchorOne(EGAs),Vista.gr))
  return(length(unique(N_confirmed_EGAs_enh$subjectHits)))
})
# 1083 EA, 237 FOCS, 1069 Genehancer, 365 JEME, 497 EN


##########################################################
# cverage of gasperini et al. enhancers


GASPERINI_enh <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Data/Gasperini/Gasperini_S2A_tested_enhancers_GRanges.rds")
Gasp_total_coverage <- sum(sum(coverage(unique(GASPERINI_enh))))
# 2,766,374 enhancers


##############################################
# n enhancers overlapping Gasperini enhancers
N_enhancer_overlap_GASP <- sapply(Datasets, function(x){
  EGAs <- readRDS(x)
  N_confirmed_EGAs_enh <- as.data.frame(findOverlaps(anchorOne(EGAs),GASPERINI_enh))
  return(length(unique(N_confirmed_EGAs_enh$queryHits)))
})
N_enhancer_overlap_GASP



# previous information: hard-coded 
tmp <- cbind(VISTA_enhancer_COVERAGE_binart,
             GASP_enhancer_COVERAGE_binart)


rownames(tmp) <- str_replace(str_replace(row.names(tmp),"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/",""),".rds","")

tmp <- cbind(tmp)
colnames(tmp) <- c("VISTA","crisprQTL")


table.plot <- melt(tmp)
colnames(table.plot) <- c("Publication","Type","value")

colorScheme <- c("black","darkorange","blue","darkcyan")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")


#############################
# version 1 plot
 png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/Publications/Enhancers_VISTA_Gasperini_enhOverlaps.png",
     width = 800,height = 650)


ggplot(table.plot, aes(fill=Publication, y=value, x=Type)) + 
  geom_bar(position="dodge", stat="identity",size=0.7)+
  ylab("Confirmed enhancers\n") +
  xlab("\n Source of enhancers") +
  theme(#legend.position = "none",
    plot.margin = unit(c(1,1,1,2.5), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90),
    axis.line = element_line(colour = "black"),
    text = element_text(size=40,face="bold"))+ 
  scale_fill_manual("Publication", values = c("FOCS" = "black", "GeneHancer" = "darkorange", "JEME" = "blue","EnhancerAtlas"="darkcyan"))+
  scale_y_continuous(limits = c(0, 5000),labels = scales::comma)




dev.off()
####################################

###################################################3




##############################################3
# what is per base coverage

VISTA_enhancer_COVERAGE <- sapply(Datasets, function(x){
  EGAs <- readRDS(x)
  
  return(sum(sum(coverage(c(reduce(unique(anchorOne(EGAs))),
                            reduce(unique(Vista.gr))))>1)))
})
names(VISTA_enhancer_COVERAGE) <- str_replace(str_replace(names(VISTA_enhancer_COVERAGE),"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/",""),".rds","")
VISTA_enhancer_COVERAGE

# general coverage per publication
EGA_enhancer_COVERAGE <- sapply(Datasets, function(x){
  EGAs <- readRDS(x)
  
  sum(sum(coverage(reduce(anchorOne(EGAs)))))
                            
})
names(EGA_enhancer_COVERAGE) <- str_replace(str_replace(names(VISTA_enhancer_COVERAGE),"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/",""),".rds","")
EGA_enhancer_COVERAGE



(100*VISTA_enhancer_COVERAGE/EGA_enhancer_COVERAGE)

Vista_enhancers_coveraged_perc <- (100*VISTA_enhancer_COVERAGE/sum(sum(coverage(unique(Vista.gr)))))




###################################

# n Gasp enhancers overlapping EGA enhancers
GASP_enhancer_COVERAGE_binart <- sapply(Datasets, function(x){
  EGAs <- readRDS(x)
  N_confirmed_EGAs_enh <- as.data.frame(findOverlaps(anchorOne(EGAs),GASPERINI_enh))
  return(length(unique(N_confirmed_EGAs_enh$subjectHits)))
})
names(GASP_enhancer_COVERAGE_binart) <- str_replace(str_replace(names(GASP_enhancer_COVERAGE_binart),"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/",""),".rds","")


# what is per base coverage

GASP_enhancer_COVERAGE <- sapply(Datasets, function(x){
  EGAs <- readRDS(x)
  
  return(sum(sum(coverage(c(reduce(unique(anchorOne(EGAs))),
                            reduce(unique(GASPERINI_enh))))>1)))
})
names(GASP_enhancer_COVERAGE) <- str_replace(str_replace(names(GASP_enhancer_COVERAGE),"/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/0_Pooled_input_Datasets/",""),".rds","")

Gasp_enhancers_coveraged_perc <- 100*GASP_enhancer_COVERAGE/Gasp_total_coverage


tmp <- cbind(Vista_enhancers_coveraged_perc,Gasp_enhancers_coveraged_perc)
colnames(tmp) <- c("VISTA","crisprQTL")


table.plot <- melt(tmp)
colnames(table.plot) <- c("Publication","Type","value")

colorScheme <- c("black","darkorange","blue","darkcyan")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")


#############################
# version 1 plot
png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/Publications/Enhancers_VISTA_Gasperini_enhcoverage.png",
    width = 800,height = 650)


ggplot(table.plot, aes(fill=Publication, y=value, x=Type)) + 
  geom_bar(position="dodge", stat="identity",size=0.7)+
  ylab("Enhancer coverage %\n") +
  xlab("\n Source of enhancers") +
  theme(#legend.position = "none",
    plot.margin = unit(c(1,1,1,2.5), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90),
    axis.line = element_line(colour = "black"),
    text = element_text(size=40,face="bold"))+ 
  scale_fill_manual("Publication", values = c("FOCS" = "black", "GeneHancer" = "darkorange", "JEME" = "blue","EnhancerAtlas"="darkcyan"))+
  scale_y_continuous(limits = c(0, 100),labels = scales::comma)




dev.off()















#######################
# overlap of our models with VISTA enhnacers
All_Tested_Models <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Data/Models_full_EGAs_9mil.rds")

elasticnet <-  All_Tested_Models[All_Tested_Models$votedelasticnet!=0]
failed_elasticnet <-  All_Tested_Models[All_Tested_Models$votedelasticnet!=0]

fullEnhCoverageEN <- sum(sum(coverage(reduce(anchorOne(All_Tested_Models)))))

EN_confirmed_EGAs_enh <- as.data.frame(findOverlaps(anchorOne(All_Tested_Models),Vista.gr))
(length(unique(EN_confirmed_EGAs_enh$subjectHits))) #988 Vista enhancers is confirmed by 
(length(unique(EN_confirmed_EGAs_enh$queryHits))) #  37,321 enhancers from EN is confirmed by VISTA enhancers
Enh.coverage <- sum(sum(coverage(c(reduce(unique(anchorOne(All_Tested_Models))),
                                   reduce(unique(Vista.gr))))>1))
#970,137 bp is covered by EN enhancers and Vista.GR
100*Enh.coverage/fullEnhCoverageEN
100*Enh.coverage/sum(sum(coverage(unique(Vista.gr))))


# elastic net models
# 
EN.Enh.coverage <- sum(sum(coverage(c(reduce(unique(anchorOne(elasticnet))),
reduce(unique(Vista.gr))))>1))





###################################################
EN_confirmed_EGAs_Gaspenh <- as.data.frame(findOverlaps(anchorOne(All_Tested_Models),GASPERINI_enh))
(length(unique(EN_confirmed_EGAs_Gaspenh$subjectHits))) #3850 GASP enhancers is confirmed by 
(length(unique(EN_confirmed_EGAs_Gaspenh$queryHits))) #  37,321 enhancers from EN is confirmed by VISTA enhancers
GaspEnh.coverage <- sum(sum(coverage(c(reduce(unique(anchorOne(All_Tested_Models))),
                                   reduce(unique(GASPERINI_enh))))>1))
#970,137 bp is covered by EN enhancers and Vista.GR
100*Enh.coverage/fullEnhCoverageEN
100*GaspEnh.coverage/sum(sum(coverage(unique(Vista.gr))))


# elastic net models
# 
GaspEN.Enh.coverage <- sum(sum(coverage(c(reduce(unique(anchorOne(elasticnet))),
                                      reduce(unique(GASPERINI_enh))))>1))


EN_confirmed_EGAs_Gaspenh <- as.data.frame(findOverlaps(anchorOne(elasticnet),GASPERINI_enh))
(length(unique(EN_confirmed_EGAs_Gaspenh$subjectHits))) #3850 GASP enhancers is confirmed by 
(length(unique(EN_confirmed_EGAs_Gaspenh$queryHits))) 


GaspEN.Enh.coverage/sum(sum(coverage(unique(Vista.gr))))


# 