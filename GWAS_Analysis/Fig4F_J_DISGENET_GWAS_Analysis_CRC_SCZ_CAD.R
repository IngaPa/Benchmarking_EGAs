
.libPaths("~/R/x86_64-redhat-linux-gnu-library/3.5/")
library(stringr)
library(pheatmap)
library(gwascat)
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/generalPlot_functions.R")
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/projectHelpFunctions.R")
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/GWASCatalog/19_06_12_helpFunction_for_GWAS.R")



GWASCatalog="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/GWASCatalog/GWAS_Catalog_rtracklayer_19_02_05_hg19_NearestGeneAdded.rds"
GWAS_PATHS <- "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/GWAS_Catalog/annotatedGWAS_EGA/"
DISGENET <- readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Data/DISGENET/DGN_all_gene_disease_associations_noGWASC.rds") 


###################
#------CRC--------#
###################

##############################
# search GWAS

    GWASCatalogGeneSets_CRC <- getGWASCatalogGeneSets(terms=c("EFO:0005842"),
                                                  GWASCatalog=GWASCatalog,
                                                  SNPannotations=GWAS_PATHS)
    
    GWASCatalogGeneSets_CRC <- lapply(GWASCatalogGeneSets_CRC,
                                  function(x){return((unique(x)))})
    
    
    GWASCatalogGeneSets_CRC <- GWASCatalogGeneSets_CRC[names(GWASCatalogGeneSets_CRC)%in%
                                                         c("FOCS","JEME","EnhancerAtlas","GeneHancer")]

#############################
# how many genes in DISGENE for 1) CRC and 2) CRC ancestor terms
    DISGENETGeneSets_CRC <- getDISGENETGeneSets(efo.terms =c("EFO:0005842"),
                                            umls.terms = c("C1527249"),
                                            DISGENET=DISGENET)
    
    
    DISGENETGeneSets_CRC <- lapply(DISGENETGeneSets_CRC,unique)
    
    
    
    
    # 1676 DISGENET GENES and 2096 DISGENET ANCESTOR GENES
    
    
    stat_CRC <- do.call("rbind.data.frame",lapply(GWASCatalogGeneSets_CRC,function(set1){
      
      #set1 <- GWASCatalogGeneSets[1]
      lapply(DISGENETGeneSets_CRC,function(x){
        
        return(sum(set1%in%x))
      })
    }))
    
    stat_CRC <- rbind(stat_CRC,sapply(DISGENETGeneSets_CRC,length))
    
    stat_CRC$fullSize <- c(sapply(GWASCatalogGeneSets_CRC,length),NA)


##############################
# add overlap with reported genes

CRC_genes <- sort(unique(c("APC","MUTYH","FAP","AFAP","MSH2","MLH1",
                           "POLE","POLD1","STK11","SMAD4","PTEN","BMPR1A",
                           "CHEK2",'VTI1A',"EIF3H0","CCND2","MYOB1","DUSP10",
                           "CDKN1A","NKX2","GREM1","CTNNM1","TCF7L2","NOS1",
                           "CCND2","MYC","BMP2","BMP5","PREX1","C11orf93","DIP2B",
                           "SLC22A2","HOA1","NXN","SMAD7","SHROOM2","CDH1","SH283",
                           "LAMC1","LAMAS","TGFB1","POLD2","PITX1","AS1","GATA3",
                           "LRIG1","MYRF","TPD52L3","ATF1","CD9","MYNN","RHPN2",
                           "CTNNB1","BMP4","MYO1B","SH2B","PITX","POLD3","FEN1",
                           "MLH1","EIF3H","CD9","PLCB1","TERC","DNMT3B","ITIH")))

inCRC_genes <- sapply(GWASCatalogGeneSets_CRC,function(x){
  
  sum(unique(x)%in%CRC_genes)
  
})

totalGenes <- sapply(GWASCatalogGeneSets_CRC,length)
stat_CRC$paper <- c(inCRC_genes,length(CRC_genes))  




#############################
# genes that are confirmed
con.CRC_genes <- sapply(GWASCatalogGeneSets_CRC,function(x){
  
  unique(x)[(unique(x)%in%CRC_genes)]
  
})

table(sort(table(unlist(con.CRC_genes))))
# 75% genes is unique 52/70

# how much genes is confirmed by other methods

df <- do.call("rbind.data.frame",lapply(con.CRC_genes,function(X){
  #
  lapply(con.CRC_genes,function(Y) return(sum(Y%in%X)))
  
}))
pheatmap(df/apply(df,1,max),display_numbers = T)


#############################
# pheatmap of confirmed genes -CRC
#############################
df <- do.call("rbind.data.frame",lapply(con.CRC_genes,function(X){
  #
  lapply(con.CRC_genes,function(Y) return(sum(Y%in%X)))
  
}))
pheatmap(df/apply(df,1,max),display_numbers = T)






df <- do.call("rbind.data.frame",lapply(GWASCatalogGeneSets_CRC,function(X){
  #
  lapply(GWASCatalogGeneSets_CRC,function(Y) return(sum(Y%in%X)))
  
}))
pheatmap(df/apply(df,1,max),display_numbers = T)



#--------------------------------------------------------

###################
#------SCZ--------#
###################

##############################
# search GWAS

    GWASCatalogGeneSets_SCZ <- getGWASCatalogGeneSets(terms=c("EFO:0000692"),
                                                  GWASCatalog=GWASCatalog,
                                                  SNPannotations=GWAS_PATHS)
    
    GWASCatalogGeneSets_SCZ <- lapply(GWASCatalogGeneSets_SCZ,
                                  function(x){return((unique(x)))})
    
    
    GWASCatalogGeneSets_SCZ <- GWASCatalogGeneSets_SCZ[names(GWASCatalogGeneSets_SCZ)%in%c("FOCS","JEME","EnhancerAtlas","GeneHancer")]


#############################
# how many genes in DISGENE for 1) CRC and 2) CRC ancestor terms
    DISGENETGeneSets_SCZ <- getDISGENETGeneSets(efo.terms =c("EFO:0000692"),
                                            umls.terms = "C0036341",
                                            DISGENET=DISGENET)
    
    DISGENETGeneSets_SCZ <- lapply(DISGENETGeneSets_SCZ,unique)
    
    sapply(DISGENETGeneSets_SCZ,length)
    


#############################
# STATISTICS
    DISGNE.SCZ_genes <- sapply(GWASCatalogGeneSets_SCZ,function(x){
      
      unique(x)[(unique(x)%in%DISGENETGeneSets$DISGENET_Genes)]
      
    })
    
    (table(sort(table(unlist(DISGNE.SCZ_genes)))))
    # 69% genes is unique 61/89
    
    
    do.call("rbind.data.frame",lapply(DISGNE.SCZ_genes,function(X){
     #
      lapply(DISGNE.SCZ_genes,function(Y) return(sum(Y%in%X)))
      
    }))
    
    # what this says is taht FOCS is covered by GeneHancer
    # what about JEME, can FE Genehancer and EnhancerAtlas cover JEME????
    
    
    stat_SCZ <- do.call("rbind.data.frame",lapply(GWASCatalogGeneSets_SCZ,function(set1){
      
      #set1 <- GWASCatalogGeneSets[1]
      lapply(DISGENETGeneSets,function(x){
        
        return(sum(set1%in%x))
      })
    }))
    
 
    stat_SCZ <- rbind(stat_SCZ,sapply(GWASCatalogGeneSets_SCZ,length))
    
    stat_SCZ$fullSize <- c(sapply(GWASCatalogGeneSets_SCZ,length),NA)
    

#############################
# add overlap with reported genes

    SCZ_genes <- read.table("/data/akalin/Base/AccessoryData/Fromer_Schizo_CMC_2016/S3_FORMER.txt",
                            he=T,sep = "\t")
    
    SCZ_gene <- as.character(unique(SCZ_genes$Gene.Symbol))
    
    inSCZ_genes <- sapply(GWASCatalogGeneSets_SCZ,function(x){
      
      sum(unique(x)%in%SCZ_gene)
      
    })
    
    stat_SCZ$paper <- c(inSCZ_genes,length(SCZ_gene))


    
    
#############################
# genes that are confirmed across methods
    con.SCZ_genes <- sapply(GWASCatalogGeneSets_SCZ,function(x){
      
      unique(x)[(unique(x)%in%SCZ_gene)]
      
    })
    table(sort(table(unlist(con.SCZ_genes))))
    # 75% genes is unique 52/70
    
    # how much genes is confirmed by other methods
    
    df <- do.call("rbind.data.frame",lapply(con.SCZ_genes,function(X){
      #
      lapply(con.SCZ_genes,function(Y) return(sum(Y%in%X)))
      
    }))
    #pheatmap(df/apply(df,1,max),display_numbers = T)

    df <- round(100*df/apply(df,1,max))
    
       df.m <- melt(df)
       df.m$Publication <- colnames(df) 
       

       
# pheatmap plot for SCZ condirmed genes
       
   path <- "/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/fig4H_SCZ_Confirmed_Genes_Pheatmap.png"
    png(path,height = 570, width = 870)
    
    ggplot(data = df.m, aes(variable,
                          Publication,
                          fill = value))+
      
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 50, limit = c(0,100), space = "Lab", 
                           name="Overlap ratio (%)") +
      ylab("")+
      xlab("")+
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 40, hjust = 1,face="bold"),
            axis.title = element_text(size=50,face="bold"),
            axis.text.y =element_text(size=40,face="bold"),
            legend.text = element_text(size=28,face="bold"),
            legend.title =  element_text(size=30,face="bold"))+
      geom_text(aes(variable,Publication,  label = value), color = "black", size =13) +
      coord_fixed()
    
    
   dev.off() 
   
   
# pheatmap plot for all SCZ genes annotated using different EGAs  
    
    df <- do.call("rbind.data.frame",lapply(GWASCatalogGeneSets_SCZ,function(X){
      #
      lapply(GWASCatalogGeneSets_SCZ,function(Y) return(sum(Y%in%X)))
      
    }))
    df <- round(100*df/apply(df,1,max))
     df.m <- melt(df)
    df.m$Publication <- colnames(df) 
    
    
    
      path <- "/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/fig4I_SCZ_Genes_Pheatmap.png"
      png(path,height = 570, width = 870)
      
      ggplot(data = df.m, aes(variable,
                              Publication,
                              fill = value))+
        
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 50, limit = c(0,100), space = "Lab", 
                             name="Overlap ratio (%)") +
        ylab("")+
        xlab("")+
        theme_minimal()+ 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 40, hjust = 1,face="bold"),
              axis.title = element_text(size=50,face="bold"),
              axis.text.y =element_text(size=40,face="bold"),
              legend.text = element_text(size=28,face="bold"),
              legend.title =  element_text(size=30,face="bold"))+
        geom_text(aes(variable,Publication,  label = value), color = "black", size =13) +
        coord_fixed()
      
      
      dev.off() 
    
    
# 98% of FOCS genes was present in other methods  as compared to 51% of GeneHancer genes,
# 49% of JEME genes, and 31% EnhancerAtlas SCZ genes
      
    FOCS_SCZ_genes <- sum(GWASCatalogGeneSets_SCZ$FOCS%in%unlist(GWASCatalogGeneSets_SCZ[names(GWASCatalogGeneSets_SCZ)!="FOCS"]))
    FOCS_SCZ_genes/length(GWASCatalogGeneSets_SCZ$FOCS)
    
    
    FOCS_SCZ_genes <- sum(GWASCatalogGeneSets_SCZ$FOCS%in%unlist(GWASCatalogGeneSets_SCZ[names(GWASCatalogGeneSets_SCZ)!="FOCS"]))
    FOCS_SCZ_genes/length(GWASCatalogGeneSets_SCZ$FOCS)
    
    
   confirmedGenes <- sapply(names(GWASCatalogGeneSets_SCZ),function(x){
      
      res <- sum(GWASCatalogGeneSets_SCZ[[x]]%in%unlist(GWASCatalogGeneSets_SCZ[names(GWASCatalogGeneSets_SCZ)!=x]))
      return(100*res/length(GWASCatalogGeneSets_SCZ[[x]]))
      
    })
    


   


###################
#------CAD--------#
###################

##############################
# search GWAS

GWASCatalogGeneSets_CAD <- getGWASCatalogGeneSets(terms=c("EFO:0001645"),
                                              GWASCatalog=GWASCatalog,
                                              SNPannotations=GWAS_PATHS)

    GWASCatalogGeneSets_CAD <- lapply(GWASCatalogGeneSets_CAD,
                              function(x){return((unique(x)))})


    GWASCatalogGeneSets_CAD <- GWASCatalogGeneSets_CAD[names(GWASCatalogGeneSets_CAD)%in%
                                                         c("FOCS","JEME","EnhancerAtlas","GeneHancer")]

#############################
# how many genes in DISGENE for 1) CRC and 2) CRC ancestor terms
# NOTE: UMLS MANUALLY SEARCHED: https://www.ebi.ac.uk/ols/ontologies/efo/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0001645
DISGENETGeneSets_CAD <- getDISGENETGeneSets(efo.terms =c("EFO:0001645"),
                                           umls.terms = "C1956346",
                                        DISGENET=DISGENET)

    DISGENETGeneSets_CAD <- lapply(DISGENETGeneSets_CAD,unique)

sapply(DISGENETGeneSets_CAD,length)





stat_CAD <- do.call("rbind.data.frame",lapply(GWASCatalogGeneSets_CAD,function(set1){
  
  #set1 <- GWASCatalogGeneSets[1]
  lapply(DISGENETGeneSets_CAD,function(x){
    
    return(sum(set1%in%x))
  })
}))


stat_CAD <- rbind(stat_CAD,sapply(GWASCatalogGeneSets_CAD,length))

stat_CAD$fullSize <- c(sapply(GWASCatalogGeneSets_CAD,length),NA)



##############################
# add overlap with reported genes

CAD_genes <- unlist(read.table("/data/akalin/Base/AccessoryData/Mäkinen_2014_PlosGenetics_CAD/CAD_genes.txt",sep=",",stringsAsFactors = F))
inCAD_genes <- sapply(GWASCatalogGeneSets_CAD,function(x){
  
  sum(unique(x)%in%c(CAD_genes))
  
})

totalGenes <- sapply(GWASCatalogGeneSets_CAD,length)
stat_CAD$paper <- c(inCAD_genes,length(CAD_genes))

stat_CAD[,c("DISGENET_AncestorGenes","CAD_paper")]



stat_SCZ$Disease <- "SCZ"
stat_CRC$Disease <- "CRC"
stat_CAD$Disease <- "CAD"



Overlap <- cbind("Publication"=rownames(stat_SCZ),
                      rbind(stat_SCZ,
                            stat_CRC,
                            stat_CAD))

    percentageGenesConirmed <- round(100*Overlap$DISGENET_AncestorGenes/Overlap$fullSize,2)
    percentageGenesConirmed_paper <- round(100*Overlap$paper/rep(c(650,63,191),each=5),2)
    
    Overlap$GenesConf_DISGENET <-percentageGenesConirmed
    Overlap$GenesConf_publ <-percentageGenesConirmed_paper


Overlap.N <- melt(cbind("Publication"=rownames(stat_SCZ)[-5],
                      Overlap [Overlap$Publication!="5",c("DISGENET_AncestorGenes",
                         "paper",
                         "Disease")]))

    Overlap.N.paper <- Overlap.N[Overlap.N$variable=="paper",]
    Overlap.N.DISGENET <- Overlap.N[Overlap.N$variable!="paper",]

Overlap.perc <- melt(cbind("Publication"=rownames(stat_SCZ)[-5],
                        Overlap [Overlap$Publication!="5",
                                 c("GenesConf_DISGENET",
                                   "GenesConf_publ",
                                   "Disease")]))

    Overlap.perc.paper <- Overlap.perc[Overlap.perc$variable=="GenesConf_publ",]
    Overlap.perc.DISGENET <- Overlap.perc[Overlap.perc$variable!="GenesConf_publ",]


###########

png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/Bench_perc_GWAS_PUBL_genes_CRC_SCZ_CAD.png",
    width = 500,height = 400)


      colorScheme <- c("black","darkorange","blue","darkcyan")
      names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")
      
      
      ggplot(Overlap.perc.paper,aes(x=Disease,
                         y=as.integer(as.character(value)),
                         fill=Publication,order=Disease), color=colorScheme) +  
        stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
        scale_fill_manual(values =  alpha(colorScheme[sort(names(colorScheme))], .9))+
        ylab("% [conf. genes] \n") +
        xlab("\n Disease") + 
        ggtitle("Selected publications") + 
        #ylim(0,85) +
        theme(text = element_text(size=30,face="bold"),
              axis.text =element_text(size=30,face="bold"),
              #axis.text.x=element_blank(),
              #axis.text.y=element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.position="none",
              axis.line = element_line(colour = "black"),
              legend.text = element_text(size=25,face="bold"))
      geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")
      
      
      dev.off()




######
png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/Bench_perc_GWAS_DISGENET_genes_CRC_SCZ_CAD.png",
    width = 500,height = 400)


colorScheme <- c("black","darkorange","blue","darkcyan")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")


ggplot(Overlap.perc.DISGENET,aes(x=Disease,
                              y=as.integer(as.character(value)),
                              fill=Publication,order=Disease), color=colorScheme) +  
  stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
  scale_fill_manual(values =  alpha(colorScheme[sort(names(colorScheme))], .9))+
  ylab("% [conf. genes] \n") +
  xlab("\n Disease") + 
  #ylim(0,85) +
  ggtitle("\t DisGeNET") + 
  theme(text = element_text(size=30,face="bold"),
        axis.text =element_text(size=30,face="bold"),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=25,face="bold"))
geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()










###########

png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/Bench_N_GWAS_PUBL_genes_CRC_SCZ_CAD.png",
    width = 500,height = 400)


colorScheme <- c("black","darkorange","blue","darkcyan")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")


ggplot(Overlap.N.paper,aes(x=Disease,
                              y=as.integer(as.character(value)),
                              fill=Publication,order=Disease), color=colorScheme) +  
  stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
  scale_fill_manual(values =  alpha(colorScheme[sort(names(colorScheme))], .9))+
  ylab("N [conf. genes] \n") +
  xlab("\n Disease") + 
  ggtitle("Selected publications")+
  #ylim(0,85) +
  theme(text = element_text(size=30,face="bold"),
        axis.text =element_text(size=30,face="bold"),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=25,face="bold"))
geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()




######
png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/Bench_N_GWAS_DISGENET_genes_CRC_SCZ_CAD.png",
    width = 500,height = 400)


colorScheme <- c("black","darkorange","blue","darkcyan")
names(colorScheme) <- c("FOCS","GeneHancer","JEME","EnhancerAtlas")


ggplot(Overlap.N.DISGENET,aes(x=Disease,
                                 y=as.integer(as.character(value)),
                                 fill=Publication,order=Disease), color=colorScheme) +  
  stat_summary(fun.y=mean,position="dodge",geom="bar",width = 0.5) +  
  scale_fill_manual(values =  alpha(colorScheme[sort(names(colorScheme))], .9))+
  ylab("N [conf. genes] \n") +
  xlab("\n Disease") + 
  ggtitle("\t DisGeNET") + 
  #ylim(0,85) +
  theme(text = element_text(size=30,face="bold"),
        axis.text =element_text(size=30,face="bold"),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=25,face="bold"))
geom_text(aes(label = value), size = 15, hjust = 1, vjust = 1, position = "stack")


dev.off()




#############################
# pheatmap of confirmed genes -CAD
#############################
# genes that are confirmed
con.CAD_genes <- sapply(GWASCatalogGeneSets_CAD,function(x){
  
  unique(x)[(unique(x)%in%CAD_genes)]
  
})

table(sort(table(unlist(con.CAD_genes))))
# 75% genes is unique 52/70

# how much genes is confirmed by other methods

df <- do.call("rbind.data.frame",lapply(con.CAD_genes,function(X){
  #
  lapply(con.CAD_genes,function(Y) return(sum(Y%in%X)))
  
}))
pheatmap(df/apply(df,1,max),display_numbers = T)




# full set of genes associated with CAD - overlap statistics
df <- do.call("rbind.data.frame",lapply(GWASCatalogGeneSets_CAD,function(X){
  #
  lapply(GWASCatalogGeneSets_CAD,function(Y) return(sum(Y%in%X)))
  
}))
pheatmap(df/apply(df,1,max),display_numbers = T)


##------------
df <- do.call("rbind.data.frame",lapply(con.SCZ_genes,function(X){
  #
  lapply(con.SCZ_genes,function(Y) return(sum(Y%in%X)))
  
}))
pheatmap(df/apply(df,1,max),display_numbers = T)






df <- do.call("rbind.data.frame",lapply(GWASCatalogGeneSets_SCZ,function(X){
  #
  lapply(GWASCatalogGeneSets_SCZ,function(Y) return(sum(Y%in%X)))
  
}))
pheatmap(df/apply(df,1,max),display_numbers = T)










##--------------------
# how much do SNP-genes generally overlap?



# GWAS_PATHS <- "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/GWAS_Catalog/annotatedGWAS_EGA/"
# 
# 
# gw <- list.files(GWAS_PATHS,full.names = T) [c(1,3,4,6)]
# EA <- readRDS(gw[1])
# FOCS <- readRDS(gw[3])
# GH <- readRDS(gw[4])
# JEME <- readRDS(gw[6])
# head(JEME)
# 
# FindOverlap.gwas <- function(query,
#          subject){
# 
#   q.s <- data.frame(linkOverlaps(query,
#                         anchorOne(subject),
#                         anchorTwo(subject)))
# 
#   q.c.s <- unique(q.s$query[q.s$subject1==q.s$subject2])
#   # query confirmed by subject
#   
#  return(q.c.s) 
#   
# }

