# plot pheatmap of the mean coverage of ChromHMM defined chromatin state by different enhancer definitions
# step 1. for chromatin state calculate across cell type coverage
# step 2. plot

# this analysis is based on...

library(stringr)


chrhmm <- list.files(outdir <- "/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Mnemomics/ChromHMM_regions/",
  pattern="EnhChrHMM",full.names=T)




# step 1. for chromatin state calculate across cell type coverage

      averageCoverage <- do.call("cbind.data.frame",lapply(chrhmm, function(x){
        # across cell type coverage
        
        print(x)
        
        #x=chrhmm[1]
       
        perEnhRes <- readRDS(x)
        
        means <- colMeans(perEnhRes)
        
        
      }))
      
      
      
      colnames(averageCoverage) <- c("IHE","EnhancerAtlas","ConsensusE","FOCS",
                                     "GeneHancer","JEME")
      
      
      # remove enhancer check
      averageCoverage <- averageCoverage[!rownames(averageCoverage)=="Enhancer",]
      
      # ordering 
      averageCoverage <- averageCoverage[order(as.integer(str_extract(rownames(averageCoverage),"[0-9]*"))),]
      averageCoverage <- round(averageCoverage*100,0)
      #rownames(averageCoverage) <- str_replace(rownames(averageCoverage),"[0-9]*_","")


      ###################################
      # png for 4 publications
      # png("/data/akalin/Projects/AAkalin_Catalog_RI/AAkalin_CatalogRI/Results/Plots/PAPER/ChrommHMMcoverageEnhancers_4publ.png", width = 900, height = 750)
      # 
      # 
      # 
      # pheatmap::pheatmap((averageCoverage[c("FOCS","JEME","GeneHancer","EnhancerAtlas")]),
      #                    display_numbers = T,
      #                    number_format = "%.2f",  
      #                    fontsize = 28,
      #                    cluster_rows =  F,
      #                    cluster_cols = F)
      # 
      # dev.off()


     
     
      
      
      # Melt the correlation matrix
      library(reshape2)
      melted_cormat <- melt(averageCoverage[,c("FOCS","JEME","GeneHancer","EnhancerAtlas")])
          melted_cormat$ChromHMM <- rownames(averageCoverage)
          
      
          melted_cormat <- transform(melted_cormat, 
                                     ChromHMM  = factor(
                                ChromHMM ,
                                levels=unique(melted_cormat$ChromHMM),
                                ordered =TRUE))
          
          
      # Heatmap
      library(ggplot2)
      
      path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/fig1G_ChromHMMcoverage.png"
      
      png(path,height =1100, width = 1100)
      
      ggplot(data = melted_cormat, aes(variable , ChromHMM, fill = value))+
        
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 50, limit = c(0,100), space = "Lab", 
                             name="Coverage \n ChromHMM (%)") +
        ylab("")+
        xlab("")+
        theme_minimal()+ 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 40, hjust = 1,face="bold"),
              axis.title = element_text(size=50,face="bold"),
              axis.text.y =element_text(size=40,face="bold"),
              legend.text = element_text(size=28,face="bold"),
              legend.title =  element_text(size=30,face="bold"))+
        geom_text(aes(variable, ChromHMM, label = value), 
                  color = "black", size =10) +
        coord_fixed()
      
      
      dev.off()




      
      
      
      # Heatmap
      library(ggplot2)
      
      path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/fig1G_ChromHMMcoverage_vhoriz.png"
      
      png(path,height =400, width = 1100)
      
      ggplot(data = melted_cormat, aes( ChromHMM,variable, fill = value))+
        
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 50, limit = c(0,100), space = "Lab", 
                             name="Coverage \n ChromHMM (%)") +
        ylab("")+
        xlab("")+
        theme_minimal()+ 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 30, hjust = 1,face="bold"),
              axis.title = element_text(size=50,face="bold"),
              
              axis.text.y  =element_text(size=40,face="bold"),
              legend.text = element_text(size=28,face="bold"),
              legend.title =  element_text(size=30,face="bold"))+
        geom_text(aes(ChromHMM,variable, label = value), 
                  color = "black", size =9) +
        coord_fixed()
      
      
      dev.off()
      
      
      
      
      melted_cormat$ChromHMM <- str_replace(melted_cormat$ChromHMM,'_'," ")
      melted_cormat <- transform(melted_cormat, 
                                 ChromHMM  = factor(
                                   ChromHMM ,
                                   levels=unique(melted_cormat$ChromHMM),
                                   ordered =TRUE))
      
      # Heatmap
      library(ggplot2)
      
      path="/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Downstream/Plots/figure1_publicationE/fig1G_ChromHMMcoverage_vhoriz2.png"
      
      png(path,height =400, width = 1100)
      
      ggplot(data = melted_cormat, aes( ChromHMM,variable, fill = value))+
        
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 50, limit = c(0,100), space = "Lab", 
                             name="Coverage \n ChromHMM (%)") +
        ylab("")+
        xlab("ChromHMM chromatin states")+
        theme_minimal()+ 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 30, hjust = 1,face="bold"),
              axis.title = element_text(size=50,face="bold"),
              
              axis.text.y  =element_text(size=40,face="bold"),
              legend.text = element_text(size=28,face="bold"),
              legend.title =  element_text(size=30,face="bold"))+
    
        coord_fixed()
      
      
      dev.off()