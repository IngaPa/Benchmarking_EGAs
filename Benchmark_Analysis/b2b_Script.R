# ---
#   title: "ComparingEP_links"
# output: html_document
# ---
#   
#   

#```{r setup, include=FALSE}
library(stringr)
library(pheatmap)
library(GenomicRanges)
library(GenomicInteractions)
library(InteractionSet)
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/downstreamFunctions.R") # extract_b2b is saved here!
source("/data/akalin/Projects/AAkalin_reg2gene/reg2gene/R/projectHelpFunctions.R") # extract_and_plot_Bench iis here
source("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Scripts/Downstream/plots/generalPlot_functions.R")


#```


#```{r}

benchResults <- sort(list.files("/data/akalin/Projects/AAkalin_Catalog_RI/final_Version/Results/Benchmarking", 
                                "Val.rds",full.names = T))

benchResults <- benchResults[str_detect(benchResults,"JEME|GeneHancer|EnhancerAtlas|FOCS|GTEx|Westra|HiC|CCSI")]
benchList <- lapply(benchResults,readRDS) # the last one is full reg2gene object

names(benchList)=str_replace(basename(benchResults),"BenchVal2BenchVal.rds","")




#```

#How many EP links across different publications:
  
  #```{r}
sapply(benchList,length)

sapply(benchList,function(x)length(unique(x)))


#Getting results of beench2bench analysis
#```{r , fig.height=15, fig.width=15}




extractB2B <- function(benchList,
                       names=c("JEME","GeneHancer","EnhancerAtlas","FOCS","GTEx","PCHiC","WestraeQTLs","CCSI"),
                       names2=c("JEME","GeneHancer","EnhancerAtlas","FOCS","GTEx","PC HiC","Westra eQTLs","CCSI")){
  
  require(stringr)
 
 
  matrixEd <-  lapply(names,function(x){
    #x="CCSI"
    #names=paste("Filter",names(benchList),sep="_")
    print(x)
    # keeping function open for both non-filtered and filtered datasets
    
        # adjusting names 
  #  if (any(str_detect(names,'Filter'))){ names=str_replace(str_replace(names,' ',"."),"-",".")}
    
    metadata <- mcols(benchList[[which(names(benchList)==x)]])
    metadata <-  as.data.frame(metadata[which(colnames(metadata)%in%names2)])
    metadata.binary <- apply(metadata,2,function(y)y!=0)
    
    return(colSums(metadata.binary))
    
  })
  
  matrixEd.df <- do.call("rbind.data.frame",matrixEd)
  rownames(matrixEd.df) <- names
  colnames(matrixEd.df) <-  names(matrixEd[[1]])
  # check if column naming is ok!
  
  return(matrixEd.df)
}



b2b <- extractB2B(benchList,
                  names=c("JEME","GeneHancer","EnhancerAtlas","FOCS","GTEx","PCHiC","WestraeQTLs","CCSI"))

b2filter <- extractB2B(benchList,
                  names2=paste("Filter",
                        names2 <- c("JEME","GeneHancer","EnhancerAtlas","FOCS","GTEx","PC.HiC","Westra.eQTLs","CCSI"),
                        sep="_"),
                  names=c("JEME","GeneHancer","EnhancerAtlas","FOCS","GTEx","PCHiC","WestraeQTLs","CCSI"))


# # arranging names
b2b <- b2b[order(rownames(b2b)),]
b2b <- b2b[,order(colnames(b2b))]
colnames(b2b) <-  rownames(b2b) 

b2filter <- b2filter[order(rownames(b2filter)),]
b2filter <- b2filter[,order(colnames(b2filter))]

percentagers <- round(100*b2b/b2filter,2)


png("/data/akalin/Projects/AAkalin_Catalog_RI/post_PhD_analyses/Results/plots/SupplemTable_Benchmark.png", 
    height=200, width=575)
library(gridExtra)
p<-tableGrob(percentagers)
grid.arrange(p)
dev.off()





