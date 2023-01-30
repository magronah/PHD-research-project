Control_Treat_Split <- function(data_list, metadata_list ){
  
  index <- 1; L =treat = list(); control = treat_mean = cont_mean <- list()
  for(n in names(metadata_list)){
    
    data <- data_list[[n]];  meta <- metadata_list[[n]]
    m <-  rownames(data)[1]; p <- colnames(data)[1]
    
    if (nchar(m) < nchar(p) ) {
      data <- t(data)
    } else {
      data <- (data)
    }
    
    asd <- which(meta$Groups == "ASD" )
    treat_data <- data[asd,]; control_data <- data[-asd,]
    stopifnot(rownames(treat_data) == meta$Samples[asd])
    stopifnot(rownames(control_data) == meta$Samples[-asd])
    
    L[[index]] <-  (list(treatment = treat_data, control =control_data))
    #cont_mean[[index]] <- colMeans(control_data)
    #treat_mean[[index]] <- colMeans(treat_data)
    #effects <- log2(colMeans(treat_mean)/colMeans(treat_mean))  
    index <- index + 1
    
  }
  names(L) <- names(metadata_list)
  return(L)
}

########################################################################

deseq_res<- function(data_list, metadata_list,plot_path, plotname, filter=10){
  
  DeseqResults <- list(); Control_Treat  <- list(); index <- 1
  
  for(n in names(data_list)){
    
    data  = data_list[[n]]; metadata = metadata_list[[n]]
    
    ##Check data is taxa  by samples
    m<-  rownames(data)[1]; p <- colnames(data)[1]
    if (nchar(m) < nchar(p) ) {
      data <- data
    } else {
      data <- t(data)
    }
    
    ##Deseq processes
    dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
    
    keep <- rowSums(counts(dds)) >= filter
    dds <- dds[keep,]; data <- data[keep,] 
    
    #' Remove samples with all zero (deseq needs this to run effectively) 
    r <- which(colSums(data) == 0)
    if (length(r) > 0){
      data <- data[,-r]; metadata <- metadata[-r,]
      dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
    }
    dds$Groups <- relevel(dds$Groups, ref = "NT")
    
    ddsDE <- DESeq(dds, sfType ="poscounts",minReplicatesForReplace=Inf) #4 was good 
    res <- results(ddsDE, cooksCutoff=TRUE, independentFiltering=TRUE) 
    
    reslt <- lfcShrink(ddsDE, res = res, coef="Groups_ASD_vs_NT", type="apeglm")
    result <- data.frame(reslt)
    DeseqResults[[index]] <- result
    
    
    file_name = paste0(plot_path, "DeseqplotMA_", n , ".pdf", sep="")
    pdf(file_name)
    DESeq2::plotMA(reslt, main="normal")
    dev.off()
    
    ##save boxplots
    file_name = paste0(plot_path, plotname, n , ".pdf", sep="")
    pdf(file_name)
    boxplot(log10(assays(ddsDE)[["cooks"]]), range=0, las=2)
    dev.off()
    
    ##Control Treatment splits
    k <- which(metadata$Groups == "ASD")
    control <- data[,-k]; treatment <- data[,k] 
    stopifnot(colnames(control) == metadata$Samples[-k])
    stopifnot(colnames(treatment) == metadata$Samples[k])
    
    cont_Means <- rowMeans(control); treat_Means <- rowMeans(treatment)
    logcontrol = log2(cont_Means);   logtreat = log2(treat_Means)
    
    Control_Treat[[index]] <- list(Control_Data = control, Control_Means = cont_Means,
                                   LogControl  =  logcontrol, Treat_Data =  treatment,
                                   Treat_Means  = treat_Means,
                                   LogTreat  = logtreat)
    
    index <- index + 1
    
  }
  names(Control_Treat) <- names(data_list)
  names(DeseqResults) <- names(data_list)
  
  save(DeseqResults, file = paste0(path, "DeseqResults.RData"))
  save(Control_Treat, file = paste0(path, "Control_Treat_Data.RData"))
  
}

########################################################
read_data <- function(Dataset_Lists, Extract){
  
  Extract_data <- list(); index <- 1
  for(n in names(Dataset_Lists)){
    data <- Dataset_Lists[[n]]
    data <- data[[Extract]]
    Extract_data[[index]] <-  data
    index <- index + 1
  }
  names(Extract_data) <- names(Dataset_Lists)
  return(Extract_data)
} 

########################################################
Plot_ContEffects <-function(deseq_list,logcontrol,logtreatment, plot_path){
  
  stopifnot(names(deseq_list) == names(logcontrol))
  stopifnot(names(deseq_list) == names(logtreatment))
  
  L = CT = list(); original <- c(); save_effect <- list() 
  index <- 1
  for(n in names(deseq_list)){  
    
    res<-deseq_list[[n]]
    logcont <- logcontrol[[n]];  logtreat <- logtreatment[[n]];  
    foldchange <- res$log2FoldChange
    
    #original[index] <- length(logcont)
    #inf[index] <- length( which( is.infinite(logcont)) )
    
    finit <- which(is.finite(logcont))
    foldchange <- foldchange[finit]
    control <- logcont[finit]
    treatment <- logtreat[finit]
    save_effect[[index]] <- foldchange
    
    CT[[index]] <- data.frame(treatment, control)
    FC <- length(res$log2FoldChange); fltFC <- length(foldchange)
    diff <- FC - fltFC
    L[[index]] <- data.frame(FoldChange= FC,
                             FoldChange_filtered = fltFC, Difference = diff)
    
    data <- data.frame(control, foldchange)
    
    file_name = paste0(plot_path, "control_effects_", n , ".pdf", sep="")
    pdf(file_name)
    plts <-  ggplot(data, aes(control, foldchange)) +
      geom_point()  + 
      geom_smooth() +
      #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
      ggtitle(paste0("Control against Fold Changes")) +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
    
    print(plts)
    dev.off()
    
    print(plts)
    index <- index + 1
  }
  names(L) <- names(deseq_list); names(save_effect) <- names(deseq_list)
  names(CT) <- names(deseq_list)
  names(save_effect) <- names(deseq_list)
  return(list(control_treat =CT,effects =save_effect))
}
########################################################


##' @Transform transforms control and treatment data from log scale to 
##'                 original scale
##'@param DataList list of control and treatment dataframes

Transform <- function(DataList){
  
  index <- 1; tranform = impute =list()
  
  for(n in names(DataList)){
    d <- DataList[[n]]
    
    #replace all infinity with 0 and log scale x to 2^x  
    df <- d %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, 2^x))
    
    tranform[[index]] <- df
    
    #do some imputation  
    impute[[index]] <- df %>% 
      mutate_if(is.numeric, function(x) ifelse(x == 0, 1e-02, x))
    
    index <-  index + 1
  }
  
  names(tranform) <- names(DataList)
  names(impute) <- names(DataList)
  return(list(transform_data = tranform, imputed_data =impute))
}


###################################################################

Impute_plot <- function(Imputed_data, plot_path,plotname){
  
  effects_list <- list(); index <- 1
  
  for(n in names(Imputed_data)){
    dd <- Imputed_data[[n]]
    
    Fold_change = log2(dd$treatment/dd$control); Control_abundance = log2(dd$control) 
    data <- data.frame(Fold_change, Control_abundance) 
    
    effects_list[[index]] <- Fold_change
    file_name = paste0(plot_path, plotname, n , ".pdf", sep="")
    pdf(file_name)
    plt <- ggplot(data, aes(Control_abundance, Fold_change)) +
      geom_point()  + 
      geom_smooth() +
      #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
      ggtitle(paste0("Control against Fold Changes")) +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
    print(plt)
    dev.off()
    print(plt)
    index <- index + 1
  }
  names(effects_list) <- names(Imputed_data)
  return(effects_list)
} 


##################################################
##' @Plot_ContEffects function to plot control verse effect sizes 
##' @deseq_list list of results from deseq2
##' @logcontrol list of log control data
##' @thresh     filter all effect sizes above this threshold

Plot_ContEffects <-function(deseq_list,logcontrol,logtreatment, plot_path, thresh=20){
  
  stopifnot(names(deseq_list) == names(logcontrol))
  stopifnot(names(deseq_list) == names(logtreatment))
  
  L = CT = list();  index <- 1
  
  for(n in names(deseq_list)){  
    
    res<-deseq_list[[n]]; logcont <- logcontrol[[n]];  logtreat <- logtreatment[[n]]
    
    foldchange <- res$log2FoldChange
    
    #finit <- which(is.finite(logcont))
    #foldchange <- foldchange[finit]
    #control <- logcont[finit]
    #treatment <- logtreat[finit]
    # k <- which(abs(foldchange) < thresh)
    # foldchange <- foldchange[k]
    # control <- control[k]
    #save_effect[[index]] <- foldchange
    #print(length((which(is.infinite(control)))))
    #CT[[index]] <- data.frame(treatment, control)
    #FC <- length(res$log2FoldChange); fltFC <- length(foldchange); diff <- FC - fltFC
    #L[[index]] <- data.frame(FoldChange= FC,
    #           FoldChange_filtered = fltFC, Difference = diff)
    
    data <- data.frame(control, foldchange)
    file_name = paste0(plot_path, "control_effects_", n , ".pdf", sep="")
    pdf(file_name)
    plts <-  ggplot(data, aes(control, foldchange)) +
      geom_point()  + 
      geom_smooth() +
      #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
      ggtitle(paste0("Control against Fold Changes")) +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
    
    print(plts)
    dev.off()
    
    print(plts)
    index <- index + 1
  }
  names(L) <- names(deseq_list)
  names(save_effect) <- names(deseq_list)
  names(CT) <- names(deseq_list)
  names(save_effect) <- names(deseq_list)
  
  return(list(control_treat =CT,effects =save_effect))
}

