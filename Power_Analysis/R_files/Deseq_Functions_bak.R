####################Computing the Standard Errors############################
Standard_Errors <- function(logcontrol,logtreatment, param.disp){
  a1 =length(logcontrol); b1 =length(logtreatment) 
  indx1 = which(is.finite(logcontrol))
  logcontrol = logcontrol[indx1]; logtreatment = logtreatment[indx1]
  indx2 = which(is.finite(logtreatment))
  logcontrol = logcontrol[indx2]; logtreatment = logtreatment[indx2]
  
  a2 =length(logcontrol); b2 =length(logtreatment) 
  #print(c("logcontrol_non_zeros_prop =",a2/a1) )
  #print(c("logtreatment__non_zeros_prop=",b2/b1))
  
  control = 2^logcontrol;  treatment =  2^logtreatment
  c0 = param.disp[["asymptDisp"]]; c1 = param.disp[["extraPois"]]
  ##################################################################
  disp1= c0 + c1/control 
  disp2= c0 + c1/treatment 
  Var_cont = control + (control^2)*disp1
  Var_treat = treatment + (treatment^2)*disp2
  cont_deriv=  1/(control*log(2)); treat_deriv = 1/(treatment*log(2))
  ##################################################################
  Var_lcont = cont_deriv*Var_cont
  Var_ltreat = treat_deriv*Var_treat
  V = Var_lcont+Var_ltreat
  p=list(sqrt(V),a2/a1,b2/b1)
  names(p) = c(sd_error,logcontrol_non_zeros_prop,logtreatment__non_zeros_prop)
  p
}

######### data
##' @param data_sim function to simulate data
##' @param param.effect parameters from fitting effects
##' @param param.disp   dispersion parameters
##' @param param.cont   parameters from fitting control 
##' @param n_sim        number of simulations
##' @param n_samples    number of samples
##' @param n_otu        number of otu
##' @param filter       filtering threshold

data_sim<-function(param.effect,param.disp,param.cont,n_sim,n_samples,n_otu,name,filter=10){
  
  Sim_data = Sim_metadata = dispersions =  True_data = SD_Error= list()
  true_treatment = true_effect =  true_control = list()
  
  new_otu <- 2*n_otu
  mu= param.cont$mu; sigma = param.cont$sigma; lambda = param.cont$lambda
  for(i in 1:n_sim){
    
    #simulate control and effects 
    myMix <- UnivarMixingDistribution(Norm(mean=mu[1], sd =sigma[1]),
                                      Norm(mean=mu[2], sd =sigma[2]),
                                      mixCoeff=c(lambda[1],lambda[2]))
    
    rmyMix <- r(myMix)
    
    control <- rmyMix(new_otu)
    cont <- unlist(control)
    
    ## scale parameters values for cauchy 
    scal = data.frame(x =  exp(param.effect[1]+param.effect[2]*cont + 
                                 param.effect[3]*cont^2))
    
    effect <-apply(scal,1,function(x){
      rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})
    
    ## treatment simulation
    eff <- unlist(effect)
    treatment <-   eff + cont; treat <- unlist(treatment)
    
    mean_abund <- data.frame(treat = 2^treat, cont = 2^cont)
    mean_abund.data <- data.frame(means = rowMeans(mean_abund))
    
    dispers <- function(mean,c0,c1){c0 + c1/mean}
    disp <- drop(1/apply(mean_abund.data, 2, dispers,  c0 = param.disp[["asymptDisp"]],
                         c1 = param.disp[["extraPois"]]))
    
    
    dispersions <- disp
    df <- data.frame(mean_abund.data,disp)
    
    ## Simulate data
    data <-matrix(apply(df,1,function(x,n=n_samples){rnbinom(n=n,mu=x[1],size=x[2])}),
                  nrow = new_otu,ncol = n_samples)
    
    ## Filter data 
    keep <- rowSums(data) >= filter
    data <- data[keep,]
    effect <- effect[keep]; control <- control[keep]
    treatment <- treatment[keep]; dispersions <- disp[keep]
    
    #Select n_otus after filtering
    select <- sample.int(n_otu)
    data <- data[select,]
    true_effect[[i]] <- effect[select]; true_control[[i]] <- control[select]
    true_treatment[[i]] <- treatment[select]; dispersions <- disp[select] 
    
    logcont =  control[select]; logtreat = treatment[select]
    
    SD_Error[[i]] = Standard_Errors(logcont,
                               logtreat, 
                               param.disp)
      
    
    Sim_metadata[[i]] <- data.frame(Groups = sort(rep(c("ASD","NT"),each = n_samples/2)),
                                    row.names=paste0("Sample",1:n_samples))
    
    rownames(data) <- paste0("otu",1:n_otu)
    colnames(data) <- paste0("Sample",1:n_samples)
    
    Sim_data[[i]] <- data
  } 
  names(true_effect) = names(true_control) = names(true_treatment) =  paste0("Sim_data",1:n_sim)
  names(dispersions) =  names(Sim_metadata) = names(Sim_data) = paste0("Sim_data",1:n_sim)
  names(SD_Error) = paste0("Sim_data",1:n_sim)
  
  save(SD_Error, file = paste0(path,"Standard_error_", name, ".RData")) 
  #save(dispersions, file = paste0(path,"Dispersions_", name, ".RData")) 
  save(true_control, file = paste0(path,"True_Control_", name, ".RData")) 
  save(true_effect, file = paste0(path,"True_Effect_", name, ".RData")) 
  save(true_treatment, file = paste0(path,"True_Treatment_", name, ".RData")) 
  
  save(Sim_metadata, file = paste0(path,"Simulated_Metadata_", name, ".RData")) 
  save(Sim_data, file = paste0(path,"Simulated_Data_", name, ".RData")) 
}


Deseq <-function(data_list,metadata_list,name,path){
  result  = list()
  index = 1
  
  for(n in names(data_list)){
    data = data_list[[n]]; metadata <- metadata_list[[n]]
    dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
    
    r <- which(colSums(data) == 0)
    if (length(r) > 0){
      data <- data[,-r]; metadata <- metadata[-r,]
      dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
    }
    dds$Groups <- relevel(dds$Groups, ref = "NT")
    dds <- DESeq(dds,sfType ="poscounts", minReplicatesForReplace=Inf) 
    res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
    reslt <- lfcShrink(dds, res=res, coef=2, type="normal")  
    result[[index]] <- data.frame(reslt)
    index = index + 1
  }
  names(result) <- paste0("result",1:n_sim)
  save(result, file = paste0(path, "Sim_Deseq_Results_",name,".RData")) 
} 
