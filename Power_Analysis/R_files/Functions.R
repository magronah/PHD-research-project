TruncCauchy <- function(x, loc, scale, low, high, log = TRUE) {
  log.a <- log(pcauchy(high, loc, scale) - pcauchy(low, loc, scale))
  ## compute log-density
  logdens <- ifelse(x < low | x > high, 0,
                    dcauchy(x, loc, scale, log = TRUE) - log.a)
  ## return log-density or density depending on 'log' argument
  if (log) logdens else exp(logdens)
}

#############################################################################
contour_plots <-function(combined_data,standardised = FALSE,
                         deseqstandardised = FALSE,
                         interval=c(0.001,0.001,0.001,0.001)){ 
  
  index = 1; v = plts= list()
  for(n in names(combined_data)){
    
    data = combined_data[[n]]
    comb1 <- with(data,tibble(lcontrol = true_control, abs_lfc = abs(true_effect),
                              power = as.numeric(power), meancounts = 2^(true_control)))
    comb2 <- with(data,tibble(lcontrol = true_control, 
                              abs_lfc = abs(true_effect/deltaStandard_error),
                              power = as.numeric(power), meancounts = 2^(true_control) ))
    
    comb3 <- with(data,tibble(lcontrol = true_control, 
                              abs_lfc = abs(true_effect/deseqStandard_error),
                              power = as.numeric(power), meancounts = 2^(true_control) ))
  
    if (standardised) comb = comb2 
    else if (deseqstandardised) comb = comb3
    else comb = comb1
    
    gg_2dim <- (ggplot(comb)
                + aes(lcontrol, abs_lfc)
                + geom_point(aes(color = power), alpha = 0.1)
                + labs(x = "log(control abundance)",
                       y = "absolute(log(fold change))")
    )
    
    fit_2d <- bam(power ~ te(lcontrol, abs_lfc),
                  data = comb,
                  family = binomial) 
    
    pp <- with(comb,
               expand.grid(lcontrol = seq(min(lcontrol),
                                              max(lcontrol),
                                              length = 50),
                           abs_lfc = seq(min(abs_lfc),
                                         max(abs_lfc),
                                         length = 50)))
    
    pp$power <- predict(fit_2d, newdata = pp,
                        type = "response")
    ################################################################################
    int = interval[index];     brkvec <- seq(0,1,int)
    
    gg_2dimc <- (gg_2dim
                  + geom_contour(data = pp,
                                 aes(z=power),
                                 breaks = brkvec)
                  + geom_label_contour(data = pp, aes(z= power),
                                       breaks = brkvec) + 
                   ggtitle(n) +
                   theme(legend.position = "bottom")
    )
    plts[[index]] = list(point_plot = gg_2dim,
                         tibble_data=pp,
                         bam_predict_plot = fit_2d,
                         contour_plot = gg_2dimc)
    index = index + 1
  }
  names(plts)  = names(combined_data)
  plts
}


#########Function to plot heatmap for power#####################
Power_Heatmap <-function(all_data_list, nTaxa =  100, 
                         standardised = FALSE,
                         deseqstandardised = FALSE,
                         xblocks=8, yblocks=8,showText=TRUE,txtSize=3,
                         heatmap.low="lightgreen",heatmap.high="orangered")
{
  index = 1; plts1 = plts2 = list()
  for(n in names(all_data_list)){
    
    data =  all_data_list[[n]]
    
    comb1 <- with(data,tibble(lcontrol = true_control, abs_lfc = abs(true_effect),
                             power = as.numeric(power), meancounts = 2^(true_control) ))
    
    
    comb2 <- with(data,tibble(lcontrol = true_control, 
                              abs_lfc = abs(true_effect/deltaStandard_error),
                             power = as.numeric(power), meancounts = 2^(true_control) ))
    
    comb3 <- with(data,tibble(lcontrol = true_control, 
                              abs_lfc = abs(true_effect/deseqStandard_error),
                              power = as.numeric(power), meancounts = 2^(true_control) ))
    
    if (standardised) comb = comb2 
    else if (deseqstandardised)  comb = comb3
    else comb =comb1
    # Call ddply to roll-up the data and calculate summary means, etc
    dfe.plot<-ddply(comb,
                    .(logControl=cut(comb$lcontrol,xblocks),
                      logFoldChange=cut(comb$abs_lfc,yblocks)),
                    summarize,
                    Power=mean(power),
                    Sum = sum(meancounts) )
    
    ## BUILD THE SUMMARY CHART
    g1<-ggplot(dfe.plot) +
      geom_raster(aes(logControl,logFoldChange,fill=Power),alpha=0.75) +
      scale_fill_gradient(low=heatmap.low, high=heatmap.high) +
      theme_bw() + theme(axis.text.x=element_text(angle=-90)) +
      ggtitle(paste(n,
                    #" ",
                    #xblocks,
                    #" X ",
                    #yblocks,
                    # " grid of Data \nbetween ( ",
                    # round(min(comb$lcontrol),3),
                    # " : ",
                    # round(min(comb$abs_lfc),3),
                    # " ) and ( ",
                    # round(max(comb$lcontrol),3),
                    # " : ",
                    # round(max(comb$abs_lfc),3),
                    # " )\n\n",
                    sep="")) +
      theme(legend.position = "bottom")
    
    if(showText)g1<-g1+geom_text(aes(logControl,logFoldChange,
                                     label=paste("Sum=",round(Sum,0),
                                                 "\nPow=",round(Power,3)),
                                     fontface=c("italic")),
                                 color="black",size=txtSize)
    ##############################################################################
    dfe.plot  = dfe.plot[dfe.plot$Sum >=nTaxa, ]
    
    g2<-ggplot(dfe.plot) +
      geom_raster(aes(logControl,logFoldChange,fill=Power),alpha=0.75) +
      scale_fill_gradient(low=heatmap.low, high=heatmap.high) +
      theme_bw() + theme(axis.text.x=element_text(angle=-90)) +
      ggtitle(paste(n, 
                    #" ",
                    # xblocks,
                    # " X ",
                    # yblocks,
                    # " grid of Data \nbetween ( ",
                    # round(min(comb$lcontrol),3),
                    # " : ",
                    # round(min(comb$abs_lfc),3),
                    # " ) and ( ",
                    # round(max(comb$lcontrol),3),
                    # " : ",
                    # round(max(comb$abs_lfc),3),
                    # " )\n\n",
                    sep="")) +
      theme(legend.position = "bottom")
  
    if(showText)g2<-g2+geom_text(aes(logControl,logFoldChange,
                                     label=paste("Sum=",round(Sum,0),
                                                 "\nPow=",round(Power,3)),
                                     fontface=c("italic")),
                                 color="black",size=txtSize)
    
    plts1[[index]] <- g1 
    plts2[[index]] <- g2 
    
    index <- index + 1
  }
  list(full_plot=plts1,reduced_plot2=plts2)
}


#########Function to extract deseq estimates#####################
combined_data_list <-function(true_control_list, true_effect_list,
                              deseq_data_list,standard_err_delta_list, alpha=0.1){
  
  index = 1;  combined_data = list()
  stopifnot(names(true_effect_list) == names(deseq_data_list))
  
  for(n in names(deseq_data_list)){

    res_data = deseq_data_list[[n]]
    true_cont = unlist(true_control_list[[n]])
    true_eff  =  unlist(true_effect_list[[n]])
    delta_err = unlist(standard_err_delta_list[[n]])
    
    lfc <- (res_data
            %>% setNames(paste0("padjust", 1:n_sim))
            %>% purrr::map_dfr(pull, log2FoldChange) 
    )
    
    lfc = unlist(lfc)
    baseMean <- (res_data
            %>% setNames(paste0("padjust", 1:n_sim))
            %>% purrr::map_dfr(pull, baseMean) 
    )
    baseMean = unlist(baseMean)
    
    lfcSE <- (res_data
                 %>% setNames(paste0("padjust", 1:n_sim))
                 %>% purrr::map_dfr(pull, lfcSE) 
    )
    lfcSE = unlist(lfcSE)
    
    p_adjust <- (res_data
                 ## 'map' won't work without names
                 %>% setNames(paste0("padjust", 1:n_sim))
                 ## map_dfr = run function on each element, combine the results into a data frame
                 %>% purrr::map_dfr(pull, padj) ## 'pull' = tidyverse equivalent of $/[[
    )
    p_val<- unlist(p_adjust)
    
    pow <-(!is.na(p_val) & p_val < alpha) # finding the adjusted pvalues less than alpha
    
    combined_data[[index]] =  tibble(deseqLogFoldChange = lfc, 
                       deseqPvalues = p_val,power=pow,
                       true_control=true_cont,
                       true_effect = true_eff,
                       deseqBaseMean = baseMean,
                       deseqStandard_error = lfcSE,
                       deltaStandard_error = delta_err
    )
    
    index = index + 1
  }
  
  names(combined_data) = names(deseq_data_list)
  save(combined_data, file = paste0(path,"combined_data_list.RData"))
  combined_data
}

####################Computing the Standard Errors############################
Standard_Error_Delta <- function(logcontrol,logtreatment, param.disp){
  
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
  cont_secderiv=  -1/((log(2))*control^2); treat_secderiv = -1/((log(2))*treatment^2)
  ##################################################################
  Var_lcont = (cont_deriv^2)*Var_cont + (1/4)*(cont_secderiv*Var_cont)^2
  Var_ltreat = (treat_deriv^2)*Var_treat + (1/4)*(treat_secderiv*Var_treat)^2
  
  #Its coming from how I am modelling the dispersion. 
  # print(c0)
  # print(c1)
  # Var_lcont = cont_deriv*Var_cont
  # Var_ltreat = treat_deriv*Var_treat
  error = sqrt(Var_lcont+Var_ltreat)
  #p=list(sqrt(V))
  #names(p) = c("sd_error","logcontrol_non_zeros_prop","logtreatment__non_zeros_prop")
  err=log2(error)
  err
  # print(err)
}
######################################################################

##' @Plot_ContEffects_multi_filters function to plot control- effect sizes 
##' for different filtering thresholds
##' @deseq_list list of results from deseq2 for different filtering thresholds
##' @logcontrol list of log control data for different filtering thresholds
##' 
Plot_ContEffects_multi_filters <- function(deseq_list,logcontrol){
  Ind = 1; plts = list(); DAT=list()
  
  for(n in names(deseq_list)){
    
    deseq_FiltThresh <- deseq_list[[n]]; logFiltThresh <- logcontrol[[n]]
    index = 1; data = list()
    for(i in names(deseq_FiltThresh)){
      res<-deseq_FiltThresh[[i]]; logcont <- logFiltThresh[[i]]
      foldchang <- res$log2FoldChange
      finite<-which(is.finite(logcont))
      control <- logcont[finite]; foldchange <- foldchang[finite]
      data[[index]] <- data.frame(control, foldchange)
      index = index + 1
    }
    names(data) <- names(deseq_FiltThresh)
    DAT[[Ind]] <- data
    Ind = Ind + 1
  }
  names(DAT) <- names(deseq_list)
  
  for(j in 1:length(data)){
    #pd = c("none","none","none","bottom")
    plt <- ggplot2::ggplot(DAT[[1]][[j]],aes(x = control, y=foldchange,color="filter=10")) +
      geom_point() +
      geom_point(data = DAT[[2]][[j]],aes(x = control, y=foldchange,color="filter=30")) +
      geom_point(data = DAT[[3]][[j]],aes(x = control, y=foldchange,color="filter=50")) +
      ggtitle(names(data)[j]) +
      theme()
    plts[[j]] <- plt 
  }
  plts
}






######### data
##' @param data_sim function to simulate data
##' @param param.effect parameters from fitting effects
##' @param param.dispersion   dispersion parameters
##' @param param.cont   parameters from fitting control 
##' @param n_sim        number of simulations
##' @param n_samples    number of samples
##' @param n_otu        number of otu
##' @param filter       filtering threshold

data_sim<-function(param.effects,param.dispersions,param.controls,n_sim,n_samples,n_otu,path,filter)
  {
  
  Sim_data = Sim_metadata = sim_data = sim_metadata = list()
  Dispersions = dispersions = Standard_Err_Delta = standard_error_delta = list()
  
  True_treatment = True_effect =  True_control = list()
  true_treatment = true_effect =  true_control = list()
  
  new_otu <- 10*n_otu; index = 1
  
  for(n in names(param.effects))
    { 
      param.effect = param.effects[[n]]; param.disp  = param.dispersions[[n]]
      
      param.cont = param.controls[[n]]; mn = param.cont$mean
      scale = param.cont$scale_param;  slant = param.cont$slant_param

      for(i in 1:n_sim)
        {
        #simulate control and effects 
        control = rsn(n=new_otu, xi = mn, omega = scale, alpha = slant) 
        #control=rtnorm(n=new_otu, mean = mu, sd = exp(logsd), a = lower) 
        #cont <- unlist(control)
        
        ## scale parameters values for cauchy 
        scal = data.frame(x =  exp(param.effect[,1]+param.effect[,2]*control + 
                                     param.effect[,3]*control^2))
        
        effect <-apply(scal,1,function(x){
          rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})

        ## treatment simulation
        ##eff <- unlist(effect)
        treatment <-   effect + control#; treat <- unlist(treatment)
        
        mean_abund <- data.frame(treatment = 2^treatment, control = 2^control)
        mean_abund.data <- data.frame(means = rowMeans(mean_abund))
        
        dispers <- function(mean,c0,c1){c0 + c1/mean}
        disp <- drop(1/apply(mean_abund.data, 2, dispers,  c0 = param.disp[["asymptDisp"]],
                             c1 = param.disp[["extraPois"]]))
        
        df <- data.frame(mean_abund.data,disp)
        
        ## Simulate data
        data <-matrix(apply(df,1,function(x,n=n_samples){rnbinom(n=n,mu=x[1],size=x[2])}),
                      nrow = new_otu,ncol = n_samples)
        
        ## Filter data 
        keep <- rowSums(data) >= filter
        data <- data[keep,]

        effect <- effect[keep]; control <- control[keep]
        treatment <- treatment[keep]; dispersion <- disp[keep]
        
        #Select n_otus after filtering
        select <- sample.int(n_otu)
        data <- data[select,]
        
        true_effect[[i]] <- effect[select]; true_control[[i]] <- control[select]
        true_treatment[[i]] <- treatment[select]; dispersions[[i]] <- dispersion[select] 
        
        sim_metadata[[i]] <- data.frame(Groups = sort(rep(c("ASD","NT"),each = n_samples/2)),
                                        row.names=paste0("Sample",1:n_samples))
        
        logcont =  control[select]; logtreat = treatment[select]
        
        standard_error_delta[[i]] = Standard_Error_Delta(logcont,
                                       logtreat, 
                                       param.disp)
        
        #print(list(logcont,logtreat,param.disp))
        
        rownames(data) <- paste0("otu",1:n_otu)
        colnames(data) <- paste0("Sample",1:n_samples)
        
        sim_data[[i]] <- data
        } 
  
      names(true_effect) = names(true_control) = names(true_treatment) =  paste0("Sim_data",1:n_sim)
      names(dispersions) =  names(sim_metadata) = names(sim_data) = paste0("Sim_data",1:n_sim)
      names(standard_error_delta) = paste0("Sim_data",1:n_sim)
      
      Sim_data[[index]] = sim_data; Sim_metadata[[index]] = sim_metadata 
      True_effect[[index]] =  true_effect; True_control[[index]]= true_control
      True_treatment[[index]] = true_treatment; Dispersions[[index]] = dispersions
      Standard_Err_Delta[[index]] = standard_error_delta
        
      index = index + 1
      }
   
  names(Sim_data) = names(Sim_metadata) = names(Dispersions) =  names(param.effects)
  names(True_effect) =  names(True_control) = names(True_treatment) = names(param.effects)
  names(Standard_Err_Delta) = names(param.effects)
  
  save(Standard_Err_Delta, file = paste0(path,"Standard_Err_Delta.RData"))
  save(Dispersions, file = paste0(path,"Dispersions.RData")) 
  save(True_control, file = paste0(path,"True_Control.RData")) 
  save(True_effect, file = paste0(path,"True_Effect.RData")) 
  save(True_treatment, file = paste0(path,"True_Treatment.RData")) 
  save(Sim_metadata, file = paste0(path,"Simulated_Metadata.RData")) 
  save(Sim_data, file = paste0(path,"Simulated_Data.RData")) 
}




# data_sim<-function(param.disp,param.cont,n_sim,n_samples,n_otu,filter=10){
#   #Simulate data
#     true_control <- rnorm(n=n_otu,mean = param.cont$mean, 
#                                      sd =param.cont$sd)
#     
#     true_effect <- rt(n=n_otu, df, ncp) 
#     true_treatment <-   true_effect + true_control
#     mean_abund <- c(2^true_treatment, 2^true_control)  
#     
#     dispers <- function(mean_abund,c0,c1){c0 + c1/mean_abund}
#     disp <- drop(1/apply(mean_abund, 2, dispers,  c0 = param.disp[["asymptDisp"]],
#                          c1 = param.disp[["extraPois"]]))
#     df <- data.frame(mean_abund,disp)
#     
#   Sim_data <- list()
#   for(i in 1:n_sim){
# 
#     new_otu <- 2*n_otu
#     
#     #Simulate dispersions
#     data <-matrix(apply(df,1,function(x,n=n_samples){rnbinom(n=n,mu=x[1],size=x[2])}),
#                   nrow = n_otu,ncol = n_samples)
#     
#     #keep <- rowSums(data) >= filter
#     #data <- data[keep,] 
#     #select <- sample.int(n_otu)
#     #data <- data[select,]
# 
#     rownames(data) <- paste0("otu",1:n_otu)
#     colnames(data) <- paste0("Sample",1:n_samples)
#     True_data[[i]] <- list(true_effect,true_control,true_treatment)
#     Sim_data[[i]] <- data
#   } 
#   
#   names(True_data) = paste0("True_data",1:n_sim)
#   names(Sim_data) = paste0("Sim_data",1:n_sim)
#   save(Sim_data, file = paste0(path, "Simualated_Data.RData")) 
# }

################################################################################
met <- function(metadata){
  metadata_list <- list(); rownames(metadata) <- metadata$Samples
  for(i in 1:n_sim){
    metadata_list[[i]] <- metadata
  }
  names(metadata_list) <- names(data_list)#paste0("Sim_data",1:n_sim)
  return(metadata_list)
}
###############################################################################
Deseq1 <-function(Data_List,Metadata_List,path){
  
  result = dispersion_coef = RESULT= list();ind= 1
  for(i in names(Data_List)){
    data_list <- Data_List[[i]]
    metadata_list = Metadata_List[[i]]
    
    index = 1
    for(n in names(data_list)){
      data = data_list[[n]]; metadata <- metadata_list[[n]]
      
      #View(data); View(metadata)
      
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
    names(result) <- paste0("result",1:length(data_list))
    RESULT[[ind]] <- result
    ind = ind + 1
  }
  names(RESULT) <- names(Data_List)
  save(RESULT, file = paste0(path, "Sim_Deseq_RESULT.RData")) 
}
#######################################################################################
#Fitting Effect sizes
Fit_Effects <-function(effect_list,logcontrol_list,plot_path,plotname){
  
  Parameters = plts = list(); index <- 1
  for(n in names(effect_list)){
    logcont <- logcontrol_list[[n]]; effect <- effect_list[[n]]
    k <- which(is.finite(logcont)) 
    logcont <- logcont[k]; effect <- effect[k]
    dat <- data.frame(logcont,effect)
    
    mm <- mle2(effect ~ TruncCauchy(loc=0,scale=exp(a+b*logcont + c*logcont^2),
                               low =-5, high = 5),
               start =list(a=0.1, b=0.1, c= 0.01), data=dat)
    
    # mm<-mle2(effect ~ dcauchy(0, scale=exp(a+b*logcont + c*logcont^2)),
    #          start =list(a=0.1, b=0.1, c= 0.001), data=dat)
    
    Parameters[[index]] <- data.frame(a = coef(mm)[1], b = coef(mm)[2],c = coef(mm)[3])

    # scal =  exp(coef(mm)[1]+coef(mm)[2]*logcont + coef(mm)[3]*logcont^2)
    # sim_effect <- sapply(scal,
    #             \(s) rtrunc(
    #               n = 1, spec = 'cauchy', a = -5, b = 5,
    #               location = 0, scale = s))
      # 
      # effect <-apply(scal,1,function(x){
      #   rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})
    #   
    # df = data.frame(Effect = unlist(list(effect,sim_effect)),
    #                 group = rep(c("observations","simulations"),
    #                             each=nrow(dat))
    # )
    # df$repeated_observations=rep(df[df$group == "observations",]$Effect, 2)

    file_name = paste0(plot_path, plotname, n , ".pdf", sep="")
    
    pdf(file_name)
    
    # pp <- ggplot(df, aes(x = Effect,color=group)) + 
    #   geom_histogram(aes(x=repeated_observations,
    #                      y = (after_stat(density))), 
    #                  colour = 1, fill = "white") +
    #   ylab("density") +   geom_density() +   ggtitle(n) + 
    #   scale_color_manual(values=c("red", "blue"))
    pp <- ggplot(dat, aes(x = effect)) +
      geom_histogram(aes(y = (after_stat(density))),colour = 1, fill = "white") +
      ylab("Density") +   geom_density() +  ggtitle(n) +
      scale_color_manual(values=c("blue")) + labs(x = "LogFoldChange")
    dev.off()
    pp
    plts[[index]] <- pp
    index = index + 1
  }
  names(Parameters)= names(effect_list)
  save(Parameters, file = paste0(path, "Effect_Parameters.RData")) 
  plts
}

#######################################################################################
##' Simulating control abundance 
##' @Control_Effect_Simulate simulates control and effect sizes
##' @control_param list of parameters for estimating the sd of norm distribution
##'      for control abundance 
##' @effect_parameters list of parameters for 
##' @nsim number of simulations
Control_Effect_Simulate <- function(control_param, effect_param, nsim,plot_path, plotname){
  
  for(n in names(control_param)){
    
    DATA_control <- rnorm(nsim,mean = control_param[[n]]$mean, sd=control_param[[n]]$sd)
    
    DATA_effect <- rnorm(nsim,mean = 0, sd= effect_param[[n]]$intercept 
                         + mean(effect_param[[n]]$coefficient*DATA_control))
    
    DATA <- data.frame(Control_Abundance = DATA_control, Effect_size = DATA_effect)
    file_name = paste0(plot_path, plotname, n , ".pdf", sep="")
    pdf(file_name, width=6, height=5)
    plt <- ggplot(DATA, aes(Control_Abundance,Effect_size)) +
      geom_point()  + 
      geom_smooth() +
      #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
      ggtitle(paste0("Simulated Control and Simulated Fold Changes")) +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16))
    print(plt)
    dev.off()
  }
  
}

#########################
# difference_plot <- function(Transform_data, plot_path,plotname){
#   
#   for(n in names(Transform_data)){
#     dd <- Transform_data[[n]]
#     data <- data.frame(Fold_change = (dd$control - dd$treatment), 
#                        Control_abundance = (dd$control) ) 
#     
#     file_name = paste0(plot_path, plotname, n , ".pdf", sep="")
#     pdf(file_name)
#     plt <- ggplot(data, aes(Control_abundance, Fold_change)) +
#       geom_point()  + 
#       geom_smooth() +
#       #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
#       ggtitle(paste0("Control against Fold Changes")) +
#       theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
#     print(plt)
#     dev.off()
#     print(plt)
#     
#   }
# } 

###################################################################
###################################################################

Impute_plot <- function(Imputed_data, plot_path,plotname){
  
  effects_list <- list(); index <- 1
  
  for(n in names(Imputed_data)){
    dd <- Imputed_data[[n]]
    
    Fold_change = log2(dd$control/dd$treatment); Control_abundance = log2(dd$control) 
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
###################################################################
###################################################################

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
###################################################################



##' @HIST_EFFECT_SIMULATED plots histogram for simulated effect sizes
##' @deseq_list            deseq results containing log foldchanges
##' @logcontrol            list of log control abundances
##' @file_path             path to save files
##' @plotname              file names

HIST_EFFECT_SIMULATED <- function(deseq_list,logcontrol,file_path,plotname){
  
  Parameters <- list(); index <- 1
  
  for(n in names(logcontrol)){
    
    logcont<-logcontrol[[n]]; effect_size <-deseq_list[[n]]$log2FoldChange  
    a<- which(is.finite(logcont))
    logcont <- logcont[a]; effect_size2 <- effect_size[a]
    
    file_name = paste0(plot_path, plotname, n , ".pdf", sep="")
    pdf(file_name, width=6, height=4)
    par(mfrow= c(1,2))
    hist(effect_size2, breaks=30, prob = TRUE,
         main='Effect size',
         cex.main=1)
    lines(density(effect_size2), col="blue", lwd=2) # add a density estimate with defaults
    data <- data.frame(logcont,effect_size2)
    #my_sd <-glmmTMB(effect_size ~ Normal(0, exp(a+b*x))
    
    #my_sd <-  mle2(effect_size ~ dnorm(0, s),parameters = list(s ~ controlabund))
    
    intercept<-coef(my_sd)[1]; coefficient <- coef(my_sd)[2]
    
    dat <- data.frame(intercept, coefficient)
    row.names(dat) <- NULL
    Parameters[[index]] <- dat 
    
    p<-(intercept + coefficient*logcont)
    l<- length(p)
    hist(rnorm(l,mean = 0, sd=mean(p)), breaks=30, prob = TRUE,
         main='Simulated effect size',
         cex.main=1)
    lines(density(effect_size2), col="blue", lwd=2) # add a density estimate with defaults
    dev.off()
    index <- index + 1
    
  }
  
  names(Parameters) <- names(logcontrol)
  save(Parameters, file = paste0(file_path, "EffectSize_Parameters.RData")) 
}

##############################################################################
##############################################################################
##' @Hist_QQ plots histogram and QQ plot for fitted normal distribution 
##'          and saves the parameter estimates
##' @logcontrol list containing control abundances for each dataset     
##' @plot_path path to save plots
##' @file_path path to save parameter estimates
##' @plotname  plot names
 
##' HIST <- function(effects_list,param.effects,param.controls,plot_path, plotname)

HIST <- function(effects_list,plot_path, plotname){

  index = 1; plts = Prop_of_Zeroes=list()
  for(n in names(effects_list)){

    data <- effects_list[[n]]
    Eff <- as.numeric(data$effect)
    effects <- Eff[which(is.finite(Eff))]
    Prop_of_Zeroes[[index]] <- length(effects)/length(Eff)
    df <- data.frame(effects)
    
    ###############
    # #simulate control and effects 
    # control <- rnorm(n=length(effects), mean = mu, sd = sd)
    # #cont <- unlist(control)
    # 
    # ## scale parameters values for cauchy 
    # scal = data.frame(x =  exp(param.effect[,1]+param.effect[,2]*control + 
    #                              param.effect[,3]*control^2))
    # 
    # effect <-apply(scal,1,function(x){
    #   rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})
    # #################
    
    file_name = paste0(plot_path, plotname, n , ".pdf", sep="")
    pdf(file_name)
    pp <- ggplot(df, aes(x = effects)) + 
      geom_histogram(aes(y = (after_stat(density))),colour = 1, fill = "white") +
      ylab("density") +   geom_density()
    print(pp)
    dev.off()
    
    plts[[index]] <- pp
    index = index + 1
  }
  list(hist = plts,Prop_of_Zeroes)
}
  
###########################################################################
Fit_control <- function(logcontrol,plot_path,file_path, plotname){
  
  Parameters = plts = list();   index <- 1
  for(n in names(logcontrol)){
    
    dat <- logcontrol[[n]]; dat <- data.frame(x=dat[is.finite(dat)])
    
    skewnormal <- mle2(x~dsn(xi=mean, omega=scale_param, alpha=slant_param),
                       data = dat,
                       start = list(mean = 0, scale_param = 2, slant_param = 2))
  
    #trucnormal <- mle2(x ~ dtnorm(mean, sd = exp(logsd), a = lower),
                #data = dat,
                #method = "Nelder-Mead",
                #start = list(mean = 0, logsd = 0,lower = min(dat$x)-0.1))

    est = coef(skewnormal); Parameters[[index]] <- as.data.frame(t(est))
    sim_cont = rsn(n=nrow(dat), xi = est[[1]], omega = est[[2]], alpha = est[3]) 
    #sim_cont =rtnorm(nrow(dat), mean = est[[1]], sd = exp(est[[2]]), a = est[3]) 
    
    df = data.frame(control = unlist(list(dat$x,sim_cont)),
                    group = rep(c("observations","simulations"),
                    each=nrow(dat))
                    )
    df$repeated_observations=rep(df[df$group == "observations",]$control, 2)

    #############################################
    file_name = paste0(plot_path, plotname, n , ".png", sep="")
    
    png(file_name)
    pp <- ggplot(df, aes(x = control,color=group)) + 
      geom_histogram(aes(x=repeated_observations,
                         y = (after_stat(density))), 
                     colour = 1, fill = "white") +
                    ylab("density") +   geom_density() +   ggtitle(n) + 
                    scale_color_manual(values=c("red", "blue"))

    plts[[index]] = pp    
    dev.off()
    
    index <- index + 1
  }
  names(Parameters) <- names(logcontrol)
  save(Parameters, file = paste0(path, "Control_Parameters.RData"))
  names(plts) = names(logcontrol)
  plts
}
##############################################################################
##############################################################################
##' @DIMS computes the dimensions of dataframe in a given list
##' @data_list list containing dataframes

DIMS <- function(data_list){
  dims <- list(); count <- 1
  
  for(n in names(data_list)){
    
    data <- data_list[[n]]
    ##Check data is samples  by taxa
    m<-  rownames(data)[1]; p <- colnames(data)[1]
    if (nchar(m) < nchar(p) ) {
      data <- t(data)
    } else {
      data <- (data)
    }
    
    dims[count] <- list(dim(data))
    count <- count + 1
  }
  names(dims) <- names(data_list)
  return(dims)
}


#############################################################################
##############################################################################
##' @Plot_ContEffects function to plot control verse effect sizes 
##' @deseq_list list of results from deseq2
##' @logcontrol list of log control data
# Plot_ContEffects_multi_filters <- function(deseq_list,logcontrol, plot_path){
#   
#   for(n in names(deseq_list)){
#     
#     res<-deseq_list[[n]]; logcont <- logcontrol[[n]] 
#     foldchang <- res$log2FoldChange
#     
#     finite<-which(is.finite(logcont))
#     control <- logcont[finite]; foldchange <- foldchang[finite]
#     data <- data.frame(control, foldchange)
#     #####################################################################
#   }
#   
#   df = as.data.frame(matrix(rnorm(16),4,4))
#   colnames(df) = c("a","b","c","d")
#   plts <- list()
#   
#   df = as.data.frame(matrix(rnorm(16),4,4))
#   colnames(df) = c("a","b","c","d")
#   df1 = as.data.frame(matrix(rnorm(9),3,3))
#   colnames(df1) = c("a1","b1","c1")
#   df2 = as.data.frame(matrix(rnorm(4),2,2))
#   colnames(df2) = c("a2","b2")
#   
#   
#   
#   ggplot2::ggplot(df,aes(x = a, y=b)) +
#     geom_point(aes(x = a, y=b),color="red") +
#     geom_point(data = df1,aes(x = a1, y=b1),color="blue") +
#     geom_point(data = df2,aes(x = a2, y=b2),color="green") 
#   
#   plts[[index]] <- plt
# }
##############################################################################
##############################################################################
##' @Plot_ContEffects function to plot control verse effect sizes 
##' @deseq_list list of results from deseq2
##' @logcontrol list of log control data

Plot_ContEffects <-function(deseq_list,logcontrol, plot_path){
  
  stopifnot(names(deseq_list) == names(logcontrol))
  index <- 1; plt1 =plt2=list()
  
  for(n in names(deseq_list)){  
    
    res<-deseq_list[[n]]; logcont <- logcontrol[[n]] 
    foldchang <- res$log2FoldChange
    
    finite<-which(is.finite(logcont))
    control <- logcont[finite]; foldchange <- foldchang[finite]
    #####################################################################
    data <- data.frame(control, foldchange)
    
    file_name = paste0(plot_path, "control_effects_", n , ".png", sep="")
    pdf(file_name)
    plts1 <-  ggplot(data, aes(control, foldchange)) +
      geom_point(alpha = 0.1)  + 
      geom_smooth() +
      #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
      ggtitle(n) +
      xlab("control (log2 scale)") + ylab("log2 fold changes") +
      theme(legend.position = "none")
      #theme(plot.title = element_text(hjust = 0.5))#, text = element_text(size = 20))
    print(plts1)
    dev.off()
    
    plt1[[index]] <- plts1
    #####################################################################
    #Scale-location plots
    dat <- data.frame(control, variance=foldchange^2)
    
    file_nam = paste0(plot_path, "scale_location_", n , ".png", sep="")
    pdf(file_nam)
    plts2 <-  ggplot(dat, aes(control, variance)) +
      geom_point(alpha=0.1)  + 
      geom_smooth() +
      #geom_smooth(method = "gam") +
      geom_smooth(method = "lm", formula = y ~poly(x,2), color="red") +
      ggtitle(n) +
      xlab("control (log2 scale)") + ylab("(fold changes)^2") +
     # theme(legend.position = "none")
    theme(plot.title = element_text(hjust = 0.5))#, text = element_text(size = 20))
    print(plts2)
    dev.off()
    
    plt2[[index]] <- plts2
    
    index <- index + 1
  }
  list(plt1, plt2)
}

##############################################################################
##############################################################################
##' @deseq_res function to find the foldchanges, Pvalues, adjpvalues, basemeans,
##'            from deseq2 and write results as a list of dataframes
##' @data_list list of dataset
##' @metadata_list list of metadata 

#outlier_thresh 
deseq_res<- function(data_list, metadata_list,plot_path, plotname,filter=10){
  
  DeseqResults = dispersion_coef= list(); Control_Treat  <- list(); index <- 1; 
  Prop_zeros= Effects_list= list()
  
  for(n in names(data_list)){
    
    data  = data_list[[n]]; metadata = metadata_list[[n]]
    #############################################################################
    ##Check data is taxa  by samples
    m<-  rownames(data)[1]; p <- colnames(data)[1]
    if (nchar(m) < nchar(p) ) {
      data <- data
    } else {
      data <- t(data)
    }
  
    #############################################################################
    ##Deseq processes
    dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
    keep <- rowSums(counts(dds)) >= filter
    dds <- dds[keep,]; data <- data[keep,] 
    #############################################################################
    #' Remove samples with all zero (deseq needs this to run effectively) 
    r <- which(colSums(data) == 0)
    if (length(r) > 0){
      data <- data[,-r]; metadata <- metadata[-r,]
      dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
    }
    #############################################################################
    dds$Groups <- relevel(dds$Groups, ref = "NT")
    dds <- DESeq(dds,sfType ="poscounts", minReplicatesForReplace=Inf) 
    res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
    reslt <- lfcShrink(dds, res=res, coef=2, type="normal") 
    result <- data.frame(reslt)
    DeseqResults <- c(DeseqResults, list(result))
    Effects_list[[index]] <- result$log2FoldChange 
    #############################################################################
    ##Control Treatment splits
    k <- which(metadata$Groups == "ASD")
    control <- data[,-k]; treatment <- data[,k] 
    stopifnot(colnames(control) == metadata$Samples[-k])
    stopifnot(colnames(treatment) == metadata$Samples[k])
    #Prop_zeros[[index]] <-  length(k)/length(cont_Means)
    #############################################################################
    cont_Means <- rowMeans(control); treat_Means <- rowMeans(treatment)
    logcontrol = log2(cont_Means);   logtreat = log2(treat_Means)
    k_cont <-which(cont_Means != 0); dispers = dispersions(dds)
    
    dat <- data.frame(logcontrol =logcontrol[k_cont], dispersion = dispers[k_cont])
    mm <- nls(dispersion ~ a + b/exp(logcontrol), start =list(a=1, b= 1), 
              data=dat)
    
    dispersion_coef[[index]] <- data.frame(asymptDisp = coef(mm)[1], extraPois =  coef(mm)[2])
    
    Control_Treat[[index]] <- list(Control_Data = control, Control_Means = cont_Means,
                                   LogControl  =  logcontrol, Treat_Data =  treatment,
                                   Treat_Means  = treat_Means,
                                   LogTreat  = logtreat)
    ############################################################################
    index <- index + 1
  }
  
  names(dispersion_coef) = names(Effects_list) = names(data_list)
  names(Control_Treat) = names(DeseqResults) = names(data_list)
  
  save(Effects_list, file = paste0(path, "DeseqFoldChanges.RData")) 
  save(dispersion_coef, file = paste0(path, "Dispersion_Parameters.RData")) 
  save(Control_Treat, file = paste0(path, "Control_Treat_Data.RData"))
  save(DeseqResults, file = paste0(path, "DeseqResults.RData"))
  return(Prop_zeros) 
}

########################################################################### 
########################################################################### 
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
########################################################################### 
########################################################################### 
##' @LogMeanData function to compute the log means of control and treatment
##' for all the dataset
##' @cont_treat_list list containing control and treatment data for all dataset

LogMeanData <- function(cont_treat_list,path){
  
  cont_treat = list()
  ControlMeans = list()
  TreatmentMeans = list()
  
  for(n in names(cont_treat_list)){
    cont  = cont_treat_list[[n]][["control"]]
    treat = cont_treat_list[[n]][["treat"]]
    
    #########Sanity check: data must be sample by taxa #####
    m <-  rownames(cont)[1]; n<- colnames(cont)[1]
    if (nchar(m) < nchar(n) ){
      cont <- t(cont)
    } else {
      cont <- cont
    }
    
    m<-  rownames(treat)[1]; n<- colnames(treat)[1]
    if (nchar(m) < nchar(n) ) {
      treat <- t(treat)
    } else {
      treat <- treat
    }
    
    #########Compute means
    contmean <- (log2(colMeans(cont)))
    #stopifnot(is.finite(contmean))
    #colnames(contmean) <- "logControl"
    treatmean <- (log2(colMeans(treat)))
    #colnames(treatmean) <- "logTreatment"
    #stopifnot(is.finite(treatmean))
    
    ControlMeans <- c(ControlMeans, list(contmean))
    TreatmentMeans <- c(TreatmentMeans, list(treatmean))
    
    cont_treat  <- c(cont_treat, list(list(ContMean = contmean, TreatMean = treatmean) ))
  }
  names(cont_treat) <- names(cont_treat_list)
  names(ControlMeans) <- names(cont_treat_list)
  names(TreatmentMeans) <- names(cont_treat_list)
  
  save(cont_treat, file = paste0(path,"Control_N_TreatmentMeans.RData"))
  save(ControlMeans, file = paste0(path,"ControlMeans.RData"))
  save(TreatmentMeans, file = paste0(path,"TreatmentMeans.RData"))
}


########################################################################### 
########################################################################### 
##'@RDataFiles function to create filtered data, control data, treatment data
##' and metadata as .RData files
##' @L list containing data and original meta data
##'thresh was threshold values for filtering
RDataFiles <- function(L, path){
  #stopifnot(names(L) == colnames(thresh))
  data_list  <-  list(); meta_list <- list() 
  cont_list  <- list();  treat_list  <- list()
  index <- 1
  
  for(n in names(L)){
    dat <- L[[n]][["data"]]; meta_dat <-L[[n]][["metadata"]]
    v <- filter_data(dat,meta_dat)
    data_list[[index]] <- as.data.frame(v[["data"]])
    cont_list[[index]] <- as.data.frame(v[["treat_data"]])
    treat_list[[index]] <- as.data.frame(v[["cont_data"]])
    meta_list[[index]]  <- as.data.frame(v[["meta_data"]])
    index = index + 1
  }
  
  names(data_list) <- names(L)
  names(cont_list) <- names(L)
  names(treat_list) <- names(L); names(meta_list) <- names(L)
  
  #save(data_list, file = paste0(path,"data_new.RData"))
  save(data_list, file = paste0(path,"data_filtered.RData"))
  save(cont_list, file = paste0(path,"control.RData"))
  save(treat_list, file = paste0(path,"treat.RData")) 
  save(meta_list, file = paste0(path,"metadata.RData"))
}

##############################################################################
##############################################################################
##'@read_datasets function to read datasets as list 
##'@data_list list of data
##'@metadata_list list of metadata
read_datasets <- function(data_list, metadata_list){
  L <- list()
  for (n in names(data_list)){ 
    L <- c(L, list(list( data = as.data.frame(data_list[[n]]), 
                         metadata = as.data.frame(metadata_list[[n]]))))
  } 
  names(L) <- names(data_list) 
  return(L)
}

##############################################################################
##############################################################################
##' @read_ContTreat function to combine control and treatment data in a list 
##' @cont_list list of control data
##' @treat_list list of treatment data

read_ContTreat <- function(cont_list, treat_list){
  L <- list()
  for (n in names(cont_list)){ 
    L <- c(L, list(list( control = as.data.frame(cont_list[[n]]), 
                         treat = as.data.frame(treat_list[[n]]))))
  } 
  names(L) <- names(cont_list) 
  return(L)
}

##############################################################################
##############################################################################
##' @filter_data function to filter based on counts per million (CPM)
##' @cmp_threshold threshold for counts per million (CPM)
##' @nOTU_thesholds threshold for number of OTUs to apply cpm on 

filter_data <- function(data, meta_data, CMP_threshold = 100, nOTU_theshold = 2){
  
  ##Previous filtering Method used
  #filter_data <- function(data, meta_data,filter){
  ##Check data is given as samples by taxa 
  ##this means that rows < columns since samples is < taxa
  # if (nrow(data) < ncol(data) ) {
  #    data <- data
  # } else {
  #   data <- t(data)
  # }
  
  #stopifnot(rownames(data)==meta_data$Samples)
  #weights <-colSums(data)/sum(data)
  
  #keep <- weights > filter
  # data <- data[,keep, drop = FALSE]
  
  ## Method 2
  ##Check data is given as OTU by samples 
  if (nrow(data) < ncol(data) ) {
    data <- t(data)
  } else {
    data <- data
  }

  stopifnot(colnames(data)==meta_data$Samples)
  
  #' Remove samples with all zero (cpm needs this to run effectively) 
  k <- which(colSums(data) == 0)
  if (length(k) > 0){
    data <- data[,-k]; meta_data <- meta_data[-k,]
  }

  #'Sanity check
  stopifnot(meta_data$Samples == colnames(data))
 
  y <- DGEList(counts=data,group=meta_data$Groups)
  keep <- rowSums(cpm(y)>CMP_threshold) >= nOTU_theshold
  y <- y[keep,]
  data  <- y$counts; meta_data <- y$sample
  
  #print(dim(data))
  
  #' Remove samples with all zero (deseq needs this to run effectively) 
  k <- which(colSums(data) == 0)
  if (length(k) > 0){
    data <- data[,-k]; meta_data <- meta_data[-k,]
  }

  ## Extract control and treatment samples
  colnames(meta_data)[1] <- "Groups"
  asd <- which(meta_data$Groups== "ASD")
  treat_data <- (t(data))[asd,]
  cont_data <- (t(data))[-asd,]
  
  stopifnot(rownames(meta_data)[asd] == rownames(treat_data))
  stopifnot(rownames(meta_data)[-asd] == rownames(cont_data))
  
  k <- which(colSums(treat_data) == 0) #; print(length(k))
  r <- which(colSums(cont_data) == 0) #; print(length(r))
  
  stopifnot((meta_data$Samples)[asd] == rownames(treat_data))
  stopifnot((meta_data$Samples)[-asd] == rownames(cont_data))
  
  v <-list(t(data), treat_data, cont_data, meta_data)
  names(v) <- c("data","treat_data", "cont_data", "meta_data")
  
  return(v)
}

##############################################################################
##############################################################################
##' @param rdatafile the name of a data file
##' @param new a vector of names of new objects to add  
##' @example add_data("mydata.RData", new = c("extra1", "extra2"))
add_data <- function(rdatafile,new) {
  n <- load(rdatafile)
  save(list = c(n, new), file = rdatafile) 
}


##############################################################################
#MontiCarlo sim function for computing the power for various samples and read depths
Powers <-function(treat_data,cont_data,numMC,sample_size,Reads){
  
  #######Dichlet parameters
  treat_param <- DM.MoM(treat_data)
  cont_param  <- DM.MoM(cont_data)
  pi0 <- cont_param$pi
  ####################################
  pval <- matrix(sample_size, length(sample_size),length(Reads) + 1) 
  ################Loop#################
  for (j in 1:length(sample_size)){
    
    for (i in 1:length(Reads)) {
      nrsGrp1 <- rep(Reads[i], sample_size[j])
      nrsGrp2 <- rep(Reads[i], sample_size[j])   
      groups <- list(nrsGrp1, nrsGrp2)
      
      group_theta <- c(treat_param$theta, cont_param$theta)
      group_p <- rbind(treat_param$pi, cont_param$pi)
      pval2 <- MC.Xmc.statistics(groups, numMC, pi0, group.pi=group_p,   group.theta=group_theta)
      pval[j,i + 1] <-pval2
    }
  } 
  
  return(pval)
}

##############################################################################
################Compute overall effect size#######################
effect_all <- function(treat_data,cont_data){
  
  
  ##1. Ensure data is given as samples by taxa
  ## My taxa names are shorter than my sample names so I can check as follows
  row <- rownames(full_data)[1]
  col <- colnames(full_data)[1]
  
  if (nchar(row) < nchar(col) ) {
    full_data <- t(full_data)
  } else {
    full_data <- full_data
  }
  # treat_param <- DM.MoM(treat_data)
  # cont_param  <- DM.MoM(cont_data)
  group_data <- list(treat_data, cont_data)
  effect <- Xmcupo.effectsize(group_data)
  return(effect)
}

###Order by row sum
Order_taxa <- function(treat_data,cont_data){
  treat_ord <- treat_data[order(rowSums(treat_data),decreasing=T),]
  cont_ord <- cont_data[order(rowSums(cont_data),decreasing=T),]
  return(list(treat_ord,cont_ord))
}

##############################################################################
##############################################################################
Mean_Lists <- function(Data_list){
  mean<-vector(mode='list')
  for (i in 1:length(Data_list)) {
    mean[[i]]    <- colSums(Data_list[[i]])/sum(Data_list[[i]])
  }
  return(mean)
}

##############################################################################
##############################################################################
Data_simulate <-function(path,treat_means,cont_means,dispersion,n_sim,sample_per_Grp){ 
  
  n_taxa <- length(treat_means)
  sim_treat <- matrix(0,n_taxa,sample_per_Grp)
  sim_cont <- matrix(0,n_taxa,sample_per_Grp)
  
  for (j in 1: n_sim){
    
    for (i in 1:n_taxa) {
      sim_treat[i,] <- rnegbin(n=sample_per_Grp,mu = treat_means[i], 
                               theta = dispersion[i])
      sim_cont[i,] <- rnegbin(n=sample_per_Grp,mu = cont_means[i],
                              theta = dispersion[i])
    }
    
    
    sim_data <-data.frame(sim_treat,sim_cont)
    rownames(sim_data) <- paste0("sp",1:nrow(sim_data))
    colnames(sim_data) <- paste0("sample",1:ncol(sim_data))
   write.csv(sim_data, file= file.path(path, paste0("sim",j,".csv",sep="")))
    
  }
}
##############################################################################
##############################################################################
results_data<-function(path, my_data,conditions){
  
  for (i in 1:length(my_data)){
    dds <- DESeqDataSetFromMatrix(my_data[[i]], conditions, ~Groups)
    ddsDE <- DESeq(dds, sfType ="poscounts") 
    res <- results(ddsDE, alpha = 0.05, 
                   contrast = c("Groups", "ASD","Healthy"))
    write.csv(res, file= file.path(path, paste0("result",i,".csv",sep="")))
  }
}

##############################################################################
##############################################################################
# deseq_res<- function(data, meta_data){
#   ##Check data is taxa  by samples 
#   m<-  rownames(data)[1]; n<- colnames(data)[1]
#   if (nchar(m) < nchar(n) ) {
#     data <- data
#   } else {
#     data <- t(data)
#   }
#   
#   #res_ <- deseq_res(data_new, meta_dat_)
#   
#   dds <- DESeqDataSetFromMatrix(data, meta_data, ~Groups)
#   ddsDE <- DESeq(dds, sfType ="poscounts") 
#   res <- results(ddsDE, alpha = 0.05, 
#                  contrast = c("Groups", "ASD","NT"))
#   return(res)
# }

##############################################################################
##############################################################################
##' @save_data_filt stores the filtered data into rdatafiles
##' @v list containing filtered data, treatment data, control data and metadata
##' @name name of the dataset
##' @path list of path names to rdatafiles
save_data_filt <- function(v,name, path)
{
  n<-length(v)
  stopifnot(names(v)==c("data","treat_data","cont_data","meta_data"))
  stopifnot(n==length(path))
  
  names(v) <- rep(name,n)
  for(i in 1:n){
    name <- v[[i]]
    add_data(path[i],name)
  }
}
##############################################################################
##############################################################################