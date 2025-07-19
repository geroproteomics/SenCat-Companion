
#Function that performs cross validation of groups of features using 10 iterations of 90% train, 10% test to compare median CV R-Squared
en_cv_test_fast <- function(clin_df, protein_list, cat_control, cont_control, trait_list, alpha, num_ensp, composite=FALSE, direction="up",study="B",translate=TRUE,scale_iv=TRUE){
  
  #Setup backend to use many processors
  totalCores = detectCores()
  #Leave one core to avoid overload your computer
  cluster <- makeCluster(totalCores[1]-1)
  registerDoParallel(cluster)
  
  ##make master list of top ensps by tissue and trait
  df_master_ensp <- vector(mode = "list", length = 0)
  
  ##make df to hold all Rsquared for 10 trials for all traits and tissues
  df_master_Rsquared <- as.data.frame(matrix(NA,nrow=0,ncol=3))
  colnames(df_master_Rsquared) <- c("Rsquared","Tissue","Trait")
  
  ##loop through sets of proteins, find ensp, cv predictive potential 
  for(val in names(protein_list)){
    ##make df to hold ENSPs for each trait
    df_ensps_per_trait <- as.data.frame(matrix(NA,nrow=num_ensp,ncol=length(trait_list)))
    colnames(df_ensps_per_trait) <- trait_list
    
    ##determine ENSPs for all traits
    #check if translation is needed
    if(translate == TRUE){
      #convert protein list to BLSA/InCHIANTI nomenclature
      if(study == "B"){
        temp_proteins <- df_sencat_blsa[df_sencat_blsa$gene %in% protein_list[[val]], "unique_gene"]
      }
      if(study == "I"){
        temp_proteins <- df_sencat_inchianti[df_sencat_inchianti$gene %in% protein_list[[val]], "unique_gene"]
      }
    }
    if(translate == FALSE){
      temp_proteins <- protein_list[[val]]
    }
    temp_en_object <- en_repeat_fast(clin_df, temp_proteins, cat_control, cont_control, trait_list, alpha, heatmap=FALSE, scale_iv=scale_iv)
    #pull coef df
    temp_en_coef_subset <- temp_en_object$coef
    
    #loop through all clinical traits; find cv Rsquared
    for(i in trait_list){
      #subset to only rows with trait listed
      clin_df_i <- clin_df[!is.na(clin_df[,i]),]
      #scale trait
      clin_df_i[,i] <- as.numeric(scale(clin_df_i[,i]))
      #if composite is false, select all ENSP and use for predictive potential
      if(composite == FALSE){
        #find ensps, subset to specific number 
        temp_ensps <- temp_en_coef_subset[i,temp_en_coef_subset[i,] != 0]
        temp_ensps <- temp_ensps[!colnames(temp_ensps) %in% c(cat_control,cont_control)]
        #proceed if 1+ ENSPs identified
        if(length(temp_ensps) >0){
          temp_ensps <- t(temp_ensps)
          temp_ensps <- as.data.frame(temp_ensps)
          #add dummy column for re-ordering 
          temp_ensps$dumb <- NA
          temp_ensps <- temp_ensps[rev(order(abs(temp_ensps[,i]))),]
          if(nrow(temp_ensps) <=num_ensp){
            temp_ensps <- rownames(temp_ensps)
          }
          else{
            temp_ensps <- rownames(temp_ensps)[1:num_ensp]
          }
          #fill df of specified number of ENSPs per trait
          df_ensps_per_trait[1:length(temp_ensps),i] <- temp_ensps
          ##scale ENSPS
          if(scale_iv ==TRUE){
            clin_df_i[,temp_ensps] <- scale(clin_df_i[,temp_ensps])
          }
          #make IV list using ensps
          IV <- paste(" ~ `",paste(c(temp_ensps),collapse="` + `"),"`",sep="")
          #create forumula
          form <- as.formula(paste("`",i,"`",IV, sep=""))
          ##carot CV
          #train the control
          train_control <- trainControl(method="cv", number=10, seeds=seeds_reproducible)
          #train model
          model <- train(form,
                         data = clin_df_i,
                         method = "lm",
                         trControl = train_control,
                         metric="Rsquared")
          #find 10 iterations of Rsquared
          cv_Rsquared <- model$resample$Rsquared
          #fill master DF with Rsquared, tissue, trait
          temp_df <- cbind(cv_Rsquared,rep(val,10),rep(i,10))
          colnames(temp_df) <- c("Rsquared","Tissue","Trait")
          df_master_Rsquared <- rbind(df_master_Rsquared,temp_df)
        }
      }
      if(composite == TRUE){
        #find ensps, subset to specific number 
        if(direction == "up"){
          temp_ensps <- temp_en_coef_subset[i,temp_en_coef_subset[i,] > 0]
        }
        if(direction == "down"){
          temp_ensps <- temp_en_coef_subset[i,temp_en_coef_subset[i,] < 0]
        }
        temp_ensps <- temp_ensps[!colnames(temp_ensps) %in% c(cat_control,cont_control)]
        #proceed if 1+ ENSPs identified
        if(length(temp_ensps) >1){
          temp_ensps <- t(temp_ensps)
          temp_ensps <- as.data.frame(temp_ensps)
          #add dummy column for re-ordering 
          temp_ensps$dumb <- NA
          temp_ensps <- temp_ensps[rev(order(abs(temp_ensps[,i]))),]
          if(nrow(temp_ensps) <=num_ensp){
            temp_ensps <- rownames(temp_ensps)
          }
          else{
            temp_ensps <- rownames(temp_ensps)[1:num_ensp]
          }
          #fill df of specified number of ENSPs per trait
          df_ensps_per_trait[1:length(temp_ensps),i] <- temp_ensps
          ##make composite score
          clin_df_i$composite_score <- scale(rowMeans(clin_df_i[,temp_ensps]))
          #make IV list using ensps
          IV <- paste(" ~ composite_score")
          #create forumula
          form <- as.formula(paste("`",i,"`",IV, sep=""))
          ##carot CV
          #train the control
          train_control <- trainControl(method="cv", number=10, seeds=seeds_reproducible)
          #train model
          model <- train(form,
                         data = clin_df_i,
                         method = "lm",
                         trControl = train_control,
                         metric="Rsquared")
          #find 10 iterations of Rsquared
          cv_Rsquared <- model$resample$Rsquared
          #fill master DF with Rsquared, tissue, trait
          temp_df <- cbind(cv_Rsquared,rep(val,10),rep(i,10))
          colnames(temp_df) <- c("Rsquared","Tissue","Trait")
          df_master_Rsquared <- rbind(df_master_Rsquared,temp_df)
        }
      }
    }
    df_master_ensp[[val]] <- df_ensps_per_trait
  }
  
  #Stop cluster
  stopCluster(cluster)
  
  #return list of accuracy and ensp by trait
  list <- list("Rsquared" = df_master_Rsquared, "ensp_table" = df_master_ensp)
  return(list)
}