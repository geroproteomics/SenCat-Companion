

#Function that accepts table with lists of predictors, performs cross-validation using 90% train, 10% test to compare prediction ability by group
en_cv_test_fast_binomial <- function(clin_df, protein_list, cat_control, cont_control, trait_list, alpha, num_ensp, composite=FALSE,study="B",translate=TRUE){
  
  #Setup backend to use many processors
  totalCores = detectCores()
  #Leave one core to avoid overload your computer
  cluster <- makeCluster(totalCores[1]-1)
  registerDoParallel(cluster)
  
  ##make master list of top ensps by tissue and trait
  df_master_ensp <- vector(mode = "list", length = 0)
  
  ##make df to hold all rmse for 10 trials for all traits and tissues
  df_master_auc <- as.data.frame(matrix(NA,nrow=0,ncol=3))
  colnames(df_master_auc) <- c("AUC","Tissue","Trait")
  
  ##loop through sets of proteins, find ensp, cv predictive potential 
  for(val in names(protein_list)){
    #convert protein list to BLSA nomenclature 
    temp_proteins <- protein_list[[val]]
    if(translate == TRUE){
      if(study=="B"){
        temp_proteins <- df_sencat_blsa[df_sencat_blsa$gene %in% temp_proteins, "unique_gene"]
      }
      if(study=="I"){
        temp_proteins <- df_sencat_inchianti[df_sencat_inchianti$gene %in% temp_proteins, "unique_gene"]
      }
    }
    #determine ENSPs per trait
    temp_en_object <- en_binomial_fast(clin_df, temp_proteins, cat_control, cont_control, trait_list, alpha, heatmap=FALSE)
    #pull coef df
    temp_en_coef_subset <- temp_en_object$coef
    
    ##make df to hold ENSPs for each trait
    df_ensps_per_trait <- as.data.frame(matrix(NA,nrow=num_ensp,ncol=length(trait_list)))
    colnames(df_ensps_per_trait) <- trait_list
    
    #loop through all clinical traits; find cv RMSE
    for(i in trait_list){
      #subset to only rows with trait listed
      clin_df_i <- clin_df[!is.na(clin_df[,i]),]
      #make sure that variable is numeric
      #clin_df_i[,i] <- as.numeric(clin_df_i[,i])
      #convert to Yes/No
      clin_df_i[clin_df_i[,i]==0,i] <- "No"
      clin_df_i[clin_df_i[,i]==1,i] <- "Yes"
      #convert trait to factor
      clin_df_i[,i] <- factor(clin_df_i[,i])
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
          clin_df_i[,temp_ensps] <- scale(clin_df_i[,temp_ensps])
          #make IV list using ensps
          IV <- paste(" ~ `",paste(c(temp_ensps),collapse="` + `"),"`",sep="")
          #create forumula
          form <- as.formula(paste("`",i,"`",IV, sep=""))
          ##carot CV
          #train the control
          train_control <- trainControl(method="cv", number=10, classProbs=TRUE,
                                        summaryFunction = twoClassSummary, seeds=seeds_reproducible)
          #train model
          model <- train(form,
                         data = clin_df_i,
                         method = "glm",
                         family = "binomial",
                         trControl = train_control,
                         metric="ROC")
          #find 10 iterations of AUC
          cv_rmse <- model$resample$ROC
          #fill master DF with RMSE, tissue, trait
          temp_df <- cbind(cv_rmse,rep(val,10),rep(i,10))
          colnames(temp_df) <- c("AUC","Tissue","Trait")
          df_master_auc <- rbind(df_master_auc,temp_df)
        }
      }
      if(composite == TRUE){
        #find ensps, subset to specific number 
        temp_ensps <- temp_en_coef_subset[i,temp_en_coef_subset[i,] > 0]
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
          train_control <- trainControl(method="cv", number=10, classProbs=TRUE,
                                        summaryFunction = twoClassSummary, seeds=seeds_reproducible)
          #train model
          model <- train(form,
                         data = clin_df_i,
                         method = "glm",
                         family = "binomial",
                         trControl = train_control,
                         metric="ROC")
          #find 10 iterations of RMSE
          cv_rmse <- model$resample$ROC
          #fill master DF with RMSE, tissue, trait
          temp_df <- cbind(cv_rmse,rep(val,10),rep(i,10))
          colnames(temp_df) <- c("AUC","Tissue","Trait")
          df_master_auc <- rbind(df_master_auc,temp_df)
        }
      }
    }
    df_master_ensp[[val]] <- df_ensps_per_trait
  }
  
  #Stop cluster
  stopCluster(cluster)
  
  #return list of accuracy and ensp by trait
  list <- list("AUC" = df_master_auc, "ensp_table" = df_master_ensp)
  return(list)
}