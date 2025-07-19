

#function for feature selection using elastic net modeling for continuous parameters 
en_repeat_fast <- function(clin_df, protein_list, cat_control, cont_control, trait_list, alpha, heatmap=FALSE,scale_iv=TRUE){
  #default is no heatmap
  p <- NA
  
  #df to hold optimized alpha and lambda for each clinical trait
  alpha_lambda <- as.data.frame(matrix(NA, ncol=2, nrow=length(trait_list)))
  colnames(alpha_lambda) <- c("alpha","lambda")
  rownames(alpha_lambda) <- trait_list
  alpha_lambda$alpha = alpha
  
  if(!NA %in% cont_control & !NA %in% cat_control){
    #Build matrix to hold coefficients using optimized model with best alpha and lambda for each clinical trait
    coef_glmnet <- as.data.frame(matrix(NA,nrow = length(trait_list), ncol = length(c(cat_control,cont_control,protein_list))))
    rownames(coef_glmnet) <- trait_list
    colnames(coef_glmnet) <- c(cat_control,cont_control,protein_list)
  }
  if(NA %in% cont_control & !NA %in% cat_control){
    #Build matrix to hold coefficients using optimized model with best alpha and lambda for each clinical trait
    coef_glmnet <- as.data.frame(matrix(NA,nrow = length(trait_list), ncol = length(c(cat_control,protein_list))))
    rownames(coef_glmnet) <- trait_list
    colnames(coef_glmnet) <- c(cat_control,protein_list)
  }
  if(!NA %in% cont_control & NA %in% cat_control){
    #Build matrix to hold coefficients using optimized model with best alpha and lambda for each clinical trait
    coef_glmnet <- as.data.frame(matrix(NA,nrow = length(trait_list), ncol = length(c(cont_control,protein_list))))
    rownames(coef_glmnet) <- trait_list
    colnames(coef_glmnet) <- c(cont_control,protein_list)
  }
  if(NA %in% cont_control & NA %in% cat_control){
    #Build matrix to hold coefficients using optimized model with best alpha and lambda for each clinical trait
    coef_glmnet <- as.data.frame(matrix(NA,nrow = length(trait_list), ncol = length(c(protein_list))))
    rownames(coef_glmnet) <- trait_list
    colnames(coef_glmnet) <- c(protein_list)
  }
  
  #for each clinical trait: run 1 iteration of glmnet to find optimized alpha, lambda: calculate BAge, for each trait, fill tables
  for (val in 1:length(trait_list)){
    ##make coef_glmnet with train data
    #remove all rows missing the clinical trait
    values_subset <- clin_df[complete.cases(clin_df[,trait_list[val]]),]
    #isolate predictor proteins and controls only
    if(NA %in% cont_control & NA %in% cat_control){
      proteins_subset <- as.data.frame(values_subset[,c(protein_list)])
    }
    if(!NA %in% cont_control & NA %in% cat_control){
      proteins_subset <- as.data.frame(values_subset[,c(cont_control,protein_list)])
    }
    if(NA %in% cont_control & !NA %in% cat_control){
      proteins_subset <- as.data.frame(values_subset[,c(cat_control,protein_list)])
    }
    if(!NA %in% cont_control & !NA %in% cat_control){
      proteins_subset <- as.data.frame(values_subset[,c(cont_control,cat_control,protein_list)])
    }
    if(scale_iv == TRUE){
      #scale proteins and continuous controls
      if(!NA %in% cont_control){
        for(val2 in c(cont_control,protein_list)){
          proteins_subset[,val2] <- scale(proteins_subset[,val2])
        }
      }
      if(NA %in% cont_control){
        for(val2 in protein_list){
          proteins_subset[,val2] <- scale(proteins_subset[,val2])
        }
      }
    }
    if(scale_iv == FALSE){
      #scale proteins and continuous controls
      if(!NA %in% cont_control){
        for(val2 in c(cont_control)){
          proteins_subset[,val2] <- scale(proteins_subset[,val2])
        }
      }
    }
    #subset to one clinical trait
    clin_subset <- as.data.frame(scale(values_subset[,trait_list[val]]))
    colnames(clin_subset) <- trait_list[val]
    #set seed for reproducibility
    set.seed(1)
    ##set foldid for reproducibility
    foldid_cont <- sample(1:10, size = nrow(clin_subset), replace = TRUE)
    #find selection of proteins using elastic net
    cvfit_subset <- cv.glmnet(as.matrix(proteins_subset),unlist(clin_subset), parallel=TRUE, standardize=FALSE, alpha=alpha,folid=foldid_cont)
    #using optimized lambda, find coefficients and fill df
    if(NA %in% cat_control & NA %in% cont_control){
      coef_glmnet[trait_list[val],] <- as.numeric(coef(cvfit_subset, s = "lambda.min")[c(1:length(c(protein_list))+1)])
    }
    if(!NA %in% cat_control & NA %in% cont_control){
      coef_glmnet[trait_list[val],] <- as.numeric(coef(cvfit_subset, s = "lambda.min")[c(1:length(c(cat_control,protein_list))+1)])
    }
    if(NA %in% cat_control & !NA %in% cont_control){
      coef_glmnet[trait_list[val],] <- as.numeric(coef(cvfit_subset, s = "lambda.min")[c(1:length(c(cont_control,protein_list))+1)])
    }
    if(!NA %in% cat_control & !NA %in% cont_control){
      coef_glmnet[trait_list[val],] <- as.numeric(coef(cvfit_subset, s = "lambda.min")[c(1:length(c(cat_control,cont_control,protein_list))+1)])
    }
  }
  coef_glmnet <- as.data.frame(coef_glmnet)
  
  if(heatmap == TRUE){
    #heatmap plot for all glmnet coefficients
    p <- heatmaply(coef_glmnet,
                   dendrogram = "column",
                   k_col = 4,
                   distfun = "spearman",
                   xlab = "", ylab = "",
                   main = "",
                   margins = c(60,100,40,20),
                   grid_color = "white",
                   grid_width = 0.00001,
                   titleX = FALSE,
                   hide_colorbar = FALSE,
                   branches_lwd = 0.1,
                   label_names = c("Clinical Trait", "Protein", "BTrait"),
                   fontsize_row = 16, fontsize_col = 12,
                   labCol = colnames(coef_glmnet),
                   labRow = rownames(coef_glmnet),
                   title = "Protein Beta by Clincal Trait",
                   heatmap_layers = theme(axis.line=element_blank())
    )
  }
  
  #count IV proteins per trait
  #table to hold IV sums per clinical trait
  IV_sums <- as.data.frame(matrix(NA,nrow=length(trait_list),ncol=1))
  rownames(IV_sums) <- trait_list
  colnames(IV_sums) <- c("IV_Count")
  for(i in trait_list){
    IV_sums[i,] <- length(coef_glmnet[i,coef_glmnet[i,]!=0])
  }
  #make bar graph with IV_sums
  #add clin to column for x label
  IV_sums$clin <- trait_list
  plot_ivsum <- ggplot(data=IV_sums, aes(x = fct_inorder(clin), y = IV_Count)) +
    theme_classic() +
    geom_bar(stat = "identity", position = position_dodge(), color = "blue", fill = "grey")+
    ggtitle("Sum Elastic Net Selected Proteins") +
    ylab("Predictor Variable Sum") +
    xlab("Clinical Trait") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16)) +
    geom_hline(yintercept=length(c(protein_list,control_list)), color = "red", linetype = "dashed", linewidth= 0.75)
  
  #order optimized lambdas
  alpha_lambda <- alpha_lambda[rev(order(alpha_lambda$lambda)),]
  #return beta coefficients, optimized lambda values, heatmap of coef, bar graph of IV sums
  list <- list("coef" = coef_glmnet, "lambda" = alpha_lambda, "heatmap" = p, "ivsum" = plot_ivsum)
  return(list)
}