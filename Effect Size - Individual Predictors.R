

##function that accepts a list of predictor variables, finds individual effect size for list of parameters
effectsize_individual <- function(ensp_list,clin_df,trait_list,cat_control,cont_control,type="lm",study="B",translate=TRUE,scale_iv=TRUE){
  
  ##DF to hold effect sizes with confidence intervals
  df_effectsizes_master <- as.data.frame(matrix(NA,nrow=0,ncol=7))
  colnames(df_effectsizes_master) <- c("tissue","trait","protein","beta","low","high","pval")
  
  #cycle through tissue
  for(i in names(ensp_list)){
    ##cycle through traits
    for(val in trait_list){
      #pull ensp list
      temp_ensps <- ensp_list[[i]][!is.na(ensp_list[[i]])]
      #tralnsate if needed
      if(translate==TRUE){
        #select study
        if(study=="B"){
          temp_ensps <- df_sencat_blsa[df_sencat_blsa$gene %in% temp_ensps,"unique_gene"]
        }
        if(study=="I"){
          temp_ensps <- df_sencat_inchianti[df_sencat_inchianti$gene %in% temp_ensps,"unique_gene"]
        }
      }
      #iterate through genes, find association with trait
      for(val2 in temp_ensps){
        #subset to rows with trait present
        clin_df_i <- clin_df[!is.na(clin_df[,val]),]
        #scale selected ensp
        if(scale_iv==TRUE){
          clin_df_i[,val2] <- as.numeric(scale(clin_df_i[,val2]))
        }
        #scale cont covariates
        if(!NA %in% cont_control){
          for(z_iv_cont in cont_control){
            clin_df_i[,z_iv_cont] <- scale(clin_df_i[,z_iv_cont])
          }
        }
        #factor cat control
        if(!NA %in% cat_control){
          for(z_iv_cat in cat_control){
            clin_df_i[,z_iv_cat] <- factor(clin_df_i[,z_iv_cat])
          }
        }
        #scale DV if lm
        if(type == "lm"){
          clin_df_i[,val] <- as.numeric(scale(clin_df_i[,val]))
        }
        #factor DV if glm
        if(type == "glm"){
          clin_df_i[,val] <- factor(clin_df_i[,val])
        }
        #building formula for modeling
        if(!NA %in% cont_control){
          IV <- paste(" ~ `",paste(c(val2,cat_control,cont_control),collapse="` + `"),"`",sep="")
        }
        if(NA %in% cont_control){
          IV <- paste(" ~ `",paste(c(val2,cat_control),collapse="` + `"),"`",sep="")
        }
        temp_form <- formula(paste("`",val,"`",IV,sep=""))
        #build model
        if(type == "lm"){
          temp_model <- lm(temp_form, data=clin_df_i)
        }
        if(type == "glm"){
          temp_model <- glm(temp_form, data=clin_df_i, family="binomial")
        }
        #add data to master df, heatmap df
        if(length(strsplit(val2,"\\-")[[1]]) >1){
          temp_row <- c(i,val,val2,
                        coef(summary(temp_model))[paste("`",val2,"`",sep=""),1],
                        coef(summary(temp_model))[paste("`",val2,"`",sep=""),1] - coef(summary(temp_model))[paste("`",val2,"`",sep=""),2],
                        coef(summary(temp_model))[paste("`",val2,"`",sep=""),1] + coef(summary(temp_model))[paste("`",val2,"`",sep=""),2],
                        coef(summary(temp_model))[paste("`",val2,"`",sep=""),4])
        }
        if(length(strsplit(val2,"\\-")[[1]])<2){
          temp_row <- c(i,val,val2,
                        coef(summary(temp_model))[val2,1],
                        coef(summary(temp_model))[val2,1] - coef(summary(temp_model))[val2,2],
                        coef(summary(temp_model))[val2,1] + coef(summary(temp_model))[val2,2],
                        coef(summary(temp_model))[val2,4])
        }
        names(temp_row) <- c("tissue","trait","protein","beta","low","high","pval")
        df_effectsizes_master <- rbind(df_effectsizes_master,temp_row)
      }
    }
  }
  colnames(df_effectsizes_master) <- c("tissue","trait","protein","beta","low","high","pval")
  #FDR Correction for ES DF
  df_effectsizes_master <- df_effectsizes_master[order(df_effectsizes_master$pval),]
  df_effectsizes_master$FDR <- p.adjust(df_effectsizes_master$pval, method = "fdr", n = nrow(df_effectsizes_master))
  ##translate uniprot ID to gene symbol for easy graphing later
  if(study=="B"){
    df_effectsizes_master$symbol <- df_sencat_blsa$unique_symbol[match(df_effectsizes_master$protein,
                                                                       df_sencat_blsa$unique_gene)]
  }
  if(study=="I"){
    df_effectsizes_master$symbol <- df_sencat_inchianti$unique_symbol[match(df_effectsizes_master$protein,
                                                                            df_sencat_inchianti$unique_gene)]
  }
  return(list("effectsize_master" = df_effectsizes_master))
}