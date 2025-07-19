

##Function that accepts a list with groups of previously selected features, 
#determines either a composite score or a seno-age gap, than finds effect size for a list of parameters
effectsize_composite <- function(ensp_list,clin_df,trait_list,cat_control, cont_control,type="lm",composite=TRUE,interaction=FALSE){
  
  #create df to hold effect sizes for all tissues/traits (heat map)
  df_heatmap_master <- as.data.frame(matrix(NA,nrow=length(names(ensp_list)),ncol=length(trait_list)))
  colnames(df_heatmap_master) <- trait_list
  rownames(df_heatmap_master) <- names(ensp_list)
  
  #df for p-values of effect sizes (for heat map)
  df_heatmap_master_p <- as.data.frame(matrix(NA,nrow=length(names(ensp_list)),ncol=length(trait_list)))
  colnames(df_heatmap_master_p) <- trait_list
  rownames(df_heatmap_master_p) <- names(ensp_list)
  
  ##DF to hold effect sizes with confidence intervals
  df_effectsizes_master <- as.data.frame(matrix(NA,nrow=0,ncol=6))
  colnames(df_effectsizes_master) <- c("tissue","trait","beta","low","high","pval")
  
  #cycle through tissue
  for(i in names(ensp_list)){
    clin_df_i <- clin_df
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
    #if composite is false, use seno-age gap
    if(composite == FALSE){
      #find ensps
      temp_ensp_table <- ensp_list[[i]]
      temp_ensps <- temp_ensp_table[,"age"][!is.na(temp_ensp_table[,"age"])]
      ##make seno-age gap from age-ensps
      clin_df_i$age <- as.numeric(scale(clin_df_i[,"age"]))
      clin_df_i[,temp_ensps] <- scale(clin_df_i[,temp_ensps])
      form <- as.formula(paste("age ~ `",paste(temp_ensps,collapse="` + `"),"`",sep=""))
      model <- lm(form,data=clin_df_i)
      clin_df_i[,i] <- as.numeric(predict(model, newdata = clin_df_i))
      #find regression line for age and seno-age
      form <- as.formula(paste("`",i,"` ~ age",sep=""))
      temp_model <- lm(form,data=clin_df_i)
      #convert seno-age to age gap
      clin_df_i[,i] <- resid(temp_model)
    }
    ##cycle through traits
    for(val in trait_list){
      #subset to rows with trait present
      clin_df_i <- clin_df_i[!is.na(clin_df_i[,val]),]
      #scale if lm
      if(type == "lm"){
        clin_df_i[,val] <- as.numeric(scale(clin_df_i[,val]))
      }
      #factor if glm
      if(type == "glm"){
        clin_df_i[,val] <- factor(clin_df_i[,val])
      }
      #if composite == True, add composite score of each tissue onto df
      if(composite == TRUE){
        temp_table <- ensp_list[[i]]
        temp_ensp <- temp_table[,val][!is.na(temp_table[,val])]
        #check for length of ensp
        if(length(temp_ensp)>1){
          #add composite score for that tissue
          clin_df_i[,i] <- scale(rowMeans(clin_df_i[,temp_ensp]))
        }
      }
      #proceed if composite or seno-age-gap has been added (ie there were enough ensps)
      if(i %in% colnames(clin_df_i)){
        #building formula for modeling
        if(!is.na(cont_control)){
          #if finding an interaction term
          if(interaction == TRUE){
            IV <- paste(" ~ age:`",i,"` + `",paste(c(i,cat_control,cont_control),collapse="` + `"),"`",sep="")
          }
          if(interaction == FALSE){
            #if finding the effect size of the predictor on the dependent variable
            IV <- paste(" ~ `",paste(c(i,cat_control,cont_control),collapse="` + `"),"`",sep="")
          }
        }
        if(is.na(cont_control)){
          #if finding an interaction term
          if(interaction == TRUE){
            IV <- paste(" ~ age:`",i,"` + `",paste(c(i,cat_control),collapse="` + `"),"`",sep="")
          }
          #if finding the effect size of the predictor on the dependent variable 
          if(interaction == FALSE){
            IV <- paste(" ~ `",paste(c(i,cat_control),collapse="` + `"),"`",sep="")
          }
        }
        temp_form <- formula(paste("`",val,"`",IV,sep=""))
        #build model 
        #for continuous outcomes 
        if(type == "lm"){
          temp_model <- lm(temp_form, data=clin_df_i)
        }
        #for binary outcomes 
        if(type == "glm"){
          temp_model <- glm(temp_form, data=clin_df_i, family="binomial")
          #temp_model <- glm(temp_form, data=clin_df_i, family=gaussian(link="identity"))
        }
        #add data to master df, heatmap df
        #if using the effect size of the predictor on dependent variable 
        if(interaction == FALSE){
          if(length(strsplit(i,"\\ ")[[1]]) >1){
            temp_row <- c(i,val,
                          coef(summary(temp_model))[paste("`",i,"`",sep=""),1],
                          coef(summary(temp_model))[paste("`",i,"`",sep=""),1] - coef(summary(temp_model))[paste("`",i,"`",sep=""),2],
                          coef(summary(temp_model))[paste("`",i,"`",sep=""),1] + coef(summary(temp_model))[paste("`",i,"`",sep=""),2],
                          coef(summary(temp_model))[paste("`",i,"`",sep=""),4])
            #add data to heatmaps 
            df_heatmap_master[i,val] <- coef(summary(temp_model))[paste("`",i,"`",sep=""),1]
            df_heatmap_master_p[i,val] <- coef(summary(temp_model))[paste("`",i,"`",sep=""),4]
          }
          if(length(strsplit(i,"\\ ")[[1]])<2){
            temp_row <- c(i,val,
                          coef(summary(temp_model))[i,1],
                          coef(summary(temp_model))[i,1] - coef(summary(temp_model))[i,2],
                          coef(summary(temp_model))[i,1] + coef(summary(temp_model))[i,2],
                          coef(summary(temp_model))[i,4])
            df_heatmap_master[i,val] <- coef(summary(temp_model))[i,1]
            df_heatmap_master_p[i,val] <- coef(summary(temp_model))[i,4]
          }
        }
        #if finding the interaction term
        if(interaction == TRUE){
          if(length(strsplit(i,"\\ ")[[1]])<2){
            temp_row <- c(i,val,
                          coef(summary(temp_model))[paste("age:",i,sep=""),1],
                          coef(summary(temp_model))[paste("age:",i,sep=""),1] - coef(summary(temp_model))[paste("age:",i,sep=""),2],
                          coef(summary(temp_model))[paste("age:",i,sep=""),1] + coef(summary(temp_model))[paste("age:",i,sep=""),2],
                          coef(summary(temp_model))[paste("age:",i,sep=""),4])
            #add data to heatmaps 
            df_heatmap_master[i,val] <- coef(summary(temp_model))[paste("age:",i,sep=""),1]
            df_heatmap_master_p[i,val] <- coef(summary(temp_model))[paste("age:",i,sep=""),4]
          }
          if(length(strsplit(i,"\\ ")[[1]])>1){
            temp_row <- c(i,val,
                          coef(summary(temp_model))[paste("age:`",i,"`",sep=""),1],
                          coef(summary(temp_model))[paste("age:`",i,"`",sep=""),1] - coef(summary(temp_model))[paste("age:`",i,"`",sep=""),2],
                          coef(summary(temp_model))[paste("age:`",i,"`",sep=""),1] + coef(summary(temp_model))[paste("age:`",i,"`",sep=""),2],
                          coef(summary(temp_model))[paste("age:`",i,"`",sep=""),4])
            #add data to heatmaps 
            df_heatmap_master[i,val] <- coef(summary(temp_model))[paste("age:`",i,"`",sep=""),1]
            df_heatmap_master_p[i,val] <- coef(summary(temp_model))[paste("age:`",i,"`",sep=""),4]
          }
        }
        names(temp_row) <- c("tissue","trait","beta","low","high","pval")
        df_effectsizes_master <- rbind(df_effectsizes_master,temp_row)
      }
    }
  }
  colnames(df_effectsizes_master) <- c("tissue","trait","beta","low","high","pval")
  #FDR Correction for ES DF
  df_effectsizes_master <- df_effectsizes_master[order(df_effectsizes_master$pval),]
  df_effectsizes_master$FDR <- p.adjust(df_effectsizes_master$pval, method = "fdr", n = nrow(df_effectsizes_master))
  #add star column for heatmap
  df_effectsizes_master_plot <- df_effectsizes_master
  df_effectsizes_master_plot$star <- ifelse(df_effectsizes_master_plot$FDR < 0.05, "*","")
  colnames(df_effectsizes_master_plot) <- c("Tissue","Trait","beta","low","high","pval","fdr","star")
  ##plot heatmap using geom_tile
  heatmap_graph <- ggplot(df_effectsizes_master_plot, aes(x=Tissue, y=Trait, fill=as.numeric(beta)))+
    geom_tile()+
    scale_fill_gradient2(low=unname(brads_heatmap_palette)[1],mid="white",high=unname(brads_heatmap_palette)[3],midpoint=0, na.value = "white")+
    geom_text(aes(label=star), color="black", size=6)+
    theme_minimal(base_size=16)+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.title.align=0.5,
          panel.background=element_rect(fill="white", color=NA),
          plot.background=element_rect(fill="white", color=NA),
          panel.grid=element_blank())+
    labs(fill="Effect\nSize", x="Cell Type")
  ##FDR Correction for heatmap p-values
  for(i in rownames(df_heatmap_master_p)){
    if(i %in% df_effectsizes_master$tissue){
      for(val in colnames(df_heatmap_master_p)){
        if(val %in% df_effectsizes_master[df_effectsizes_master$tissue ==i,"trait"]){
          df_heatmap_master_p[i,val] <- df_effectsizes_master[df_effectsizes_master$tissue == i & df_effectsizes_master$trait ==val,"FDR"]
        }
      }
    }
  }
  return(list("effectsize_master" = df_effectsizes_master,"heat" = df_heatmap_master,"heat_p"=df_heatmap_master_p, heat_graph=heatmap_graph))
}