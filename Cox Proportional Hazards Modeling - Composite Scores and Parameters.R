

##function that accepts a list of features, performs cox proportional hazards modeling using composite scores and parameters
cox_modeling_trait_celltype <- function(df, sap_list, disease_list, ensps = 25) {
  
  # Ensure df is a tibble for compatibility with dplyr
  df <- as_tibble(df)
  
  # Initialize an empty data frame to store Cox model results
  results_df <- tibble(
    tissue = character(),
    disease = character(),
    effect_size = numeric(),
    low_ci = numeric(),
    high_ci = numeric(),
    pval = numeric()
  )
  
  # Initialize an empty data frame to store selected features for each trait and tissue
  selected_ensps_df <- tibble(
    tissue = character(),
    disease = character(),
    ensps = character()
  )
  
  #make df to hold shoenfeld stats per composite score/parameter combination 
  ensp_shoenfeld_stats <- tibble(
    tissue = character(),
    disease = character(),
    chisq = numeric(),
    pval = numeric()
  )
  
  # Iterate through each tissue in the sap_list
  for (tissue in names(sap_list)) {
    
    # Pull the proteins specific to the current tissue
    proteins <- sap_list[[tissue]]
    
    #translate proteins to BLSA nomenclature
    proteins <- df_sencat_blsa[df_sencat_blsa$gene %in% proteins,"unique_gene"]
    
    #add controls
    proteins <- c(proteins,c("sex","race"))
    
    # Iterate through each disease in the disease_list
    for (disease in disease_list) {
      
      # Use sym() around the disease name to handle names with spaces
      disease_col <- sym(disease)
      
      # Filter out cases where disease occurred prior to proteomics measurement
      incident_cases <- df %>%
        dplyr::filter(!is.na(!!disease_col)) %>%
        group_by(IDNo) %>%
        filter(!any(Visit <= SAP_Visit & !!disease_col == 1)) %>%
        summarize(
          age_followup = ifelse(any(!!disease_col == 1), min(age[!!disease_col == 1]), max(age)),
          event = ifelse(any(!!disease_col == 1), 1, 0)
        ) %>%
        ungroup()
      
      #select the first instance of the proteome for each before joining with incident_cases
      tissue_data <- df %>%
        dplyr::select(IDNo, SAP_Age, all_of(proteins)) %>%
        filter(!duplicated(IDNo))
      
      #merge incident cases with proteins and controls
      incident_cases <- merge(incident_cases, tissue_data, by="IDNo", all.x=TRUE)
      
      #filter out individuals who's last visit was their SAP measurement
      incident_cases <- incident_cases %>%
        filter(age_followup > SAP_Age)
      
      # Prepare matrix X and response y for cv.glmnet
      X <- as.matrix(incident_cases %>% dplyr::select(all_of(proteins)))
      y <- Surv(incident_cases$SAP_Age,incident_cases$age_followup, incident_cases$event)
      
      #set seed, define folds for reproducibility
      set.seed(1)
      foldid_cox <- sample(1:10, size = nrow(y), replace = TRUE)
      
      # Run Elastic Net using cv.glmnet
      cv_model <- cv.glmnet(X, y, family = "cox", alpha = 0.5, standardize=TRUE, parallel=TRUE, folid=foldid_cox)
      
      # Extract the coefficients of selected proteins
      best_lambda <- cv_model$lambda.min
      selected_proteins <- coef(cv_model, s = best_lambda)
      
      #reorder in descending order
      selected_proteins <- as.data.frame(as.matrix(selected_proteins))
      selected_proteins$gene <- rownames(selected_proteins)
      selected_proteins <- selected_proteins[selected_proteins[,1] > 0,]
      
      #remove controls
      selected_proteins <- selected_proteins[!selected_proteins$gene %in% c("age","sex","race"),]
      
      #skip iteration if no ensps
      if(nrow(selected_proteins) <2) next 
      
      #add dummy column for reordering 
      selected_proteins <- selected_proteins[rev(order(selected_proteins[,1])),]
      
      # Limit to the top ensps
      selected_proteins <- head(rownames(selected_proteins), ensps)
      
      # Check if more than one selected protein for composite score, otherwise skip
      if (length(selected_proteins) <2) next
      
      # Add selected ensps for this tissue and disease to the ensps data frame
      selected_ensps_df <- bind_rows(selected_ensps_df, tibble(
        tissue = tissue,
        disease = disease,
        ensps = paste(selected_proteins, collapse = ", ")
      ))
      
      # Calculate composite score as the scaled mean of selected proteins
      incident_cases$composite_score <- scale(rowMeans(incident_cases[,selected_proteins]))
      
      # Perform Cox model with composite score
      cox_model <- coxph(Surv(SAP_Age,age_followup, event) ~ composite_score + strata(factor(sex)) + strata(factor(race)), 
                         data = incident_cases)
      summary_cox <- summary(cox_model)
      # Find hazard ratio of composite score
      effect_size <- summary_cox$coefficients["composite_score", "exp(coef)"]
      pval <- summary_cox$coefficients["composite_score", "Pr(>|z|)"]
      ci_low <- summary_cox$conf.int["composite_score", "lower .95"]
      ci_high <- summary_cox$conf.int["composite_score", "upper .95"]
      # Update results data frame
      results_df <- bind_rows(results_df, tibble(
        tissue = tissue,
        disease = disease,
        effect_size = effect_size,
        low_ci = ci_low,
        high_ci = ci_high,
        pval = pval
      ))
      #make and fill shoenfeld stats 
      shoenfeld_model <- cox.zph(cox_model)
      ensp_shoenfeld_stats <- bind_rows(ensp_shoenfeld_stats, tibble(
        tissue = tissue,
        disease = disease,
        chisq = shoenfeld_model$table["composite_score","chisq"],
        pval = shoenfeld_model$table["composite_score","p"]
      ))
      
    }
  }
  
  # Return both data frames as a list
  return(list(results = results_df, selected_ensps = selected_ensps_df,shoenfeld=ensp_shoenfeld_stats))
}