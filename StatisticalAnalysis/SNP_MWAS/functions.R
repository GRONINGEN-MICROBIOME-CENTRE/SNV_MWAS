library(tidyverse)
library(haven)
library(readxl)
library(scales)
library(caret)
library(data.table)
library(pROC)
library(car)
library(tidymodels)
library(bonsai)
library(lightgbm)
library(plyr)
library(magrittr)
library(coin)
library(permute)
library(beepr)
library(plyr)
library(basicPlotteR)
library(psych)
library(ppcor)
library(coda4microbiome)
library(mediation)
library(robustbase)

library(wesanderson)
library(ggsci)
library(ggridges)
library(ggpubr)
library(ggrepel)
library(ggpie)
library(ggalluvial)
library(ggseqlogo)
library(cowplot)
library(ComplexHeatmap)
library(ComplexUpset)
library(viridis)
library(circlize)
library(RColorBrewer)
library(gridBase)
#library(rgl)
#library(scatterplot3d)
#library(plot3D)

# Color settings
mycolor2_blue_red<-c("#05B2DC", "#F9627D")
mycolor2_green_blue  <- c("#2EC4B6","#235789")
viridis<-c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725","white") %>% rev

# Cutoffs
studyWideP <- 3.01E-11
metaWideP <- 7.34E-9
speWideP   <- 6.47E-6

sigLevel <- c("StudyWideSig", "MetaWideSig", "SpeWideSig")
sigCol   <- c("red", "orange", "grey50")

# Annotations
annotationType <- c("stop_gained","start_lost","stop_lost&splice_region_variant",
                    "missense_variant", "splice_region_variant&stop_retained_variant", "start_retained_variant","initiator_codon_variant",
                    "upstream_gene_variant", "downstream_gene_variant", "synonymous_variant")
annotationLevel <- c("High", "High", "High",
                     "Moderate", "Moderate","Moderate","Moderate",
                     "Low","Low","Low")
annotationImpact <- c("High", "High", "High",
                      "Moderate", "Moderate","Moderate","Moderate",
                      "Low","Low","Low") 
annotationCol <- c("#E64B35", "#E64B35", "#E64B35",
                   "#00A087", "#00A087","#00A087","#00A087",
                   "#3C5488","#3C5488","#3C5488") 
impactLevelCol <- c("#E64B35", "#00A087", "#3C5488")
# Phenotypes
phenClassOrder <- c("Anthropology", "CurrentExposure", "Diet", "EarlyLifeExposure",
                    "Stool", "GeographySocioeconomics", "HealthDisease", "MedicalMeasurement", "MedicationUse", "Metabolite")
phenClassColor <- c("#257CB2", "#E93E3F", "#4F9F3B", "#FE8B15",
                    "#10aec2", "#A6CEE3", "#F68181", "#A5D981", "#F7F599", "#73489F")

# Cohort
cohort_order       <- c("dmp", "lld1", "lld2_apk","lld2_fsk", "500fg_apk", "500fg_fsk", "300ob", "ibd", "300tzfg")
cohort_label       <- c("DMP", "LLD1", "LLD2-APK","LLD2-FSK", "500FG-APK", "500FG-FSK", "300OB", "IBD", "300TZFG")
cohort_title_label <- c("DMP", "LLD1", "LLD2_APK","LLD2_FSK", "500FG_APK", "500FG_FSK", "300OB", "IBD", "300TZFG")

cohort_color<-c(pal_material("grey")(6)[c(6)],                 # 300TZFG
                pal_material("red")(6)[c(5)],                  # IBD
                pal_material("amber")(6)[c(5)],                # 300OOB
                pal_material("indigo")(6)[c(3,5)],             # 500FG-FSK 500FG-APK
                pal_material("cyan")(6)[c(5)],                 # LLD2-FSK
                pal_material("light-blue")(6)[c(3,5)],         # LLD2-APK LLD1
                pal_material("light-green")(6)[c(5)]) %>% rev  # DMP

cohort_order_6       <- c("dmp", "lld1",  "500fg",  "300ob", "ibd", "300tzfg")
cohort_label_6       <- c("DMP", "LLD1",  "500FG",  "300OB", "IBD", "300TZFG")
cohort_title_label_6 <- c("DMP", "LLD1",  "500FG",  "300OB", "IBD", "300TZFG")

cohort_color_6<-c(pal_material("grey")(6)[c(6)],                 # 300TZFG
                pal_material("red")(6)[c(5)],                  # IBD
                pal_material("amber")(6)[c(5)],                # 300OOB
                pal_material("indigo")(6)[c(3)],             # 500FG-FSK
                pal_material("light-blue")(6)[c(5)],         # LLD1
                pal_material("light-green")(6)[c(5)]) %>% rev  # DMP


## Calculate SE
se <- function(x) {x<-na.omit(x);sqrt(var(x)/length(x))}

## The Normal Quantile Transformation
qtrans<-function(x){
  k<-!is.na(x)
  k<-which(x!="-999")
  ran<-rank(as.numeric(x[k]))
  y<-qnorm((1:length(k)-0.5)/length(k))
  x[k]<-y[ran]
  x
}

transform_and_filter_mtb=function(x, samples_row=T, method="asin", missing_filter=0){
  #x[x=="NA"]=0
  #x[x=="NaN"]=0
  #if samples are in columns transpose
  if (!samples_row){
    
    x=as.data.frame(t(x))
    print("transposing matrix!")
    
  } 
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>100){
    
    stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 100") 
    
  }
  
  #x_filt=x[,((colSums(x !=0) / nrow(x)) *100 )>missing_filter]
  x_filt=x[,((colSums(!is.na(x)) / nrow(x)) *100 )>missing_filter]
  my_num_removed=ncol(x)-ncol(x_filt)
  print (paste(my_num_removed, "species removed due to many missing values"))
  
  if (method=="asin"){
    print ("ASIN")
    x_filt=x_filt/100
    x_filt=asin(sqrt(x_filt))
    
  } else if (method=="log"){
    print ("LOG10")
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0])/2
    x_filt=x_filt+my_min
    x_filt=log10(x_filt)
    
  }else if (method=="clr"){
    print ("CLR")
    #Adapted from Alexander Kurilshikov 
    #x_filt=x_filt/100
    #replace 0 by the half of the smallest value observed
    #my_min=min(x[x>0], na.rm=T)/2
    #x=x+my_min
    #Calculate geometric mean
    gm_mean = function(x, na.rm=TRUE){
      exp(mean(log(x),na.rm=T))
    }
    Gmean_core = apply(x, 1, gm_mean)
    x = cbind(Gmean_core,x)
    d <- t(apply(x, 1, function(b) {
      log(b/b[1])[-1]
    }))
    x=d
    x_filt=x[,colnames(x) %in%colnames(x_filt)]
  }
  return(as.data.frame(x_filt))
}

## Get density of 2-demision dataset
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

## associations between two matrices using linear model
lm_btw_mats<-function(y_mat,x_mat,cov_mat,covar,convert.x = T){
  require(reshape2)
  require(R.utils)
  
  ## test block
  # y_mat<- y_mat_i #lld_intri[,c(1:3)]
  # x_mat<- x_mat #lld_vsv[,c(1:3)]
  # cov_mat<- cov_mat# lld_covar
  # covar<- covar_i #covar
  ## test block
  
  my_lm<-function(y,x){
    # y<-y_mat[,1]
    # x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    #    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    lm_input<-as.data.frame(lm_input)
    lm_input$Y<-qtrans(lm_input$Y)
    if(convert.x){
      lm_input$X<-qtrans(lm_input$X)
    }
    
    try(lm_res <- summary(lm(Y~.,data = lm_input)), silent = T)
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 10, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}


## Machine learning model for regression task based on lightgbm and nested cross-validation
## Including model parameter optimization and feature selection, takes time
my_gbt_reg.ncv <- function(input, Sample, Outcome, Feature, treeN = 2000, outside.v = 5, outside.repeat = 5, inside.v = 5, bayes.iter = 20){
  ## test
  # input   <- train
  # Sample  <- rownames(train)
  # Outcome <- target
  # Feature <- feature
  # treeN   <- 2000
  # outside.v = 2
  # outside.repeat = 2
  # inside.v = 2
  # bayes.iter = 3
  ## test
  
  require(tidyverse)
  require(tidymodels)
  require(bonsai)
  require(lightgbm)
  require(data.table)
  
  # clean data
  input <- data.frame(Outcome = input[Sample,Outcome], input[Sample, Feature], check.names = T)
  input$Outcome<-as.numeric(input$Outcome)
  input <- input[!is.na(input$Outcome),]
  
  set.seed(42)
  folds <- nested_cv(input,
                     outside = vfold_cv(v = outside.v, repeats = outside.repeat, strata = Outcome),
                     inside  = vfold_cv(v = inside.v, strata = Outcome))
  
  res <- NULL
  for (i in 1:length(folds$splits)) {
    ## outside loop
    #    i <- 1
    cat(paste("Round",i, "started on",date(),"\n"))
    
    # Testset split
    input_train <- training(folds$splits[[i]])
    input_test  <- testing(folds$splits[[i]])
    
    # recipe
    input_train_rec <- 
      recipe(Outcome ~ ., data = folds$splits[[i]]) %>%
      step_zv(all_predictors())
    
    # Pre-training: construct model
    boost_mod <- boost_tree(trees = treeN,
                            tree_depth = tune(), min_n = tune(),
                            loss_reduction = tune(), ## first three: model complexity
                            sample_size = tune(), 
                            mtry = tune(),           ## randomness
                            learn_rate = tune()) %>% 
      set_engine("lightgbm") %>% set_mode("regression")
    boost_wf <- workflow() %>% add_model(boost_mod) %>% add_recipe(input_train_rec)
    
    # Pre-training: Train with CV and bayes optimization, tune parameters and select features
    boost_set<-extract_parameter_set_dials(boost_wf) %>% finalize(input[,-1])
    
    doParallel::registerDoParallel()
    set.seed(42)
    boost_fit_rs <- tune_bayes(
      boost_wf,
      resamples = folds$inner_resamples[[i]],
      param_info = boost_set,
      # Generate five at semi-random to start
      initial = nrow(boost_set) * 2,
      iter = bayes.iter,
      # How to measure performance?
      metrics = metric_set(rsq, rmse),
      control = control_bayes(no_improve = 5, verbose = T)
    )
    
    # Final best model
    best_boost <- select_best(boost_fit_rs, "rsq") # Best parameter within current inner CV
    final_wf <- boost_wf %>% finalize_workflow(best_boost) # Best model within cirrent inner CV
    final_fitted_wf <- final_wf %>% fit(data=input_train)
    
    # Independent test and assessment
    final_fit <- final_wf %>% last_fit(folds$splits[[i]]) 
    final_perf <- collect_metrics(final_fit)$.estimate # Final performance of current inner CV
    cat(paste("Round",i, "performance",final_perf,"\n"))
    
    # Final prediction
    final_prediction<-bind_cols(
      SampleID = rownames(input),
      SampleRole = "Test",
      ActualClass = input$Outcome,
      predict(final_fitted_wf, input)
    )
    final_prediction$SampleRole[folds$splits[[i]]$in_id]<-"Train"
    
    
    # Get feature importance
    shp <- shapviz::shapviz(extract_fit_engine(final_fitted_wf), X_pred = as.matrix(input_train[,final_fitted_wf$pre$mold$predictors%>%colnames()]))
    shp.df <- shapviz::sv_importance(shp,max_display=Inf)
    importance <- data.frame(Feature = shp.df$data$feature,
                             Importance = shp.df$data$value)
    
    res[[i]]<-list(best.para          = best_boost,
                   best.model         = final_fitted_wf, 
                   final.fit          = final_fit, 
                   final.performance  = final_perf, 
                   final_prediction   = final_prediction,
                   feature.importance = importance)
    
    cat(paste("Round",i, "finished on",date(),"\n"))
  }
  
  names(res) <- paste(folds$id, folds$id2, sep = "_")
  
  # Decide parameters for final model
  merged_best_para <- res %>%
    map(~ .x$best.para) %>%
    bind_rows() %>%
    summarise_all(median)
  merged_best_para$.config <- 'Final'
  
  # Select features for final model
  cat(paste("Feature selection", "started on",date(),"\n"))
  feature.long <- NULL
  for (i in 1:length(res)) {
    # i<-1
    feature_i <- res[[i]]$feature.importance
    feature_i$Fold <- names(res)[i]
    feature.long <- rbind(feature.long, feature_i)
  }
  feature.final.importance <- feature.long %>%
    group_by(Feature) %>%
    summarise(FinalImportance = median(Importance, na.rm = TRUE)) %>%
    dplyr::arrange(desc(FinalImportance)) %>%
    as.data.frame()
  feature.final.importance$Feature <- as.character(feature.final.importance$Feature)
  feature.final.importance$Cumsum <- cumsum(feature.final.importance$FinalImportance)
  feature.final.importance.sum <- sum(feature.final.importance$FinalImportance)
  feature.final.importance.cutoff <- feature.final.importance.sum*0.95
  feature.number <- 1:nrow(feature.final.importance) %>% .[feature.final.importance$Cumsum>feature.final.importance.cutoff] %>% .[1]
  feature.final.importance$CheckCutoff <- feature.final.importance$Cumsum > feature.final.importance.cutoff
  final.feature <- feature.final.importance$Feature[1:feature.number] %>% as.character()
  cat(paste("Feature selection", "finished on",date(),"\n"))
  
  
  # Final model
  cat(paste("Train final model", "started on",date(),"\n"))
  
  final_input_train_rec <- 
    recipe(Outcome ~ ., data = input[,c("Outcome", final.feature)]) %>%
    step_zv(all_predictors())
  final_boost_wf <- workflow() %>%
    add_model(boost_mod) %>% 
    add_recipe(final_input_train_rec) %>%
    finalize_workflow(best_boost) # Best model within cirrent inner CV
  final2_fitted_wf <- final_boost_wf %>% fit(data=input)
  cat(paste("Train final model", "finished on",date(),"\n"))
  
  # Final result
  final_res <- list(cv = res,
                    ready.model = final2_fitted_wf,
                    features = feature.final.importance,
                    parameters = merged_best_para)
  return(final_res)
}


## Machine learning model for classification task based on lightgbm and nested cross-validation
## Including model parameter optimization and feature selection, takes time
my_gbt_class.ncv <- function(input, Sample, Outcome, Feature, treeN = 2000, outside.v = 5, outside.repeat = 5, inside.v = 5, bayes.iter = 20){
  ## test
  # input   <- train
  # Sample  <- rownames(train)
  # Outcome <- target
  # Feature <- feature
  # treeN   <- 2000
  # outside.v = 2
  # outside.repeat = 2
  # inside.v = 2
  # bayes.iter = 3
  ## test
  
  require(tidyverse)
  require(tidymodels)
  require(bonsai)
  require(lightgbm)
  require(data.table)
  
  # clean data
  input <- data.frame(Outcome = input[Sample,Outcome], input[Sample, Feature], check.names = T)
  input$Outcome<-as.numeric(input$Outcome)
  input <- input[!is.na(input$Outcome),]
  input$Outcome <- factor(input$Outcome, levels = c(0, 1))
  
  
  set.seed(42)
  folds <- nested_cv(input,
                     outside = vfold_cv(v = outside.v, repeats = outside.repeat, strata = Outcome),
                     inside  = vfold_cv(v = inside.v, strata = Outcome))
  
  res <- NULL
  for (i in 1:length(folds$splits)) {
    ## outside loop
    #    i <- 1
    cat(paste("Round",i, "started on",date(),"\n"))
    
    # Testset split
    input_train <- training(folds$splits[[i]])
    input_test  <- testing(folds$splits[[i]])
    
    # recipe
    input_train_rec <- 
      recipe(Outcome ~ ., data = folds$splits[[i]]) %>%
      step_zv(all_predictors())
    
    # Pre-training: construct model
    boost_mod <- boost_tree(trees = treeN,
                            tree_depth = tune(), min_n = tune(),
                            loss_reduction = tune(), ## first three: model complexity
                            sample_size = tune(), 
                            mtry = tune(),           ## randomness
                            learn_rate = tune()) %>% 
      set_engine("lightgbm") %>% set_mode("classification")
    boost_wf <- workflow() %>% add_model(boost_mod) %>% add_recipe(input_train_rec)
    
    # Pre-training: Train with CV and bayes optimization, tune parameters and select features
    boost_set<-extract_parameter_set_dials(boost_wf) %>% finalize(input[,-1])
    
    doParallel::registerDoParallel()
    set.seed(42)
    boost_fit_rs <- tune_bayes(
      boost_wf,
      resamples = folds$inner_resamples[[i]],
      param_info = boost_set,
      # Generate five at semi-random to start
      initial = nrow(boost_set) * 2,
      iter = bayes.iter,
      # How to measure performance?
      metrics = metric_set(roc_auc),
      control = control_bayes(no_improve = 5, verbose = T)
    )
    
    # Final best model
    best_boost <- select_best(boost_fit_rs, "roc_auc") # Best parameter within current inner CV
    final_wf <- boost_wf %>% finalize_workflow(best_boost) # Best model within cirrent inner CV
    final_fitted_wf <- final_wf %>% fit(data=input_train)
    
    # Independent test and assessment
    final_fit <- final_wf %>% last_fit(folds$splits[[i]]) 
    final_perf <- collect_metrics(final_fit)$.estimate # Final performance of current inner CV
    cat(paste("Round",i, "performance",final_perf,"\n"))
    
    # Final prediction
    final_prediction<-bind_cols(
      SampleID = rownames(input),
      SampleRole = "Test",
      ActualClass = input$Outcome,
      predict(final_fitted_wf, input)
    )
    final_prediction$SampleRole[folds$splits[[i]]$in_id]<-"Train"
    
    
    # Get feature importance
    shp <- shapviz::shapviz(extract_fit_engine(final_fitted_wf), X_pred = as.matrix(input_train[,final_fitted_wf$pre$mold$predictors%>%colnames()]))
    shp.df <- shapviz::sv_importance(shp,max_display=Inf)
    importance <- data.frame(Feature = shp.df$data$feature,
                             Importance = shp.df$data$value)
    
    res[[i]]<-list(best.para          = best_boost,
                   best.model         = final_fitted_wf, 
                   final.fit          = final_fit, 
                   final.performance  = final_perf, 
                   final_prediction   = final_prediction,
                   feature.importance = importance)
    
    cat(paste("Round",i, "finished on",date(),"\n"))
  }
  
  names(res) <- paste(folds$id, folds$id2, sep = "_")
  
  # Decide parameters for final model
  merged_best_para <- res %>%
    map(~ .x$best.para) %>%
    bind_rows() %>%
    summarise_all(median)
  merged_best_para$.config <- 'Final'
  
  # Select features for final model
  cat(paste("Feature selection", "started on",date(),"\n"))
  feature.long <- NULL
  for (i in 1:length(res)) {
    # i<-1
    feature_i <- res[[i]]$feature.importance
    feature_i$Fold <- names(res)[i]
    feature.long <- rbind(feature.long, feature_i)
  }
  feature.final.importance <- feature.long %>%
    group_by(Feature) %>%
    summarise(FinalImportance = median(Importance, na.rm = TRUE)) %>%
    dplyr::arrange(desc(FinalImportance)) %>%
    as.data.frame()
  feature.final.importance$Feature <- as.character(feature.final.importance$Feature)
  feature.final.importance$Cumsum <- cumsum(feature.final.importance$FinalImportance)
  feature.final.importance.sum <- sum(feature.final.importance$FinalImportance)
  feature.final.importance.cutoff <- feature.final.importance.sum*0.95
  feature.number <- 1:nrow(feature.final.importance) %>% .[feature.final.importance$Cumsum>feature.final.importance.cutoff] %>% .[1]
  feature.final.importance$CheckCutoff <- feature.final.importance$Cumsum > feature.final.importance.cutoff
  final.feature <- feature.final.importance$Feature[1:feature.number] %>% as.character()
  cat(paste("Feature selection", "finished on",date(),"\n"))
  
  # Final model
  cat(paste("Train final model", "started on",date(),"\n"))
  
  final_input_train_rec <- 
    recipe(Outcome ~ ., data = input[,c("Outcome", final.feature)]) %>%
    step_zv(all_predictors())
  final_boost_wf <- workflow() %>%
    add_model(boost_mod) %>% 
    add_recipe(final_input_train_rec) %>%
    finalize_workflow(best_boost) # Best model
  final2_fitted_wf <- final_boost_wf %>% fit(data=input)
  cat(paste("Train final model", "finished on",date(),"\n"))
  
  # Final result
  final_res <- list(cv = res,
                    ready.model = final2_fitted_wf,
                    features = feature.final.importance,
                    parameters = merged_best_para)
  return(final_res)
}

## Machine learning model for regression task based on lightgbm and nested cross-validation
## Only do feature selection,  model parameters are pre-defined by experience
my_gbt_reg.ncv.predef <- function(input, Sample, Outcome, Feature, 
                                  treeN = 2000, tree_depth = 15, min_n = 30,sample_size = 0.6,mtry = round(0.2*length(Feature)),learn_rate = 0.01,
                                  outside.v = 5, outside.repeat = 5, inside.v = 5){
  ## test
  # input   <- train.mbio.input
  # Sample  <- rownames(train.mbio.input)
  # Outcome <- "Age"
  # Feature <- mbio.feature.list
  # 
  # treeN   <- 2000
  # tree_depth = 15
  # min_n = 30
  # sample_size = 0.6
  # mtry = round(0.2*length(mbio.feature.list))
  # learn_rate = 0.01
  # 
  # outside.v = 5
  # outside.repeat = 5
  # inside.v = 5
  ## test
  
  require(tidyverse)
  require(tidymodels)
  require(bonsai)
  require(lightgbm)
  require(data.table)
  
  # clean data
  input <- data.frame(Outcome = input[Sample,Outcome], input[Sample, Feature], check.names = T)
##  input <- data.frame(Outcome = input[match(Sample, rownames(input)), match(Outcome, colnames(input))], 
#                      input[match(Sample, rownames(input)), match(Feature, colnames(input))], 
#                      check.names = T)
  
  input$Outcome<-as.numeric(input$Outcome)
  input <- input[!is.na(input$Outcome),]
  
  ### Pre-filter features
  # Feature.p <- NULL
  # for (i in 1:length(Feature)) {
  #   cor_i <- cor.test(input[,Feature[i]], input[,"Outcome"],method = "spearman")
  #   Feature.p <- c(Feature.p, cor_i$p.value)
  # }
  # input <- data.frame(Outcome = input$Outcome,input[, Feature[Feature.p < 0.05]])
  ###
  
  set.seed(42)
  folds <- nested_cv(input,
                     outside = vfold_cv(v = outside.v, repeats = outside.repeat, strata = Outcome),
                     inside  = vfold_cv(v = inside.v, strata = Outcome))
  
  res <- NULL
  for (i in 1:length(folds$splits)) {
    ## outside loop
    #    i <- 1
    cat(paste("Round",i, "started on",date(),"\n"))
    
    # Testset split
    input_train <- training(folds$splits[[i]])
    input_test  <- testing(folds$splits[[i]])
    
    # recipe
    input_train_rec <- 
      recipe(Outcome ~ ., data = folds$splits[[i]]) %>%
      step_zv(all_predictors())
    
    # Pre-training: construct model
    boost_mod <- boost_tree(trees = treeN,
                            tree_depth = tree_depth, min_n = min_n,
                            # loss_reduction = tune(),                     ## first three: model complexity
                            sample_size = sample_size, 
                            mtry = mtry,         ## randomness
                            learn_rate = learn_rate) %>% 
      set_engine("lightgbm") %>% set_mode("regression")
    boost_wf <- workflow() %>% add_model(boost_mod) %>% add_recipe(input_train_rec)
    
    # Skiped parameter optimazition here
    final_fitted_wf <- boost_wf %>% fit(data=input_train)
    
    # Independent test and assessment
    final_fit <- final_fitted_wf %>% last_fit(folds$splits[[i]]) 
    final_perf <- collect_metrics(final_fit)$.estimate # Final performance of current inner CV
    cat(paste("Round",i, "performance",final_perf,"\n"))
    
    # Final prediction
    final_prediction<-bind_cols(
      SampleID = rownames(input),
      SampleRole = "Test",
      ActualClass = input$Outcome,
      predict(final_fitted_wf, input)
    )
    final_prediction$SampleRole[folds$splits[[i]]$in_id]<-"Train"
    
    
    # Get feature importance
    shp <- shapviz::shapviz(extract_fit_engine(final_fitted_wf), X_pred = as.matrix(input_train[,final_fitted_wf$pre$mold$predictors%>%colnames()]))
    shp.df <- shapviz::sv_importance(shp,max_display=Inf)
    importance <- data.frame(Feature = shp.df$data$feature,
                             Importance = shp.df$data$value)
    
    res[[i]]<-list(
      best.model         = final_fitted_wf, 
      final.fit          = final_fit, 
      final.performance  = final_perf, 
      final_prediction   = final_prediction,
      feature.importance = importance)
    
    cat(paste("Round",i, "finished on",date(),"\n"))
  }
  
  names(res) <- paste(folds$id, folds$id2, sep = "_")
  
  # Select features for final model
  cat(paste("Feature selection", "started on",date(),"\n"))
  feature.long <- NULL
  for (i in 1:length(res)) {
    # i<-1
    feature_i <- res[[i]]$feature.importance
    feature_i$Fold <- names(res)[i]
    feature.long <- rbind(feature.long, feature_i)
  }
  feature.final.importance <- feature.long %>%
    group_by(Feature) %>%
    dplyr::summarise(FinalImportance = median(Importance, na.rm = TRUE)) %>%
    dplyr::arrange(desc(FinalImportance)) %>%
    as.data.frame()
  feature.final.importance$Feature <- as.character(feature.final.importance$Feature)
  feature.final.importance$Cumsum <- cumsum(feature.final.importance$FinalImportance)
  feature.final.importance.sum <- sum(feature.final.importance$FinalImportance)
  feature.final.importance.cutoff <- feature.final.importance.sum*0.95
  feature.number <- 1:nrow(feature.final.importance) %>% .[feature.final.importance$Cumsum>feature.final.importance.cutoff] %>% .[1]
  feature.final.importance$CheckCutoff <- feature.final.importance$Cumsum > feature.final.importance.cutoff
  final.feature <- feature.final.importance$Feature[1:feature.number] %>% as.character()
  cat(paste("Feature selection", "finished on",date(),"\n"))
  
  
  # Final model
  cat(paste("Train final model", "started on",date(),"\n"))
  
  final_input_train_rec <- 
    recipe(Outcome ~ ., data = input[,c("Outcome", final.feature)]) %>%
    step_zv(all_predictors())
  final_boost_wf <- workflow() %>%
    add_model(boost_mod) %>% 
    add_recipe(final_input_train_rec)
  final2_fitted_wf <- final_boost_wf %>% fit(data=input)
  cat(paste("Train final model", "finished on",date(),"\n"))
  
  # Final result
  final_res <- list(cv = res,
                    ready.model = final2_fitted_wf,
                    features = feature.final.importance)
  return(final_res)
}

## Machine learning model for classification task based on lightgbm and nested cross-validation
## Only do feature selection,  model parameters are pre-defined by experience
my_gbt_class.ncv.predef <- function(input, Sample, Outcome, Feature, 
                                    treeN = 1000,tree_depth = 15, min_n = 30,sample_size = 0.6,mtry = round(0.2*length(Feature)),learn_rate = 0.01,
                                    outside.v = 5, outside.repeat = 5, inside.v = 5){
  # test
  # input   <- inData
  # Sample  <- trainSample$V1
  # Outcome <- target
  # Feature <- feature$V1
  # 
  # treeN   <- 2000
  # tree_depth = 10
  # min_n = 30
  # sample_size = 0.8
  # mtry = round(0.2*ncol(input))
  # learn_rate = 0.001
  # 
  # outside.v = 2
  # outside.repeat = 2
  # inside.v = 2
  # bayes.iter = 3
  ## test
  
  require(tidyverse)
  require(tidymodels)
  require(bonsai)
  require(lightgbm)
  require(data.table)
  
  # clean data
  input <- data.frame(Outcome = input[Sample,Outcome], input[Sample, Feature], check.names = T)
  input$Outcome<-as.numeric(input$Outcome)
  input <- input[!is.na(input$Outcome),]
  input$Outcome <- factor(input$Outcome, levels = c(0, 1))
  
  
  set.seed(42)
  folds <- nested_cv(input,
                     outside = vfold_cv(v = outside.v, repeats = outside.repeat, strata = Outcome),
                     inside  = vfold_cv(v = inside.v, strata = Outcome))
  
  res <- NULL
  for (i in 1:length(folds$splits)) {
    ## outside loop
    #    i <- 1
    cat(paste("Round",i, "started on",date(),"\n"))
    
    # Testset split
    input_train <- training(folds$splits[[i]])
    input_test  <- testing(folds$splits[[i]])
    
    # recipe
    input_train_rec <- 
      recipe(Outcome ~ ., data = folds$splits[[i]]) %>%
      step_zv(all_predictors())
    
    # Pre-training: construct model
    boost_mod <- boost_tree(trees = treeN,
                            tree_depth = tree_depth, min_n = min_n,
                            # loss_reduction = tune(),                     ## first three: model complexity
                            sample_size = sample_size, 
                            mtry = mtry,         ## randomness
                            learn_rate = learn_rate) %>%
      set_engine("lightgbm") %>% set_mode("classification")
    boost_wf <- workflow() %>% add_model(boost_mod) %>% add_recipe(input_train_rec)
    final_fitted_wf <- boost_wf %>% fit(data=input_train)
    
    # Independent test and assessment
    final_fit <- final_fitted_wf %>% last_fit(folds$splits[[i]]) 
    final_perf <- collect_metrics(final_fit)$.estimate # Final performance of current inner CV
    cat(paste("Round",i, "performance",final_perf,"\n"))
    
    # Final prediction
    final_prediction<-bind_cols(
      SampleID = rownames(input),
      SampleRole = "Test",
      ActualClass = input$Outcome,
      predict(final_fitted_wf, input)
    )
    final_prediction$SampleRole[folds$splits[[i]]$in_id]<-"Train"
    
    
    # Get feature importance
    shp <- shapviz::shapviz(extract_fit_engine(final_fitted_wf), X_pred = as.matrix(input_train[,final_fitted_wf$pre$mold$predictors%>%colnames()]))
    shp.df <- shapviz::sv_importance(shp,max_display=Inf)
    importance <- data.frame(Feature = shp.df$data$feature,
                             Importance = shp.df$data$value)
    
    res[[i]]<-list(
      final.fit          = final_fit, 
      final.performance  = final_perf, 
      final_prediction   = final_prediction,
      feature.importance = importance)
    
    cat(paste("Round",i, "finished on",date(),"\n"))
  }
  
  names(res) <- paste(folds$id, folds$id2, sep = "_")
  
  # Select features for final model
  cat(paste("Feature selection", "started on",date(),"\n"))
  feature.long <- NULL
  for (i in 1:length(res)) {
    # i<-1
    feature_i <- res[[i]]$feature.importance
    feature_i$Fold <- names(res)[i]
    feature.long <- rbind(feature.long, feature_i)
  }
  feature.final.importance <- feature.long %>%
    group_by(Feature) %>%
    dplyr::summarise(FinalImportance = median(Importance, na.rm = TRUE)) %>%
    dplyr::arrange(desc(FinalImportance)) %>%
    as.data.frame()
  feature.final.importance$Feature <- as.character(feature.final.importance$Feature)
  feature.final.importance$Cumsum <- cumsum(feature.final.importance$FinalImportance)
  feature.final.importance.sum <- sum(feature.final.importance$FinalImportance)
  feature.final.importance.cutoff <- feature.final.importance.sum*0.95
  feature.number <- 1:nrow(feature.final.importance) %>% .[feature.final.importance$Cumsum>feature.final.importance.cutoff] %>% .[1]
  feature.final.importance$CheckCutoff <- feature.final.importance$Cumsum > feature.final.importance.cutoff
  final.feature <- feature.final.importance$Feature[1:feature.number] %>% as.character()
  cat(paste("Feature selection", "finished on",date(),"\n"))
  
  # Final model
  cat(paste("Train final model", "started on",date(),"\n"))
  
  final_input_train_rec <- 
    recipe(Outcome ~ ., data = input[,c("Outcome", final.feature)]) %>%
    step_zv(all_predictors())
  final_boost_wf <- workflow() %>%
    add_model(boost_mod) %>% 
    add_recipe(final_input_train_rec)
  final2_fitted_wf <- final_boost_wf %>% fit(data=input)
  cat(paste("Train final model", "finished on",date(),"\n"))
  
  # Final result
  final_res <- list(cv = res,
                    ready.model = final2_fitted_wf,
                    features = feature.final.importance)
  return(final_res)
}

## Machine learning model for regression task based on lightgbm and cross-validation (not nested)
## Only do feature selection,  model parameters are pre-defined by experience, the fastest option
my_gbt_reg.cv.predef <- function(input_df, Sample, Outcome, Feature, treeN = 1000,
                                   tree_depth = 15, min_n = 30,sample_size = 0.6,mtry = round(0.5*length(Feature)),
                                   learn_rate = 0.01,  inside.v = 5){
  ## test
  # input_df   <- inData
  # Sample  <- trainSample$V1
  # Outcome <- "BMI"
  # Feature <- feature$V1
  # treeN = 2000
  # tree_depth = 10
  # min_n = 3
  # sample_size = 0.8
  # mtry = round(0.2*length(Feature))
  # learn_rate = 0.001
  # inside.v = 3
  
  ## test
  require(tidyverse)
  require(tidymodels)
  require(bonsai)
  require(lightgbm)
  
  # clean data
  input <- data.frame(Outcome = input_df[Sample,Outcome], input_df[Sample, Feature], check.names = T)
  input$Outcome<-as.numeric(input$Outcome)
  input <- input[!is.na(input$Outcome),]
  
  #  input$Outcome <- factor(input$Outcome, levels = c(0, 1))
  
  set.seed(42)
  folds <- vfold_cv(input, v = inside.v)
  
  res <- NULL
  for (i in 1:length(folds$splits)) {
    ## outside loop
    #    i <- 1
    cat(paste("Round",i, "started on",date(),"\n"))
    # Testset split
    input_train <- training(folds$splits[[i]])
    input_test  <- testing(folds$splits[[i]])
    
    # recipe
    input_train_rec <- 
      recipe(Outcome ~ ., data = folds$splits[[i]]) %>%
      step_zv(all_predictors())
    
    # Pre-training: construct model
    boost_mod <- boost_tree(trees = treeN,
                            tree_depth = tree_depth, min_n = min_n,
                            # loss_reduction = tune(),                     ## first three: model complexity
                            sample_size = sample_size, 
                            mtry = mtry,         ## randomness
                            learn_rate = learn_rate) %>% 
      set_engine("lightgbm") %>% set_mode("regression")
    boost_wf <- workflow() %>% add_model(boost_mod) %>% add_recipe(input_train_rec)
    final_fitted_wf <- boost_wf %>% fit(data=input_train)
    
    # Independent test and assessment
    final_fit <- final_fitted_wf %>% last_fit(folds$splits[[i]]) 
    final_perf <- collect_metrics(final_fit)$.estimate # Final performance of current inner CV
    cat(paste("Round",i, "performance",final_perf,"\n"))
    
    # Final prediction
    final_prediction<-bind_cols(
      SampleID = rownames(input),
      SampleRole = "Test",
      ActualClass = input$Outcome,
      predict(final_fitted_wf, input)
    )
    final_prediction$SampleRole[folds$splits[[i]]$in_id]<-"Train"
    
    # Get feature importance
    shp <- shapviz::shapviz(extract_fit_engine(final_fitted_wf), X_pred = as.matrix(input_train[,final_fitted_wf$pre$mold$predictors%>%colnames()]))
    shp.df <- shapviz::sv_importance(shp,max_display=Inf)
    importance <- data.frame(Feature = shp.df$data$feature,
                             Importance = shp.df$data$value)
    
    
    res[[i]]<-list(best.model = final_fitted_wf, final.fit = final_fit, 
                   final.performance = final_perf, final_prediction = final_prediction,
                   feature.importance = importance)
    cat(paste("Round",i, "finished on",date(),"\n"))
  }
  
  names(res) <- folds$id
  
  # Select features for final model
  cat(paste("Feature selection", "started on",date(),"\n"))
  feature.long <- NULL
  for (i in 1:length(res)) {
    feature_i <- res[[i]]$feature.importance
    feature_i$Fold <- names(res)[i]
    feature.long <- rbind(feature.long, feature_i)
  }
  feature.final.importance <- feature.long %>%
    group_by(Feature) %>%
    dplyr::summarise(FinalImportance = median(Importance, na.rm = TRUE)) %>%
    dplyr::arrange(desc(FinalImportance)) %>%
    as.data.frame()
  feature.final.importance$Feature <- as.character(feature.final.importance$Feature)
  feature.final.importance$Cumsum <- cumsum(feature.final.importance$FinalImportance)
  feature.final.importance.sum <- sum(feature.final.importance$FinalImportance)
  feature.final.importance.cutoff <- feature.final.importance.sum*0.95
  feature.number <- 1:nrow(feature.final.importance) %>% .[feature.final.importance$Cumsum>feature.final.importance.cutoff] %>% .[1]
  feature.final.importance$CheckCutoff <- feature.final.importance$Cumsum > feature.final.importance.cutoff
  final.feature <- feature.final.importance$Feature[1:feature.number] %>% as.character()
  cat(paste("Feature selection", "finished on",date(),"\n"))
  
  # Final model
  cat(paste("Train final model", "started on",date(),"\n"))
  
  final_input_train_rec <- 
    recipe(Outcome ~ ., data = input[,c("Outcome", final.feature)]) %>%
    step_zv(all_predictors())
  final_boost_wf <- workflow() %>% add_model(boost_mod) %>% add_recipe(final_input_train_rec)
  final2_fitted_wf <- boost_wf %>% fit(data=input)
  cat(paste("Train final model", "finished on",date(),"\n"))
  
  final_res <- list(cv = res,
                    ready.model = final2_fitted_wf,
                    features = feature.final.importance)
  return(final_res)
}


my_gbt_class.cv.predef <- function(input_df, Sample, Outcome, Feature, treeN = 1000,
                                 tree_depth = 15, min_n = 30,sample_size = 0.6,mtry = round(0.5*length(Feature)),
                                 learn_rate = 0.01,  inside.v = 5){
  ## test
  # input_df   <- inData
  # Sample  <- trainSample$V1
  # Outcome <- "BMI"
  # Feature <- feature$V1
  # treeN = 2000
  # tree_depth = 10
  # min_n = 3
  # sample_size = 0.8
  # mtry = round(0.2*length(Feature))
  # learn_rate = 0.001
  # inside.v = 3
  
  ## test
  require(tidyverse)
  require(tidymodels)
  require(bonsai)
  require(lightgbm)
  
  # clean data
  input <- data.frame(Outcome = input_df[Sample,Outcome], input_df[Sample, Feature], check.names = T)
  input$Outcome<-as.numeric(input$Outcome)
  input <- input[!is.na(input$Outcome),]
  input$Outcome <- factor(input$Outcome, levels = c(0, 1))
  
  
  #  input$Outcome <- factor(input$Outcome, levels = c(0, 1))
  
  set.seed(42)
  folds <- vfold_cv(input, v = inside.v)
  
  res <- NULL
  for (i in 1:length(folds$splits)) {
    ## outside loop
    #    i <- 1
    cat(paste("Round",i, "started on",date(),"\n"))
    # Testset split
    input_train <- training(folds$splits[[i]])
    input_test  <- testing(folds$splits[[i]])
    
    # recipe
    input_train_rec <- 
      recipe(Outcome ~ ., data = folds$splits[[i]]) %>%
      step_zv(all_predictors())
    
    # Pre-training: construct model
    boost_mod <- boost_tree(trees = treeN,
                            tree_depth = tree_depth, min_n = min_n,
                            # loss_reduction = tune(),                     ## first three: model complexity
                            sample_size = sample_size, 
                            mtry = mtry,         ## randomness
                            learn_rate = learn_rate) %>% 
      set_engine("lightgbm") %>% set_mode("classification")
    boost_wf <- workflow() %>% add_model(boost_mod) %>% add_recipe(input_train_rec)
    final_fitted_wf <- boost_wf %>% fit(data=input_train)
    
    # Independent test and assessment
    final_fit <- final_fitted_wf %>% last_fit(folds$splits[[i]]) 
    final_perf <- collect_metrics(final_fit)$.estimate # Final performance of current inner CV
    cat(paste("Round",i, "performance",final_perf,"\n"))
    
    # Final prediction
    final_prediction<-bind_cols(
      SampleID = rownames(input),
      SampleRole = "Test",
      ActualClass = input$Outcome,
      predict(final_fitted_wf, input)
    )
    final_prediction$SampleRole[folds$splits[[i]]$in_id]<-"Train"
    
    # Get feature importance
    shp <- shapviz::shapviz(extract_fit_engine(final_fitted_wf), X_pred = as.matrix(input_train[,final_fitted_wf$pre$mold$predictors%>%colnames()]))
    shp.df <- shapviz::sv_importance(shp,max_display=Inf)
    importance <- data.frame(Feature = shp.df$data$feature,
                             Importance = shp.df$data$value)
    
    
    res[[i]]<-list(best.model = final_fitted_wf, final.fit = final_fit, 
                   final.performance = final_perf, final_prediction = final_prediction,
                   feature.importance = importance)
    cat(paste("Round",i, "finished on",date(),"\n"))
  }
  
  names(res) <- folds$id
  
  # Select features for final model
  cat(paste("Feature selection", "started on",date(),"\n"))
  feature.long <- NULL
  for (i in 1:length(res)) {
    feature_i <- res[[i]]$feature.importance
    feature_i$Fold <- names(res)[i]
    feature.long <- rbind(feature.long, feature_i)
  }
  feature.final.importance <- feature.long %>%
    group_by(Feature) %>%
    dplyr::summarise(FinalImportance = median(Importance, na.rm = TRUE)) %>%
    dplyr::arrange(desc(FinalImportance)) %>%
    as.data.frame()
  feature.final.importance$Feature <- as.character(feature.final.importance$Feature)
  feature.final.importance$Cumsum <- cumsum(feature.final.importance$FinalImportance)
  feature.final.importance.sum <- sum(feature.final.importance$FinalImportance)
  feature.final.importance.cutoff <- feature.final.importance.sum*0.95
  feature.number <- 1:nrow(feature.final.importance) %>% .[feature.final.importance$Cumsum>feature.final.importance.cutoff] %>% .[1]
  feature.final.importance$CheckCutoff <- feature.final.importance$Cumsum > feature.final.importance.cutoff
  final.feature <- feature.final.importance$Feature[1:feature.number] %>% as.character()
  cat(paste("Feature selection", "finished on",date(),"\n"))
  
  # Final model
  cat(paste("Train final model", "started on",date(),"\n"))
  
  final_input_train_rec <- 
    recipe(Outcome ~ ., data = input[,c("Outcome", final.feature)]) %>%
    step_zv(all_predictors())
  final_boost_wf <- workflow() %>% add_model(boost_mod) %>% add_recipe(final_input_train_rec)
  final2_fitted_wf <- boost_wf %>% fit(data=input)
  cat(paste("Train final model", "finished on",date(),"\n"))
  
  final_res <- list(cv = res,
                    ready.model = final2_fitted_wf,
                    features = feature.final.importance)
  return(final_res)
}


batch_coda_glmnet<-function(x_input, y_input, covar_input = NULL, batch_size = 100){
  ## test
  # x_input <- apk.test.metab.assoc # data frame of compositional dataset
  # y_input <- input_snv_age # Vector of Y
  # covar_input <- as.data.frame(input_age) # data frame containing all covariates, if no covariate, set to NULL
  # batch_size <- 100
  ##
  
  require(coda4microbiome)
  
  # Set batches
  batch_n <- ceiling(ncol(x_input)/batch_size)
  batch_index <- split(1:ncol(x_input), 
                       ceiling(seq_along(1:ncol(x_input)) / batch_size))
  preSelect_feature <- NULL
  
  # Run batches
  for (i in 1:batch_n) {
    #  i <- 1
    cat(paste(i, "\n"))
    set.seed(42)
    coda_glmnet_res_i <- coda_glmnet(x=x_input[, batch_index[[i]] ], 
                                     y=y_input,
                                     covar = covar_input,
                                     nfolds = 10, alpha = 0.9,lambda = "lambda.min",
                                     showPlots= FALSE)
    preSelect_feature <- c(preSelect_feature, coda_glmnet_res_i$taxa.name)
  }
  
  # Final model
  set.seed(42)
  coda_glmnet_res <- coda_glmnet(x=x_input[, preSelect_feature ], 
                                 y=y_input,
                                 covar = covar_input,
                                 nfolds = 10, alpha = 0.9,lambda = "lambda.min",
                                 showPlots= FALSE)
  return(coda_glmnet_res)
}


## Mediation analysis for linear model
my_lm_mediation<-function(input.inv,input.med, input.dv,  covDf){
  #input.inv<-lld_exp$melatonine
  #input.med<-lld_vsv$`Faecalibacterium cf. prausnitzii KLE1255:1373_1377`
  #input.dv <-lld_ba$CA_dehydro_deconju_ratio
  #covDf<-lld_covar[,c('Gender','Age','BMI','Reads_number','s__Faecalibacterium_prausnitzii')]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=lm(input.dv~.,input.df.rmna[,-3])
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=lm(input.med~., input.df.rmna[,-2]) # input.inv
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=lm(input.dv~.,input.df.rmna)
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
    }
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
  }
  
  return(res_list)
}

## Bidirectional mediation analysis for linear model
lm_bimediation<-function(inVec, indvDf, dvDf1, dvDf2, covDf, covar ){
  ## test block
  # inVec<-as.character(exp_snp_intrin_df[1,])
  # indvDf<-coh.ptDS.dmp
  # dvDf1<-phen.dmp
  # dvDf2<-phen.dmp
  # covDf<-phen.dmp
  covar<-c('Age','Sex')
  # test block
  
  # if(is.na(inVec[4])){
  #   covar<-covar
  # }else{
  #   covar<-c(covar, colnames(covDf)[grep(inVec[4],colnames(covDf))])
  # }
  
  indv<-indvDf[,match(inVec[1], colnames(indvDf))]
  dv1<-dvDf1[,match(inVec[2], colnames(dvDf1))]
  dv2<-dvDf2[,match(inVec[3], colnames(dvDf2))]
  
  dir1_res <- my_lm_mediation(indv, dv1, dv2, covDf[,covar])
  dir2_res <- my_lm_mediation(indv, dv2, dv1, covDf[,covar])
  
  names(dir1_res)<-paste("dir1.",names(dir1_res),sep = "")
  names(dir2_res)<-paste("dir2.",names(dir2_res),sep = "")
  
  MediationDirection<-"none"
  if(!is.na(dir1_res$dir1.Prop.mediated.p) & !is.na(dir2_res$dir2.Prop.mediated.p)){
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "both"}
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p>0.05){MediationDirection <- "indv_dv1_dv2"}
    if( dir1_res$dir1.Prop.mediated.p>0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "indv_dv2_dv1"}
  }
  
  bires<-list(indv=inVec[1], dv1=inVec[2], dv2=inVec[3],MediationDirection = MediationDirection)
  
  res<-c(bires,dir1_res,dir2_res)
  
  return(res)
}

## Mediation analysis for logistic model
my_lr_mediation_1<-function(input.inv, input.dv, input.med, covDf){
  #input.inv<-indv
  #input.med<-dv2
  #input.dv <-dv1
  #covDf<- covDf[,covar]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  med<-input.df.rmna$input.med
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  input.df.rmna$input.med<-med
  #input.df.rmna$input.med <- factor(input.df.rmna$input.med, levels = c(0,1))
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=lm(input.dv~.,input.df.rmna[,-3])
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=glm(input.med~., input.df.rmna[,-2], family = "binomial")
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=lm(input.dv~.,input.df.rmna) 
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
      
    }
    
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
    
  }
  
  return(res_list)
}


my_lr_mediation_2<-function(input.inv, input.dv, input.med, covDf){
  #input.inv<-lld_exp$cereals
  #input.med<-lld_ba$GDCA
  #input.dv <-lld_dsv$`[Eubacterium] hallii DSM 3353:2969_2983`
  #covDf<-lld_covar[,c('Gender','Age','BMI','Reads_number','s__Eubacterium_hallii')]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  dv<-input.df.rmna$input.dv
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  input.df.rmna$input.dv<-dv
  #input.df.rmna$input.dv <- factor(input.df.rmna$input.dv, levels = c(0,1))
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=glm(input.dv~.,input.df.rmna[,-3], family = "binomial")
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=lm(input.med~., input.df.rmna[,-2])
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=glm(input.dv~.,input.df.rmna, family = "binomial") 
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
      
    }
    
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
    
  }
  
  return(res_list)
}


lr_bimediation<-function(inVec, indvDf, dvDf1, dvDf2, covDf, covar ){
  ## test block
  #inVec<-exp_dsv_ba_df[1,]
  #indvDf<-lld_exp
  #dvDf1<-lld_dsv
  #dvDf2<-lld_ba
  #covDf<-lld_covar
  #covar<-c('Gender','Age','BMI','Reads_number')
  # test block
  
  if(is.na(inVec[4])){
    covar<-covar
  }else{
    covar<-c(covar, inVec[4])
  }
  
  indv<-indvDf[,match(inVec[1], colnames(indvDf))]
  dv1<-dvDf1[,match(inVec[2], colnames(dvDf1))]
  dv2<-dvDf2[,match(inVec[3], colnames(dvDf2))]
  
  dir1_res <- my_lr_mediation_1(indv, dv2, dv1, covDf[,covar])
  dir2_res <- my_lr_mediation_2(indv, dv1, dv2, covDf[,covar])
  
  names(dir1_res)<-paste("dir1.",names(dir1_res),sep = "")
  names(dir2_res)<-paste("dir2.",names(dir2_res),sep = "")
  
  MediationDirection<-"none"
  if(!is.na(dir1_res$dir1.Prop.mediated.p) & !is.na(dir2_res$dir2.Prop.mediated.p)){
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "both"}
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p>0.05){MediationDirection <- "indv_dv1_dv2"}
    if( dir1_res$dir1.Prop.mediated.p>0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "indv_dv2_dv1"}
  }
  
  bires<-list(indv=inVec[1], dv1=inVec[2], dv2=inVec[3],MediationDirection = MediationDirection)
  
  res<-c(bires,dir1_res,dir2_res)
  
  return(res)
}

