rm(list = ls())
gc()
setwd("/data/nas1/mengqingyao_OD/project/Program-179")
if (! dir.exists("./Custome_run_ML")){
  dir.create("./Custome_run_ML")
}
setwd("./Custome_run_ML")
options(stringsAsFactors = FALSE)

library(tidyverse)
library(randomForest)
library(glmnet)
library(xgboost)
library(Boruta)
library(pROC)

# 数据加载 --------------------------------------------------------------------
PCG_v36 <- read_delim("../00_rawdata/PCG_v36.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

expr <- read_csv("../00_rawdata/Seq_verify_dat_GSE126848_count.csv")
expr <- expr %>% column_to_rownames("...1")
expr <- expr[rownames(expr) %in% PCG_v36$gene_name, ]

pdata <- read_csv("../00_rawdata/Seq_verify_group_GSE126848.csv")
expr <- expr[, pdata$sample]
group <- data.frame(barcode = colnames(expr), 
                    sample = factor(pdata$group, 
                                    levels = c("NASH", "NAFLD")))

common <- read_csv("../04_Cad_PPI/common.csv")

expr_ml <- t(na.omit(expr[common$symbol, ]))
x <- as.matrix(expr_ml)
y <- as.numeric(group$sample == "NASH")

# 定义所有方法 ----------------------------------------------------------------
run_lasso <- function() {
  set.seed(0623)
  cvfit <- cv.glmnet(x, y, family="binomial", alpha=1, nfolds=5)
  coef_min <- coef(cvfit, s="lambda.min")
  key_genes_min <- rownames(coef_min)[which(coef_min != 0)][-1]
  return(key_genes_min)
}

run_svm <- function() {
  set.seed(0623)
  y <- as.factor(y)
  
  svmRFE <- function(X, y, k=10, verbose=FALSE) {
    features <- colnames(X)
    ranked_features <- c()
    
    while(length(features) > 0) {
      cv <- svm(x = X[, features], y = y, kernel="linear", cross=k)
      model <- svm(X[, features], y, kernel="linear")
      w <- t(model$coefs) %*% model$SV
      weights <- w * w
      ranking <- sort(weights, index.return=TRUE, decreasing=FALSE)$ix
      ranked_features <- c(features[ranking[1]], ranked_features)
      features <- features[-ranking[1]]
    }
    return(ranked_features)
  }
  
  feature_ranking <- svmRFE(x, y, k=5)
  ctrl <- trainControl(method="repeatedcv", number=5, repeats=3)
  
  results <- data.frame()
  for (i in 1:length(feature_ranking)) {
    current_features <- feature_ranking[1:i]
    x_subset <- x[, current_features, drop = FALSE]
    model <- train(x = x_subset, y = y, method = "svmLinear",
                   trControl = ctrl, preProcess = c("center", "scale"))
    results <- rbind(results, data.frame(NumFeatures = i, Accuracy = model$results$Accuracy[1]))
  }
  
  best_num <- which.max(results$Accuracy)
  best_features <- feature_ranking[1:best_num]
  return(best_features)
}

run_boruta <- function() {
  set.seed(0623)
  boruta_result <- Boruta(x = expr_ml, y = as.factor(group$sample), 
                          doTrace = 2, maxRuns = 100)
  Boruta <- attStats(boruta_result)
  boruta_geneids <- Boruta[Boruta$decision=='Confirmed',] %>% rownames(.)
  return(boruta_geneids)
}

run_rf <- function() {
  set.seed(0623)
  dat <- data.frame(expr_ml, group = as.factor(group$sample))
  
  cv_error <- c()
  ntree_range <- seq(10, 100, by = 2)
  for (ntree in ntree_range) {
    rf_model <- randomForest(group ~ ., data = dat, ntree = ntree, mtry=5, importance = TRUE)
    oob_error <- rf_model$err.rate[ntree, "OOB"]
    cv_error <- c(cv_error, oob_error)
  }
  
  optimal_ntree <- ntree_range[which.min(cv_error)]
  final_rf_model <- randomForest(group ~ ., data = dat, ntree = optimal_ntree, importance = TRUE)
  
  imp_df <- importance(final_rf_model, type = 2) %>%
    as.data.frame() %>%
    arrange(desc(.[,1])) %>%
    head(10) %>%
    rownames_to_column("Gene")
  
  return(imp_df$Gene)
}

run_xgboost <- function() {
  set.seed(0623)
  dtrain <- xgb.DMatrix(data = as.matrix(x), label = y)
  res.xgb <- xgboost(data = dtrain, max_depth = 2, eta = 0.3,
                     objective = "binary:logistic", nrounds = 10, verbose = 0)
  
  importance_matrix <- xgb.importance(feature_names = colnames(x), model = res.xgb)
  importance_matrix$Feature <- gsub(".", "-", importance_matrix$Feature, fixed = TRUE)
  return(importance_matrix$Feature)
}

# 运行所有方法并获取结果 ------------------------------------------------------
all_methods <- list(
  Lasso = run_lasso(),
  SVM = run_svm(),
  Boruta = run_boruta(),
  RF = run_rf(),
  XGBoost = run_xgboost()
)

# 生成所有可能的三方法组合 ----------------------------------------------------
method_names <- names(all_methods)
combinations <- combn(method_names, 3, simplify = FALSE)

# 验证每种组合的交集基因 -----------------------------------------------------
validate_combination <- function(methods) {
  # 获取三种方法的基因列表
  genes_list <- all_methods[methods]
  
  # 计算交集
  common_genes <- Reduce(intersect, genes_list)
  
  if (length(common_genes) < 1) {
    return(list(methods = methods, final_biomarkers_count = 0, biomarkers = NULL))
  }
  
  # ROC和wilcoxon验证
  final_genes <- common_genes
  
  # 加载验证数据
  exprs_train <- read_csv("../00_rawdata/Seq_verify_dat_GSE126848_count.csv")
  exprs_train <- exprs_train %>% column_to_rownames("...1")
  expr_train <- exprs_train[rownames(exprs_train) %in% PCG_v36$gene_name, ] %>% t()
  pdata_train <- read.csv("../00_rawdata/Seq_verify_group_GSE126848.csv")
  group_train <- as.numeric(pdata_train$group == "NASH")
  
  expr_val <- read.csv("../00_rawdata/Array_train_dat.GSE89632.csv", row.names = 1) %>% t()
  pdata_val <- read.csv("../00_rawdata/Array_train_group.GSE89632.csv")
  group_val <- as.numeric(pdata_val$group == "NASH")
  
  expr_train <- scale(expr_train[, final_genes, drop = FALSE])
  expr_val <- expr_val[, final_genes, drop = FALSE]
  
  # ROC分析
  roc_results <- data.frame(
    Gene = colnames(expr_train),
    Train_AUC = sapply(colnames(expr_train), function(g) {
      roc(group_train, expr_train[, g], quiet = TRUE)$auc
    }),
    Val_AUC = sapply(colnames(expr_val), function(g) {
      if (g %in% colnames(expr_val)) {
        roc(group_val, expr_val[, g], quiet = TRUE)$auc
      } else NA
    })
  ) %>% na.omit()
  
  biomarker_candidates <- roc_results %>%
    filter(Train_AUC > 0.7 & Val_AUC > 0.7) %>%
    mutate(AUC_diff = abs(Train_AUC - Val_AUC)) %>%
    arrange(desc(Val_AUC)) %>%
    top_n(15, Val_AUC)
  
  # 如果没有符合条件的候选基因，直接返回
  if (nrow(biomarker_candidates) == 0) {
    return(list(
      methods = methods,
      final_biomarkers_count = 0,
      biomarkers = NULL,
      common_genes_count = length(common_genes),
      common_genes = common_genes
    ))
  }
  
  # Wilcoxon检验
  prepare_expr_df <- function(expr, pdata, genes) {
    df <- as.data.frame(expr[, genes, drop = FALSE])
    df$Group <- pdata$group
    df
  }
  
  expr_train_df <- prepare_expr_df(expr_train, pdata_train, biomarker_candidates$Gene)
  expr_val_df <- prepare_expr_df(expr_val, pdata_val, biomarker_candidates$Gene)
  
  perform_wilcox <- function(gene, train_df, val_df) {
    train_test <- wilcox.test(as.formula(paste(gene, "~ Group")), data = train_df, exact = FALSE)
    val_test <- wilcox.test(as.formula(paste(gene, "~ Group")), data = val_df, exact = FALSE)
    
    train_median_diff <- median(train_df[train_df$Group == "NASH", gene]) -
      median(train_df[train_df$Group == "NAFLD", gene])
    val_median_diff <- median(val_df[val_df$Group == "NASH", gene]) -
      median(val_df[val_df$Group == "NAFLD", gene])
    
    data.frame(
      Gene = gene,
      Train_pvalue = train_test$p.value,
      Val_pvalue = val_test$p.value,
      Train_median_diff = train_median_diff,
      Val_median_diff = val_median_diff,
      Direction_consistent = sign(train_median_diff) == sign(val_median_diff)
    )
  }
  
  wilcox_results <- map_dfr(biomarker_candidates$Gene, perform_wilcox, expr_train_df, expr_val_df)
  
  # 筛选最终生物标志物
  final_biomarkers <- wilcox_results %>%
    filter(Train_pvalue < 0.05 & Val_pvalue < 0.05 & Direction_consistent) %>%
    arrange(Train_pvalue + Val_pvalue)
  
  return(list(
    methods = methods,
    final_biomarkers_count = nrow(final_biomarkers),
    biomarkers = final_biomarkers$Gene,
    biomarker_candidates = biomarker_candidates,
    common_genes_count = length(common_genes),
    common_genes = common_genes
  ))
}

# 遍历所有组合并验证 ---------------------------------------------------------
results <- list()
for (i in seq_along(combinations)) {
  combo <- combinations[[i]]
  cat("Testing combination:", paste(combo, collapse = "+"), "\n")
  res <- validate_combination(combo)
  results[[i]] <- res
  cat("Found", res$final_biomarkers_count, "final biomarkers\n\n")
}

# 筛选出符合条件的组合 -------------------------------------------------------
valid_combinations <- results[sapply(results, function(x) x$final_biomarkers_count >= 2)]



# 输出结果 -------------------------------------------------------------------
if (length(valid_combinations) > 0) {
  cat("\nValid combinations with >= 2 final biomarkers:\n")
  for (i in seq_along(valid_combinations)) {
    combo <- valid_combinations[[i]]
    cat("\nCombination:", paste(combo$methods, collapse = "+"), "\n")
    cat("Number of common genes:", combo$common_genes_count, "\n")
    cat("Number of final biomarkers:", combo$final_biomarkers_count, "\n")
    cat("Final biomarkers:", paste(combo$biomarkers, collapse = ", "), "\n")
    cat("Common genes:", paste(combo$common_genes, collapse = ", "), "\n")
  }
} else {
  cat("\nNo combinations found with >= 2 final biomarkers\n")
}

# 保存结果到文件 -------------------------------------------------------------
saveRDS(list(
  all_methods = all_methods,
  all_results = results,
  valid_combinations = valid_combinations
), "method_combinations_results.rds")

# 简化版Wilcoxon检验函数（不进行ROC分析） -----------------------------------------------
perform_wilcoxon_only <- function(methods, genes_list) {
  # 获取指定方法的基因列表交集
  common_genes <- Reduce(intersect, genes_list[methods])
  
  if (length(common_genes) < 1) {
    return(list(
      methods = methods,
      final_biomarkers_count = 0,
      biomarkers = NULL,
      common_genes_count = 0,
      common_genes = NULL
    ))
  }
  
  # 加载数据
  expr_train <- read_csv("../00_rawdata/Seq_verify_dat_GSE126848_count.csv") %>% 
    column_to_rownames("...1") %>% 
    t()
  pdata_train <- read_csv("../00_rawdata/Seq_verify_group_GSE126848.csv")
  
  expr_val <- read.csv("../00_rawdata/Array_train_dat.GSE89632.csv", row.names = 1) %>% t()
  pdata_val <- read.csv("../00_rawdata/Array_train_group.GSE89632.csv")
  
  # 准备数据
  prepare_expr_df <- function(expr, pdata, genes) {
    df <- as.data.frame(expr[, genes, drop = FALSE])
    df$Group <- pdata$group
    df
  }
  
  expr_train_df <- prepare_expr_df(expr_train, pdata_train, common_genes)
  expr_val_df <- prepare_expr_df(expr_val, pdata_val, common_genes)
  
  # Wilcoxon检验函数
  perform_wilcox <- function(gene) {
    train_test <- wilcox.test(as.formula(paste(gene, "~ Group")), 
                              data = expr_train_df, exact = FALSE)
    val_test <- wilcox.test(as.formula(paste(gene, "~ Group")), 
                            data = expr_val_df, exact = FALSE)
    
    train_median_diff <- median(expr_train_df[expr_train_df$Group == "NASH", gene]) -
      median(expr_train_df[expr_train_df$Group == "NAFLD", gene])
    val_median_diff <- median(expr_val_df[expr_val_df$Group == "NASH", gene]) -
      median(expr_val_df[expr_val_df$Group == "NAFLD", gene])
    
    data.frame(
      Gene = gene,
      Train_pvalue = train_test$p.value,
      Val_pvalue = val_test$p.value,
      Train_median_diff = train_median_diff,
      Val_median_diff = val_median_diff,
      Direction_consistent = sign(train_median_diff) == sign(val_median_diff)
    )
  }
  
  # 执行检验
  wilcox_results <- map_dfr(common_genes, perform_wilcox)
  
  # 筛选显著结果
  final_biomarkers <- wilcox_results %>%
    filter(Train_pvalue < 0.05 & Val_pvalue < 0.05 & Direction_consistent) %>%
    arrange(Train_pvalue + Val_pvalue)
  
  return(list(
    methods = methods,
    final_biomarkers_count = nrow(final_biomarkers),
    biomarkers = final_biomarkers$Gene,
    common_genes_count = length(common_genes),
    common_genes = common_genes,
    wilcox_results = wilcox_results
  ))
}

# 遍历所有组合的简化版函数
run_all_combinations_wilcoxon <- function(all_methods) {
  method_names <- names(all_methods)
  combinations <- combn(method_names, 3, simplify = FALSE)
  
  results <- list()
  for (i in seq_along(combinations)) {
    combo <- combinations[[i]]
    cat("Testing combination:", paste(combo, collapse = "+"), "\n")
    res <- perform_wilcoxon_only(combo, all_methods)
    results[[i]] <- res
    cat("Found", res$final_biomarkers_count, "final biomarkers\n\n")
  }
  
  # 筛选有效组合
  valid_combinations <- results[sapply(results, function(x) x$final_biomarkers_count >= 2)]
  
  # 输出结果
  if (length(valid_combinations) > 0) {
    cat("\nValid combinations with >= 2 final biomarkers:\n")
    for (i in seq_along(valid_combinations)) {
      combo <- valid_combinations[[i]]
      cat("\nCombination:", paste(combo$methods, collapse = "+"), "\n")
      cat("Number of common genes:", combo$common_genes_count, "\n")
      cat("Number of final biomarkers:", combo$final_biomarkers_count, "\n")
      cat("Final biomarkers:", paste(combo$biomarkers, collapse = ", "), "\n")
    }
  } else {
    cat("\nNo combinations found with >= 2 final biomarkers\n")
  }
  
  return(list(
    all_results = results,
    valid_combinations = valid_combinations
  ))
}

# 使用示例
# 1. 运行所有方法获取基因列表
all_methods <- list(
  Lasso = run_lasso(),
  SVM = run_svm(),
  Boruta = run_boruta(),
  RF = run_rf(),
  XGBoost = run_xgboost()
)

# 运行简化版组合分析
wilcoxon_results <- run_all_combinations_wilcoxon(all_methods)

# 3. 保存结果
saveRDS(wilcoxon_results, "wilcoxon_only_results.rds")
