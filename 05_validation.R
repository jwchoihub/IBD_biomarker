# 05_validation.R
# External validation with multiple feature sets (HC vs IBD framework)

library(lightgbm)
library(xgboost)
library(pROC)
library(caret)

# 1) Read feature sets
fs1 <- read.csv("HCIBD_DA_biomarker_genus.csv")
fs2 <- read.csv("HCIBD_ML_biomarker_genus.csv")
fs3 <- read.csv("HCIBD_NW_biomarker_genus.csv")
fs4 <- read.csv("HCIBD_LC_biomarker_genus.csv")

# Helper: split & train/evaluate
split_data <- function(df){
  y <- factor(ifelse(df$Health_status=="HC", 0, 1))
  X <- df[, setdiff(names(df), "Health_status")]
  idx <- createDataPartition(y, p=0.8, list=FALSE)
  list(X_train=X[idx,], y_train=y[idx], X_test=X[-idx,], y_test=y[-idx])
}
process_set <- function(d){
  # LightGBM
  dtrain <- lgb.Dataset(as.matrix(d$X_train), label=as.numeric(as.character(d$y_train)))
  lgb <- lgb.train(
    params=list(objective="binary", metric="auc"),
    data=dtrain, nrounds=150
  )
  pred1 <- predict(lgb, as.matrix(d$X_test))
  # XGBoost
  dmat <- xgb.DMatrix(as.matrix(d$X_train), label=as.numeric(as.character(d$y_train)))
  xgbm <- xgboost(
    params=list(objective="binary:logistic", eval_metric="auc"),
    data=dmat, nrounds=150, verbose=0
  )
  pred2 <- predict(xgbm, as.matrix(d$X_test))
  # Ensemble
  prob <- (pred1 + pred2)/2
  roc_res <- roc(d$y_test, prob, levels=c(0,1), percent=TRUE)
  list(auc=auc(roc_res), roc=roc_res)
}

# 2) Apply to all sets
for(fs in list(fs1, fs2, fs3, fs4)){
  res <- process_set(split_data(fs))
  cat("AUC:", res$auc, "\n")
}
