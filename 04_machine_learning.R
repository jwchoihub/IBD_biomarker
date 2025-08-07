# 04_machine_learning.R
# Recursive Feature Elimination (HC vs IBD; can swap to other comparisons)

library(randomForest)
library(pROC)
library(caret)
library(vegan)
library(dplyr)

# 1) Load pre-processed CLR table + metadata
data   <- read.csv("CLR_data_genus_HCIBD.csv", row.names=1)
labels <- factor(data$Health_status, levels=c("HC","IBD"))

# 2) Prepare X, y
X <- data[, grep("^Genus", names(data))]
y <- labels

# 3) RFE setup
set.seed(1234) # Call set.seed() at the start and before each stochastic step 
ctrl    <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10, repeats=25)
subsets <- seq_len(ncol(X))

# 4) Run RFE
rfe_res <- rfe(X, y, sizes=subsets, rfeControl=ctrl)
selected <- predictors(rfe_res)

# 5) Evaluate on training & test splits
train_idx <- createDataPartition(y, p=0.8, list=FALSE)
X_train   <- X[train_idx, selected]
y_train   <- y[train_idx]
X_test    <- X[-train_idx, selected]
y_test    <- y[-train_idx]

# Train ROC
train_pred <- predict(rfe_res$fit, X_train, type="prob")[,2]
roc_train  <- roc(y_train, train_pred, levels=c("HC","IBD"), percent=TRUE)
print(auc(roc_train))

# Test ROC
test_pred  <- predict(rfe_res$fit, X_test, type="prob")[,2]
roc_test   <- roc(y_test, test_pred, levels=c("HC","IBD"), percent=TRUE)
print(auc(roc_test))
