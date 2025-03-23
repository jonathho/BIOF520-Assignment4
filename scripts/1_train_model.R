library(tidyverse)
library(caret)
library(missForest)
library(survival)
library(survminer)

# load data
# remove NAs in Recurrence since they would not be useful anyways
uromol = readRDS("./data/UROMOL_TaLG.teachingcohort.rds")
uromol_expr = uromol |>
  pull(exprs) |>
    as.data.frame()
uromol_df = uromol |>
  select(-exprs, -UROMOL.ID) |>
  cbind(uromol_expr)
uromol_clean = uromol_df[!is.na(uromol_df$Recurrence), ]

# factor chr variables
char_cols <- names(uromol_clean)[sapply(uromol_clean, is.character)]
uromol_clean[char_cols] <- lapply(uromol_clean[char_cols], as.factor)

# split train-test data
set.seed(24)
train_index <- createDataPartition(uromol_clean$Recurrence, p = 0.8, list = FALSE)
train_data <- uromol_clean[train_index, ]
test_data  <- uromol_clean[-train_index, ]

# impute data
train_data_imputed <- missForest(train_data)$ximp
test_data_imputed <- missForest(test_data)$ximp

# confirm missing values are removed
sum(is.na(train_data_imputed))
sum(is.na(test_data_imputed))

# scale numerical variables
# num_cols <- sapply(train_data_imputed, is.numeric)
# train_data_imputed[, num_cols] <- scale(train_data_imputed[, num_cols])

# check and remove variables with little variance
nzv <- nearZeroVar(train_data_imputed)
train_data_filtered <- train_data_imputed[, -nzv]

# remove collinear variables with a 0.8 cutoff
num_cols <- sapply(train_data_filtered, is.numeric)
filtered_data <- train_data_filtered[, num_cols]
filtered_data <- filtered_data[, -findCorrelation(cor(filtered_data), cutoff = 0.8)]
# write.csv(filtered_data, "./data/train_filtered_data.csv")

# use univariate Cox regression to keep top features
univ_cox <- apply(filtered_data[, -which(names(filtered_data) == "Recurrence")], 2, function(x) {
  coxph(Surv(RFS_time, Recurrence) ~ x, data = filtered_data)$coefficients})

# select top 30 variables with the strongest effect
top_vars <- names(sort(abs(univ_cox), decreasing = TRUE)[1:30])  # Keep top 30
filtered_data <- filtered_data[, c(top_vars, "Recurrence")]
# top 30 variables already (unsurprisingly) include RFS_time

# fit cox model
cox_model <- coxph(Surv(RFS_time, Recurrence) ~ ., data = filtered_data)
summary(cox_model)

# get risk scores
filtered_data$risk_score <- predict(cox_model, type = "lp")  # Linear predictor risk score

# split into high and low risk groups
filtered_data$risk_group <- ifelse(filtered_data$risk_score > median(filtered_data$risk_score), "High Risk", "Low Risk")

# plot kaplan-meier curves
ggsurvplot(survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = filtered_data), 
           pval = TRUE, conf.int = TRUE, risk.table = TRUE, 
           legend.title = "Risk Group")
