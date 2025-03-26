# library(languageserver)
# library(httpgd)
library(tidyverse)
library(caret)
library(missForest)
library(survival)
library(survminer)

# load data
# remove NAs in Recurrence since they would not be useful anyways
uromol <- readRDS("./data/UROMOL_TaLG.teachingcohort.rds") |>
  rename(UROMOL_ID = "UROMOL.ID")
uromol_expr <- uromol |>
  pull(exprs) |>
  as.data.frame()
uromol_df <- uromol |>
  select(-exprs, -UROMOL_ID) |>
  cbind(uromol_expr)
uromol_clean <- uromol_df[!is.na(uromol_df$Recurrence), ]

# factor chr variables
char_cols <- names(uromol_clean)[sapply(uromol_clean, is.character)]
uromol_clean[char_cols] <- lapply(uromol_clean[char_cols], as.factor)

# split train-test data
set.seed(24)
train_index <- createDataPartition(uromol_clean$Recurrence, p = 0.8, list = FALSE)
train_data <- uromol_clean[train_index, ]
test_data <- uromol_clean[-train_index, ]

# impute data
train_data_imputed <- missForest(train_data)$ximp
test_data_imputed <- missForest(test_data)$ximp

# confirm missing values are removed
sum(is.na(train_data_imputed))
sum(is.na(test_data_imputed))

# check and remove variables with little variance
nzv <- nearZeroVar(train_data_imputed)
train_data_filtered <- train_data_imputed[, -nzv]

# remove collinear variables with a 0.8 cutoff
num_cols <- sapply(train_data_filtered, is.numeric)
filtered_data <- train_data_filtered[, num_cols]
filtered_data <- filtered_data[, -findCorrelation(cor(filtered_data), cutoff = 0.8)]
# write_tsv(filtered_data, "./data/train_filtered_data.tsv")
# write_tsv(test_data_imputed, "./data/test_data_imputed.tsv")

# load in data to skip steps above
filtered_data <- read_tsv("./data/train_filtered_data.tsv") |>
  as.data.frame()
test_data_imputed <- read_tsv("./data/test_data_imputed.tsv") |>
  as.data.frame()

# use univariate Cox regression to keep top features
univ_cox <- apply(filtered_data[, -which(names(filtered_data) == "Recurrence")], 2, function(x) {
  coxph(Surv(RFS_time, Recurrence) ~ x, data = filtered_data)$coefficients
})

# select variables
top_vars <- names(sort(abs(univ_cox), decreasing = TRUE)[1:9])
filtered_data <- filtered_data[, c("Recurrence", top_vars)]
# top variables already (unsurprisingly) includes RFS_time, no need to add again

# fit cox model
cox_model <- coxph(Surv(RFS_time, Recurrence) ~ ., data = filtered_data)
summary(cox_model)

# get risk scores
filtered_data$risk_score <- predict(cox_model, type = "lp") # Linear predictor risk score

test_data_filtered <- test_data_imputed[, c("Recurrence", top_vars)]
test_data_filtered$risk_score <- predict(cox_model, newdata = test_data_filtered, type = "lp")

# split into high and low risk groups
threshold_riskscore <- median(filtered_data$risk_score)
filtered_data$risk_group <- ifelse(filtered_data$risk_score > threshold_riskscore, "High Risk", "Low Risk")

test_data_filtered$risk_group <- ifelse(test_data_filtered$risk_score > threshold_riskscore, "High Risk", "Low Risk")

# plot kaplan-meier curves
ggsurvplot(survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = filtered_data),
  pval = TRUE, conf.int = TRUE, risk.table = FALSE,
  legend.title = "Risk Group", xlab = "Time (Months)"
)

ggsurvplot(survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = test_data_filtered),
  pval = TRUE, conf.int = TRUE, risk.table = FALSE,
  legend.title = "Risk Group", legend.labs = c("Low Risk", "High Risk"),
  title = "Kaplan-Meier Curves for Test Set"
)


# load in validation data set
knowles <- readRDS("./data/knowles_matched_TaLG_final.rds")
knowles_expr <- knowles |>
  pull(exprs) |>
  as.data.frame()
knowles_df <- knowles |>
  select(-exprs, -knowles_ID) |>
  cbind(knowles_expr)
knowles_clean <- knowles_df[!is.na(knowles_df$Recurrence), ]

# factor chr variables
knowles_char_cols <- names(knowles_clean)[sapply(knowles_clean, is.character)]
knowles_clean[knowles_char_cols] <- lapply(knowles_clean[knowles_char_cols], as.factor)

# impute data
knowles_imputed <- missForest(knowles_clean)$ximp

# predict
knowles_imputed$risk_score <- predict(cox_model, newdata = knowles_imputed, type = "lp")
knowles_imputed$risk_group <- ifelse(knowles_imputed$risk_score > threshold_riskscore, "High Risk", "Low Risk")

ggsurvplot(survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = knowles_imputed),
  pval = TRUE, conf.int = TRUE, risk.table = FALSE,
  legend.title = "Risk Group", legend.labs = c("Low Risk", "High Risk"), palette = c("#619CFF", "#F8766D"),
  title = "Kaplan-Meier Curves for KNOWLES"
)

# create decision matrix (confusion table)
decision_matrix <- table(
  Predicted = knowles_imputed$risk_group,
  Actual = ifelse(knowles_imputed$Recurrence == 1, "Recurred", "Not Recurred")
)

# print the decision matrix
print(decision_matrix)


library(rmda)
# extract baseline survival function at a chosen time point (e.g., 3 years = 36 months)
base_surv <- basehaz(cox_model, centered = FALSE)

# choose the time point for probability estimation
time_point <- 6
baseline_survival_at_time <- exp(-base_surv$hazard[which.min(abs(base_surv$time - time_point))])

# convert risk scores into recurrence probabilities
knowles_imputed$predicted_prob <- 1 - (baseline_survival_at_time^exp(knowles_imputed$risk_score))

# Decision Curve Analysis
dca_results <- decision_curve(Recurrence ~ predicted_prob,
  data = knowles_imputed,
  thresholds = seq(0.01, 0.5, by = 0.01),
  family = "binomial",
  policy = "opt-in"
)

# plot Decision Curve
plot_decision_curve(dca_results,
  curve.names = "Cox Model",
  col = "blue", lwd = 2,
  legend.position = "none",
  standardize = FALSE
)
