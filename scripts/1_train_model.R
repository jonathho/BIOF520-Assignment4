library(tidyverse)
library(caret)
library(missForest)

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

# Confirm missing values are removed
sum(is.na(train_data_imputed))  
sum(is.na(test_data_imputed))  
