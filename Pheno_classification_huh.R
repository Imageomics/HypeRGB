################################################################################
#                             Script README!
# Classification of leaf phenological stages using PLS-DA and LDA
################################################################################
# Packages
################################################################################
library("pls")
library("caret")
library("dplyr")
library("plsVarSel")
library("randomForest")
library("MASS")


################################################################################
# Read data
################################################################################
huh_spectra     = read_csv("FloraPalooza2026/Data/DMWhiteHUHspec1_sp25leaf560_resampled_5nm_updated.csv")

#######################################
# Process data
#######################################
# Convert "day of collection" to months
huh_spectra$month = as.numeric(format(
  as.Date(huh_spectra$doyOfCollection - 1, origin = "2025-01-01"),
  "%m"
))


# Define phenological stages
# Young leaf: Months 3-5; Mature leaf: Months 6-9; Old/Senescent leaf: Months 10-11
# Create "Phenophase" column to correspond with the Week column
huh_spectra$Phenophase  = NA
huh_spectra = huh_spectra %>%
  dplyr::mutate(
    Phenophase = case_when(
      month %in% 3:5   ~ "Young",
      month %in% 6:9   ~ "Mature",
      month %in% 10:11 ~ "Old/Senescent",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Phenophase))

write.csv(huh_spectra, "FloraPalooza2026/Data/DMWhiteHUHspec1_sp25leaf560_resampled_5nm_updated_leafDev.csv", row.names = F)


# Arrange the columns
huh_spectra = huh_spectra %>%
  dplyr::select(
    filename, simpleFilename, projectID, opticalSetupDescription, lightSourceType, measurementAreaDiameter, 
    measurementAreaDiameter, backgroundMaterialDescription, hasLowReflectanceBackground, 
    whiteCalibrationMaterialDescription, measurementQuant, collectionCode, eventDate, backgroundClass, 
    specimenIdentifier, targetTissueClass, targetTissueNumber, tissueDevelopmentalStage, measurementFlags, 
    tissueNotes, hasGlue, measurementIndex, Class, Order, Family, Genus, scientificName, collector, growthForm, 
    doyOfCollection, month, Phenophase, Age, leafKg_m2, leafThicknessMm, greenIndex, file_base, 
    integration1, integration2, integration3, scan.coadds1, scan.coadds2, scan.coadds3, scan.time, 
    temp1, temp2, temp3, battery, error, gpstime, memory.slot, inclinometer.x.offset, inclinometer.y.offset, 
    everything()
  )

# Split data by BranchID into training and testing
train_idx = createDataPartition(huh_spectra$Phenophase, p = 0.8, list = FALSE)

data_train = huh_spectra[train_idx, ]
data_test  = huh_spectra[-train_idx, ]

# Variables for prediction
predictor_columns = 54:ncol(huh_spectra)

################################################################################
# Model tuning: 
## To select the optimal number of components.
################################################################################

# Ensure Phenophase column is factor
data_train$Phenophase = as.factor(data_train$Phenophase)
data_test$Phenophase  = as.factor(data_test$Phenophase)

# Arrange phenophase column
data_train$Phenophase = factor(data_train$Phenophase, levels = c("Young", "Mature", "Old/Senescent"))
data_test$Phenophase  = factor(data_test$Phenophase, levels = c("Young", "Mature", "Old/Senescent"))

# Define Cross-validation
ctrl = trainControl(method = "repeatedcv", 
                    number = 10, repeats = 100, # 1000 reps
                    classProbs = FALSE, summaryFunction = defaultSummary,
                    sampling = NULL, savePredictions = "final")


################################################################################
# Because of sampling imbalance across phenological stages,
## We downsampled to reduce oversampling of the most sampled phenophase
## We upsampled to increase the sampling of the least sampled phenophase
## In both scenerio, we recorded equal number of sampling.
################################################################################
# Sanity check for number of sampling
table(huh_spectra$Phenophase)
#Young = 96; Mature = 1340; Old/Senescent = 33

#####################################################
# Train PLS-DA model 
#####################################################
# STEP 1: Downsample majority classes first to Choose Optimal ncomp

cal_down = downSample(
  x = data_train [, predictor_columns],
  y = data_train$Phenophase,
  yname = "Phenophase"
)

# Define grid of PLS components:
#max_comp = 15   
grid_down = expand.grid(ncomp = 1:15)

#####################################################
# Train PLS-DA model using downsampled data
#####################################################
pls_down = train(Phenophase ~ .,
                 data = cal_down, method = "pls", 
                 tuneGrid = grid_down, metric = "Kappa",
                 trControl = ctrl, preProcess = c("center", "scale"))


# Extract optimal number of components:
best_n = pls_down$bestTune$ncomp # This becomes the maximum no. of components for upsampling


# STEP 2: Upsample minority classes: (Restricted to ≤ best_n)
cal_up = upSample(
  x = data_train [, predictor_columns],
  y = data_train$Phenophase,
  yname = "Phenophase"
)

# Restrict tuning grid:
grid_up = expand.grid(ncomp = 1:best_n)


#####################################################
# Train final cross-validated PLS-DA model
#####################################################
pls_up = train(Phenophase ~ .,
               data = cal_up, method = "pls",
               tuneGrid = grid_up, metric = "Kappa",
               trControl = ctrl, preProcess = c("center", "scale"))


# Best components after upsampling:
final_n = pls_up$bestTune$ncomp

#####################################################
# Apply pls_up Model to Validation dataset
#####################################################
pred_val = predict(pls_up, newdata = data_test)

#############################
# Confusion matrix:
#############################
cm_pls = confusionMatrix(pred_val, data_test$Phenophase)

# Extract Performance Metrics
accuracy_pls = cm_pls$overall["Accuracy"]
kappa_pls    = cm_pls$overall["Kappa"]

# Save
saveRDS(cm_pls, "FloraPalooza2026/Data/PLSDA_confusion_matrix_huh_herbarium_Phenophase.rds")

###############################################################
# Extract coefficients for the optimal number of components
###############################################################
pls_coef_array = coef(pls_up$finalModel, ncomp = final_n)

# Convert to a data frame
pls_coef_df = as.data.frame(pls_coef_array[, , 1])

# Add wavelength
pls_coef_df$Wavelength = as.numeric(gsub("^X", "", rownames(pls_coef_df)))

# Reorder columns
pls_coef = pls_coef_df %>%
  dplyr::select(Wavelength, Young, Mature, `Old/Senescent`)

# Save
saveRDS(pls_coef, "FloraPalooza2026/Data/PLSDA_coefficients.rds")


##############################################
# PLS-DA Plot
# Calculate phenopase percentage accuracy
##############################################
cm_table = as.data.frame(cm_pls$table)
colnames(cm_table) = c("Prediction", "Reference", "Freq")

diag_percent = cm_table %>%
  group_by(Reference) %>%
  mutate(Total = sum(Freq)) %>%
  ungroup() %>%
  mutate(Percent = 100 * Freq / Total)


# PLSDA Plot
plsda_plot = ggplot(diag_percent, aes(x = Prediction, y = Reference, fill = Percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(Percent, 0))),
            size = 4) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal(base_size = 11) +
  labs(x = "", # "Predicted Identity",
       y = "True Identity", 
       title = "HU Herbarium Spectra - White et al. 2025",
       subtitle = "PLS-DA") +
  #ggtitle("PLS-DA Classification") +
  theme(
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 9),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.3),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#ggsave("Figs/Phenophase_Herbarium_Black_PLSDA_plot.pdf", plot = plsda_plot, height = 6.5, width = 7.5, dpi = 600)

####################################
# VIP Computation and Plotting
####################################
pls_final = pls_up$finalModel

# Compute VIP using optimal components
pls_ncomp_opt = pls_up$bestTune$ncomp

pls_vip_scores = plsVarSel::VIP(pls_final, pls_ncomp_opt)

pls_vip_df = data.frame(
  Wavelength = names(pls_vip_scores),
  VIP = pls_vip_scores
)

pls_vip_df = pls_vip_df[order(-pls_vip_df$VIP), ]
#head(pls_vip_df)

# Plot VIP
pls_vip_df$Wavelength = as.numeric(gsub("^X", "", rownames(pls_vip_df)))
#head(pls_vip_df)

plsda_vip_plot = ggplot(pls_vip_df, aes(x = Wavelength, y = VIP)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(x = "Wavelength (nm)",
       y = "VIP Score",
       title = "PLS-DA Variable Importance (VIP)",
       subtitle = paste("Components =", pls_ncomp_opt)) +
  theme(panel.background = element_rect(fill = "#EEF4FB", color = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 0.3),
        panel.grid.major = element_line(color = "grey75", size = 0.3),
        panel.grid.minor = element_line(color = "grey75", size = 0.2),
        plot.title = element_text(hjust = 0.5),
        #legend.position = "none",
        #legend.text = element_text(size = 8, face = "italic"),
        legend.title = element_text(size = 10, face = "bold"))

ggsave("FloraPalooza2026/Figs/Phenophase_herbarium_VIP_PLSDA_plot.pdf", plot = plsda_vip_plot, height = 6, width = 8.5, dpi = 600)
ggsave("FloraPalooza2026/Figs/Phenophase_herbarium_VIP_PLSDA_plot.png", plot = plsda_vip_plot, height = 6, width = 8.5, dpi = 600)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
################################################################################
# LINEAR DISCRIMINANT ANALYSIS (LDA)
################################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
###############################################
# Train LDA model
# STEP 1 — Downsampling 
## (No tuning grid needed for basic LDA)
###############################################

lda_down = train(Phenophase ~ .,
                 data = cal_down,
                 method = "lda", metric = "Kappa",
                 trControl = ctrl, preProcess = c("center", "scale")
)

#########################################
# STEP 2 — Upsampling (Final model)
#########################################

lda_up = train(Phenophase ~ .,
               data = cal_up,
               method = "lda", metric = "Kappa",
               trControl = ctrl, preProcess = c("center", "scale")
)

#################
# Validation
#################

pred_lda = predict(lda_up, newdata = data_test)
cm_lda   = confusionMatrix(pred_lda, data_test$Phenophase)

accuracy_lda = cm_lda$overall["Accuracy"]
kappa_lda    = cm_lda$overall["Kappa"]

# Save
saveRDS(cm_lda, "FloraPalooza2026/Data/LDA_confusion_matrix_HUH_herbarium_Phenophase.rds")

###############################################################
# Extract coefficients 
###############################################################
lda_coef_matrix = lda_up$finalModel$scaling

lda_coef <- as.data.frame(lda_coef_matrix) %>%
  tibble::rownames_to_column("Predictor") %>%
  mutate(
    Wavelength = as.numeric(gsub("^X", "", Predictor))
  ) %>%
  dplyr::select(Wavelength, everything(), -Predictor)

# Save
saveRDS(lda_coef, "FloraPalooza2026/Data/LDA_coefficients.rds")


#################
# LDA Plot
#################
# Convert Confusion matrix into a table
cm_lda_table = as.data.frame(cm_lda$table)
colnames(cm_lda_table) = c("Prediction", "Reference", "Freq")

## Calculate leaf development percentage accuracy
cm_lda_percent = cm_lda_table %>%
  group_by(Reference) %>%
  mutate(Total = sum(Freq)) %>%
  ungroup() %>%
  mutate(Percent = 100 * Freq / Total)


# LDA Plot
lda_plot = ggplot(cm_lda_percent, aes(x = Prediction, y = Reference, fill = Percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(Percent, 0))),
            size = 4) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal(base_size = 11) +
  labs(x = "Predicted Identity",
       y = "True Identity",
       subtitle = "LDA") +
  #ggtitle("Linear Discriminant Analysis Classification") +
  theme(
    plot.subtitle = element_text(size = 9),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.3),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#ggsave("Figs/Pheno_Class/Phenophase_herbarium_Black_LDA_plot.pdf", plot = lda_plot, height = 6.5, width = 7.5, dpi = 600)


#####################################################
# LDA Variable Importance Computation and Plotting
#####################################################
lda_final = lda_up$finalModel

# Extract scaling (coefficients of discriminant functions)
lda_coef = lda_final$scaling

# For multi-class, combine importance across LD axes
lda_importance = apply(abs(lda_coef), 1, sum)

lda_imp_df = data.frame(
  Wavelength = as.numeric(gsub("^X", "", names(lda_importance))),
  Importance = lda_importance
)

lda_imp_df = lda_imp_df[order(lda_imp_df$Wavelength), ]
rownames(lda_imp_df) = NULL

# Plot VI
lda_vi_plot = ggplot(lda_imp_df, aes(x = Wavelength, y = Importance)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_minimal(base_size = 12) +
  labs(x = "Wavelength (nm)",
       y = "Mean Decrease Accuracy",
       title = "LDA Variable Importance") +
  theme(panel.background = element_rect(fill = "#EEF4FB", color = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 0.3),
        panel.grid.major = element_line(color = "grey75", size = 0.3),
        panel.grid.minor = element_line(color = "grey75", size = 0.2),
        plot.title = element_text(hjust = 0.5),
        #legend.position = "none",
        #legend.text = element_text(size = 8, face = "italic"),
        legend.title = element_text(size = 10, face = "bold"))

ggsave("FloraPalooza2026/Figs/Phenophase_Herbarium_VI_LDA_plot.pdf", plot = lda_vi_plot, height = 6, width = 8.5, dpi = 600)
ggsave("FloraPalooza2026/Figs/Phenophase_Herbarium_VI_LDA_plot.png", plot = lda_vi_plot, height = 6, width = 8.5, dpi = 600)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
################################################################################
# RANDOM FOREST
################################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
# #############################
# # Define tuning grid
# #############################
# 
# rf_grid = expand.grid(
#   mtry = floor(sqrt(length(predictor_columns))) * c(0.5, 1, 1.5)
# )
# 
# #############################
# # STEP 1 — Downsampling
# #############################
# 
# rf_down = train(Phenophase ~ .,
#                 data = cal_down, method = "rf",
#                 tuneGrid = rf_grid, metric = "Kappa",
#                 trControl = ctrl, importance = TRUE)
# 
# best_mtry = rf_down$bestTune$mtry
# 
# #################################################
# # STEP 2 — Upsampling (restricted to best mtry)
# #################################################
# 
# rf_up = train(Phenophase ~ .,
#               data = cal_up, method = "rf",
#               tuneGrid = expand.grid(mtry = best_mtry),
#               metric = "Kappa", trControl = ctrl,
#               importance = TRUE)
# 
# #########################
# # Validation
# #########################
# 
# pred_rf = predict(rf_up, newdata = data_test)
# cm_rf   = confusionMatrix(pred_rf, data_test$Phenophase)
# 
# accuracy_rf = cm_rf$overall["Accuracy"]
# kappa_rf    = cm_rf$overall["Kappa"]
# 
# # Save
# saveRDS(cm_rf, "FloraPalooza2026/Data/RF_confusion_matrix_HUH_herbarium_Phenophase.rds")
# 
# #######################
# # RF Plot
# #######################
# # Convert Confusion matrix into a table
# cm_rf_table = as.data.frame(cm_rf$table)
# colnames(cm_rf_table) = c("Prediction", "Reference", "Freq")
# 
# ## Calculate per-species percentage accuracy
# cm_rf_percent = cm_rf_table %>%
#   group_by(Reference) %>%
#   mutate(Total = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Percent = 100 * Freq / Total)
# 
# 
# # RF Plot
# rf_plot = ggplot(cm_rf_percent, aes(x = Prediction, y = Reference, fill = Percent)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = paste0(round(Percent, 0))),
#             size = 4) +
#   scale_fill_gradient(low = "white", high = "darkgreen") +
#   theme_minimal(base_size = 11) +
#   labs(x = "Predicted Identity",
#        y = "True Identity",
#        subtitle = "RF") +
#   #ggtitle("Random Forest Classification") +
#   theme(
#     plot.subtitle = element_text(size = 9),
#     axis.text.x = element_text(angle = 0, hjust = 0.5),
#     panel.border = element_rect(color = "black", fill = NA, size = 0.3),
#     panel.grid = element_blank(),
#     legend.position = "none"
#   )
# 
# #ggsave("FloraPalooza2026/Figs/Phenophase_herbarium_RF_plot.pdf", plot = rf_plot, height = 6.5, width = 7.5, dpi = 600)
# 
# 
# ########################################
# # Importance Computation and Plotting
# ########################################
# rf_final = rf_up$finalModel
# 
# rf_imp = importance(rf_final)
# 
# rf_imp_df = data.frame(
#   Wavelength = as.numeric(gsub("^X", "", rownames(rf_imp))),
#   MeanDecreaseAccuracy = rf_imp[, "MeanDecreaseAccuracy"],
#   MeanDecreaseGini = rf_imp[, "MeanDecreaseGini"]
# )
# 
# rf_imp_df = rf_imp_df[order(rf_imp_df$Wavelength), ]
# rownames(rf_imp_df) = NULL
# 
# 
# # Plot VIP
# rf_vip_plot = ggplot(rf_imp_df, aes(x = Wavelength, y = MeanDecreaseAccuracy)) +
#   geom_line(linewidth = 0.8) +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   theme_minimal(base_size = 12) +
#   labs(x = "Wavelength (nm)",
#        y = "Mean Decrease Accuracy",
#        title = "Random Forest Variable Importance") +
#   theme(panel.background = element_rect(fill = "#EEF4FB", color = NA),
#         panel.border = element_rect(color = "black", fill = NA, size = 0.3),
#         panel.grid.major = element_line(color = "grey75", size = 0.3),
#         panel.grid.minor = element_line(color = "grey75", size = 0.2),
#         plot.title = element_text(hjust = 0.5),
#         #legend.position = "none",
#         #legend.text = element_text(size = 8, face = "italic"),
#         legend.title = element_text(size = 10, face = "bold"))
# 
# ggsave("FloraPalooza2026/Figs/Phenophase_herbarium_VIP_RF_plot.pdf", plot = rf_vip_plot, height = 6, width = 8.5, dpi = 600)
# ggsave("FloraPalooza2026/Figs/Phenophase_herbarium_VIP_RF_plot.png", plot = rf_vip_plot, height = 6, width = 8.5, dpi = 600)
# 

######################################
# Combine plots
######################################
class_models = patchwork::wrap_plots(plsda_plot, lda_plot, #rf_plot,
                                     ncol = 1, nrow = 2, guides = "collect")

ggsave("FloraPalooza2026/Figs/Phenophase_Classification_HUHerbarium_Spectra.pdf", plot = class_models,
       height = 8.5, width = 7, dpi = 600)
ggsave("FloraPalooza2026/Figs/Phenophase_Classification_HUHerbarium_Spectra.png", plot = class_models,
       height = 8.5, width = 7, dpi = 600)


################################################################################
# Save Classification Model results 
################################################################################
save_classification_model_results = function(model, cm, filename_prefix) {
  
  # Cross-validation results
  res = model$resample
  
  summary_stats <- data.frame(
    Accuracy_mean = mean(res$Accuracy),
    Accuracy_sd   = sd(res$Accuracy),
    Kappa_mean    = mean(res$Kappa),
    Kappa_sd      = sd(res$Kappa)
  )
  
  results_list <- list(
    BestTune        = model$bestTune,
    CV_Iterations   = res,
    CV_Summary      = summary_stats,
    Validation_CM   = cm$table,
    Validation_Acc  = cm$overall["Accuracy"],
    Validation_Kappa = cm$overall["Kappa"]
  )
  
  saveRDS(results_list,
          paste0("FloraPalooza2026/Data/", filename_prefix, "_full_results_herbarium_phenophase.rds"))
}

# Apply to each model
save_classification_model_results(pls_up, cm_pls, "PLSDA")
save_classification_model_results(lda_up, cm_lda, "LDA")
#save_classification_model_results(rf_up, cm_rf, "RF")

# Save all models together
all_models_summary = data.frame(
  Model = c("PLS-DA", "LDA"), #, "Random Forest"
  Accuracy = c(
    cm_pls$overall["Accuracy"],
    cm_lda$overall["Accuracy"]
#    cm_rf$overall["Accuracy"]
  ),
  Kappa = c(
    cm_pls$overall["Kappa"],
    cm_lda$overall["Kappa"]
#    cm_rf$overall["Kappa"]
  )
)

write.csv(all_models_summary, 
          "FloraPalooza2026/Data/Pheno_Classification_Model_Validation_Summary_Herbarium_Black_Phenophase.csv",
          row.names = FALSE)


################################################################################
# Compare classification models
################################################################################
# Extract CV Mean and SD
extract_classification_metrics = function(model) {
  
  res = model$resample
  
  acc_mean  = mean(res$Accuracy)
  acc_sd    = sd(res$Accuracy)
  
  kap_mean  = mean(res$Kappa)
  kap_sd    = sd(res$Kappa)
  
  cv_preds = model$pred
  
  # Handle tuning parameter filtering automatically
  best = model$bestTune
  for(n in names(best)){
    cv_preds = cv_preds[cv_preds[[n]] == best[[n]], ]
  }
  
  bal_acc = cv_preds %>%
    group_by(Resample) %>%
    do({
      cm_fold = confusionMatrix(.$pred, .$obs)
      data.frame(BalancedAccuracy = mean(cm_fold$byClass[, "Balanced Accuracy"]),
                 Precision = mean(cm_fold$byClass[, "Pos Pred Value"]))
    })
  
  data.frame(
    Accuracy_Mean = acc_mean,
    Accuracy_SD   = acc_sd,
    Kappa_Mean    = kap_mean,
    Kappa_SD      = kap_sd,
    BalAcc_Mean   = mean(bal_acc$BalancedAccuracy),
    BalAcc_SD     = sd(bal_acc$BalancedAccuracy),
    Prec_Mean     = mean(bal_acc$Precision),
    Prec_SD       = sd(bal_acc$Precision)
  )
}

# Extract CV Mean and SD
perf_pls = extract_classification_metrics(pls_up)
perf_lda = extract_classification_metrics(lda_up)
#perf_rf  = extract_classification_metrics(rf_up)

herbarium_model_comparison = rbind(
  PLSDA = perf_pls,
  LDA   = perf_lda
#  RF    = perf_rf
)

herbarium_model_comparison
write.csv(herbarium_model_comparison, 
          "FloraPalooza2026/Data/Phenophase_Herbarium_classification_models_comparison.csv", 
          row.names = F)
