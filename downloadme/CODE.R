#### ============================================= #### 
#### Geographic variation in Manx shearwater calls #### 
#### ============================================= #### 
##### == SET UP
require(dplyr)
require(tidyr)
require(ggplot2)
require(paletteer)
require(tibble)

require(factoextra)
require(lme4)
require(lmerTest)
require(ggcorrplot)
require(emmeans)
library(car)
library(randomForest)
library(caret)

# Prepare variable loadings & plot labels
acoustic_vars <- c("Length", "Peak.frequency.Mean", "Mean.frequency.Mean", "Median.frequency.Mean", 
                   "Fundamental.frequency.Mean", "Mean.frequency.change.Mean", "Fundamental.frequency.change.Mean", 
                   "Harmonicity.Mean", "Frequency.bandwidth.Mean", "Vibrato.asymmetry.Mean",
                   "Wiener.entropy.Mean", "Within.syllable.gap", "Between.syllable.gap")

other_vars <- c("Syllable.Number", "Time.Step", "Frame.Length", "Maximum.frequency", "Windowing.Method",
                "Dynamic.Range", "Dynamic.Equalization", "Echo.Tail", "Echo.Reduction", "dy",
                "Syllable.repetitions.per.phrase")

freq_change <- c("Peak.frequency.change.Mean", "Mean.frequency.change.Mean", 
                 "Median.frequency.change.Mean", "Fundamental.frequency.change.Mean")
island <- c("a", "b", "c", "d", "e")
island_names <- c("Copeland"="a", "Lundy"="b", "Rum"="c", "Scillies"="d", "Skomer"="e")
outlier_rows <- c(208, 232, 242, 298, 602, 634)
pop_lookup <- tibble(Island = c("Lundy", "Scillies", "Copeland", "Skomer", "Rum"), 
                     PopSize = c(25000, 2122, 6888, 699324, 288894))

##### == LOAD DATA - NMDS & PARAMETERS
# Load NMDS CSVs
nmds_motif <- read.csv("nmds_motif_raw.csv", header = T)

##### ============================= NMDS ============================= #####
# As the NMDS and resulting distance matrix was already calculated in LUSCINIA, this 
# section is primarily for visualizing the data.

# Extract colony grouping and sex of individual from dataframe
nmds_motif <- nmds_motif %>%
  separate(Population, into = c("Island", "Sex"), sep = "_", remove = FALSE) %>%
  relocate(Island, Sex, .after = Population) %>%
  select(-X)

# Adding individual index for plot
nmds_motif_plot <- nmds_motif %>%
  group_by(Island) %>%
  mutate(IndvIndex = factor(as.integer(factor(Individual)))) %>%
  relocate(IndvIndex, .after = Individual) %>%
  ungroup()

# Plot NMDS for the islands separately; colour-blind friendly palette
ggplot(nmds_motif_plot, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = IndvIndex, shape = Sex)) +
  labs(x = "NMDS1", y = "NMDS2", shape = "Sex") +
  stat_ellipse(aes(group = Sex, fill = Sex), geom = "polygon", level = 0.95, alpha = 0.15) + 
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
        panel.background = element_blank(), legend.position = c(.9,.1),
        text = element_text(size=16)) +
  scale_color_manual(values =   c("#9F0162", "#009F81", "#FF5AAF", "#00FCCF", "#8400CD", "#008DF9",
                                  "#00C2F9", "#FFB2FD", "#A40122", "#E20134", "#FF6E3A", "#FFC33B")) +
  guides(color = "none") +
  facet_wrap(~Island, scales = "free", labeller = as_labeller(island_names))

##### ====================== PRINCIPAL COMPONENT ANALYSIS (PCA) ====================== #####
# Load acoustic parameters for PCA
params <- read.csv("PARAMETER.csv", header = T)
head(params, 0)

# Extract island grouping and sex of individual from dataframe
# Also filtering syllable repetitions out for just the motif
params <- params %>%
  drop_na(all_of(acoustic_vars)) %>% select(-"X") %>% 
  filter(Syllable.repetitions.per.phrase == "true") %>%
  select(-all_of(other_vars)) %>%
  rename(Individual = Individual.name, Song = Song.name) %>%
  mutate(Island = sub("^([A-Za-z]).*", "\\1", `Individual`), 
         Sex = sub(".*_([A-Za-z])_.*", "\\1", `Individual`),
         Island = dplyr::recode(Island, C = "Copeland", L = "Lundy", S = "Scillies", K = "Skomer", R = "Rum")) %>%
  relocate(Island, Sex, .after = Individual) %>%
  slice(-outlier_rows)

# Looking at acoustic parameters sample size
table(params$Island, params$Sex)

# Prepare parameters for PCA analysis 
params_pca <- params %>% select(all_of(acoustic_vars))
str(params)

# === RUNNING THE PCA
pca <- prcomp(params_pca,
              center = TRUE, scale. = TRUE)
# Quick visual of PCA to check for success
summary(pca)

##### ======== INTERPRETING PCA RESULTS ========
# Extract variance from the standard deviations
pca_variance <- pca$sdev^2 / sum(pca$sdev^2) * 100

# ==== PCA Scree plot and biplot
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 40), ggtheme = theme_classic())
fviz_pca_biplot(pca, label = "var", habillage = "none", geom = "text", ggtheme = theme_classic())

# View loadings for each variable, for each PC
pca_loadings <- as.data.frame(pca$rotation) %>% rownames_to_column("Variable")

# Kaiser method to retain PCs with eigenvalues higher than 1
eigenvalues <- pca$sdev^2
pca_kaiser <- paste0("PC", seq_along(eigenvalues))[eigenvalues > 1]

# Find total variance explained for components with highest eigenvalues
pca_variance_kaiser <-sum(pca_variance[1:paste0(length(pca_kaiser))])
cat("Total variance explained by retained eigenvalues:", sum(pca_variance[1:paste0(length(pca_kaiser))]))
# Find individual PC variance for each
cat("Variance explained by retained eigenvalues (individual):", pca_variance[1:paste0(length(pca_kaiser))])

# Selecting only the highest loadings for each of the filtered PCs
pca_loadings_top <- pca_loadings %>%
  select(Variable, all_of(pca_kaiser)) %>%   # only the retained PC columns
  pivot_longer(-Variable, names_to = "PC", values_to = "Loading") %>%
  mutate(AbsLoading = abs(Loading)) %>% group_by(PC) %>% 
  slice_max(AbsLoading, n = 6) %>% relocate(PC, .before=Variable)

write.csv(pca_loadings_top, "PCA_top.csv", row.names = FALSE)

##### ======== PLOTTING PCA RESULTS ========
# Extract PCA scores and bind to original dataframe
pca_scores <- as.data.frame(pca$x)
params_complete <- bind_cols(params, pca_scores %>% select(all_of(pca_kaiser)))

# Convert Individual, Island, Sex into factor columns
cols.fac <- c("Individual", "Island", "Sex")
params_complete[cols.fac] <- lapply(params_complete[cols.fac], as.factor)
sapply(params_complete, class)

# Filtering two separate dataframes for males and females
complete_F <- filter(params_complete, Sex == "F")
complete_M <- filter(params_complete, Sex == "M")

# Filtering for outliers
complete_filter <- params_complete %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(pca_kaiser), ~ abs(.) <= 10))

# Plot PCA
ggplot(params_complete, aes(x = PC3, y = PC4)) +
  geom_point(aes(color = Island, shape = Sex)) + 
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
        panel.background = element_blank(),
        text = element_text(size=16)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_ellipse(geom = "polygon", aes(color = Island, group = Island, fill = Island), level = 0.95, alpha = 0.05) +
  labs(x = paste0("PC1 (", round(pca_variance[1], 1), "%)"), 
       y = paste0("PC2 (", round(pca_variance[2], 1), "%)"))


##### ============================= LINEAR MIXED MODELS (LMM) ============================= #####
# PC ~ Sex + Island + (1 | Individual)
# Response variable: PC#; Fixed effects: Sex, Island; Random effect: Individual

# Comparing Sex * Island and Sex + Island for interactive effects (does Sex affect it differently?)
MODEL_isl <- lmer(PC1 ~ Sex + Island + (1 | Individual), data = params_complete, REML = TRUE)
MODEL_int <- lmer(PC1 ~ Sex * Island + (1 | Individual), data = params_complete, REML = TRUE)
# Likelihood ratio test to confirm that more complex model is better
lrtest(MODEL_isl, MODEL_int)

# Result is significantly different: more complex (interactive) model is better. 
# But given the imbalanced sample size, 
# I will continue with the non-interactive linear model.

# Set lists for model results
lmm_models <- list()
lmm_results <- list()

# Establish Scillies as reference island
params_complete$Island <- relevel(params_complete$Island, ref = "Rum")

# For loop to conduct a linear model for each PCA component
for (component in pca_kaiser) {
  lmm_str <- paste0(component, " ~ Sex + Island + (1 | Individual)")
  model <- lmer(as.formula(lmm_str), data = params_complete, REML = FALSE)
  
  # Put models into list
  lmm_models[[component]] <- model
  
  # Put model coefficients into list
  coefs <- as.data.frame(summary(model)$coefficients) %>%
    rownames_to_column("Term") %>%
    mutate(Component = component) %>%
    relocate(Component)
  lmm_results[[component]] <- coefs
}

# Checks outside of the loop
summary(model)
MODEL_PC1 <- lmer(PC1 ~ Sex + Island + (1 | Individual), data = params_complete,  REML = FALSE)
MODEL_PC2 <- lmer(PC2 ~ Sex + Island + (1 | Individual), data = params_complete,  REML = FALSE)
MODEL_PC3 <- lmer(PC3 ~ Sex + Island + (1 | Individual), data = params_complete,  REML = FALSE)
MODEL_PC4 <- lmer(PC4 ~ Sex + Island + (1 | Individual), data = params_complete,  REML = FALSE)

# Summary checks outside of the loop
summary(MODEL_PC1)
summary(MODEL_PC2)
summary(MODEL_PC3)
summary(MODEL_PC4)

# Combine all model results into one data frame
lmm_results_all <- bind_rows(lmm_results)
lmm_results_all

write.csv(lmm_results_all, "lmm_results_all.csv", row.names = FALSE)

##### ======== ESTIMATED MARGINAL MEANS: ALL ISLAND PAIR-WISE TESTS (with sex-effect) ========
# Set options for emmeans results
emm_options(pbkrtest.limit = 3500, lmerTest.limit = 3500)
emm_results <- list()
means_results <- list()

for (component in pca_kaiser) {
  fit <- lmm_models[[component]]
  emm_model <- emmeans(fit, pairwise ~ Island, adjust = "tukey")
  
  # pairwise comparisons
  emm_results[[component]] <- emm_model$contrasts %>%
    as.data.frame() %>% mutate(Component = component) %>% relocate(Component)
  
  # estimated means per island
  means_results[[component]] <- emm_model$emmeans %>%
    as.data.frame() %>% mutate(Component = component) %>% relocate(Component)
}

# Writing results
emm_results_all   <- bind_rows(emm_results)
means_results_all <- bind_rows(means_results)

# Saving results as a .CSV
write.csv(emm_results_all, "emm_results_all.csv", row.names = FALSE)

# Viewing emmeans out of the loop
emmeans(MODEL_PC1, pairwise ~ Island)
emmeans(MODEL_PC2, pairwise ~ Island)
emmeans(MODEL_PC3, pairwise ~ Island)
emmeans(MODEL_PC4, pairwise ~ Island)

##### ========================== RANDOM FOREST (RF) - COLONY ========================== #####
# Set up new rf dataframe
rf_data <- params %>%
  select(Individual, Island, all_of(acoustic_vars)) %>% 
  mutate(Island = factor(Island), Individual = factor(Individual)) %>% na.omit()

# Check individual motif counts for each island
# Individual motif counts
individual_counts <- rf_data %>% count(Individual, Island)
individual_counts

# Splitting data by individual for training
set.seed(127) # For recreating my results
individuals <- unique(rf_data$Individual)

# Split 70% of individuals for training
train_individuals <- sample(individuals, size = round(0.70 * length(individuals)), replace = FALSE)
# Create training and testing datasets for colony (splitting individuals)
train_col <- rf_data %>% filter(Individual %in% train_individuals)
test_col <- rf_data %>% filter(!Individual %in% train_individuals)
rm(train_individuals)

cat("Training individuals:", n_distinct(train_col$Individual),"\n", 
    "Testing individuals:", n_distinct(test_col$Individual), "\n",
    "Training calls:", nrow(train_col), "\n",
    "Testing calls:", nrow(test_col), "\n")

# Check for no overlaps
length(intersect(unique(train_col$Individual), unique(test_col$Individual))) # 0

##### == RANDOM FOREST BY ISLAND (COLONY-LEVEL VARIATION) - TRAINED
set.seed(127)
rf_model <- randomForest(reformulate(acoustic_vars, response = "Island"),
                         data = train_col, ntree = 10000, importance = TRUE)
# View model
print(rf_model)

rf_predictions <- predict(rf_model, newdata = test_col)
rf_accuracy <- mean(rf_predictions == test_col$Island)
cat("Random Forest accuracy:", round(rf_accuracy * 100, 1), "%\n")

# Confusion matrix
confusion <- confusionMatrix(rf_predictions, test_col$Island)
print(confusion)
confusion$table

##### == RANDOM FOREST BY ISLAND (COLONY-LEVEL VARIATION) - UNTRAINED RANDOM CLASSIFYING
set.seed(81)
islands <- levels(train_col$Island)
n_sim <- 10000

random_accuracy <- replicate(n_sim, 
                             mean(sample(islands, size = nrow(test_col), 
                                         replace = TRUE) == test_col$Island))
# Check results of random accuracy
mean(random_accuracy)
sd(random_accuracy)
quantile(random_accuracy, probs = c(0.025, 0.975))

cat("Random Forest accuracy:", round(rf_accuracy * 100, 1), "%\n",
    "Random classification accuracy:", round(mean(random_accuracy) * 100, 1), "%\n",
    "Random classification accuracy SD:", round(sd(random_accuracy) * 100, 1), "%\n",
    "Random 95% range:", round(quantile(random_accuracy, 0.025) * 100, 1),
    "% - ", round(quantile(random_accuracy, 0.975) * 100, 1), "%\n")

##### == VARIABLE IMPORTANCE (COLONY)
importance(rf_model)
importance_results <- data.frame(Variable = rownames(importance(rf_model)),
                                 importance(rf_model))
importance_results

# Plot importance
varImpPlot(rf_model, main = "Random Forest variable importance (Colony)")

##### ========================== RANDOM FOREST (RF) - COLONY ========================== #####
# Check class balance (call counts per individual)
table(rf_data$Individual)

# Create training and testing datasets for colony (splitting calls)
train_indv <- rf_data %>% group_by(Individual) %>% slice_sample(prop = 0.70) %>% ungroup()
test_indv <- rf_data %>% anti_join(train_indv, by = names(rf_data))

cat("Training individuals:", n_distinct(train_indv$Individual),"\n", 
    "Testing individuals:", n_distinct(test_indv$Individual), "\n",
    "Training calls:", nrow(train_indv), "\n",
    "Testing calls:", nrow(test_indv), "\n")

# Calls per individual, training vs testing
train_indv %>% count(Individual)
test_indv %>% count(Individual)

##### == RANDOM FOREST BY ISLAND (INDIVIDUAL-LEVEL VARIATION) - TRAINED
set.seed(127)
rf_model_indv <- randomForest(reformulate(acoustic_vars, response = "Individual"),
                              data = train_indv, ntree = 10000, importance = TRUE)

# View model
print(rf_model_indv)

rf_predictions_indv <- predict(rf_model_indv, newdata = test_indv)
rf_accuracy_indv <- mean(rf_predictions_indv == test_indv$Individual)
cat("Random Forest accuracy:", round(rf_accuracy_indv * 100, 1), "%\n")

# Confusion matrix
confusion_indv <- confusionMatrix(rf_predictions_indv, test_indv$Individual)
confusion_indv$table

##### == RANDOM FOREST BY ISLAND (INDIVIDUAL-LEVEL VARIATION) - UNTRAINED RANDOM CLASSIFYING
set.seed(18)
individuals_levels <- levels(train_indv$Individual)

random_accuracy_indv <- replicate(n_sim, 
                                  mean(sample(individuals_levels, size = nrow(test_indv), 
                                              replace = TRUE) == test_indv$Individual))
# Check results of random accuracy
mean(random_accuracy_indv)
sd(random_accuracy_indv)
quantile(random_accuracy_indv, probs = c(0.025, 0.975))

cat("Random Forest accuracy:", round(rf_accuracy_indv * 100, 1), "%\n",
    "Random classification accuracy:", round(mean(random_accuracy_indv) * 100, 1), "%\n",
    "Random classification accuracy SD:", round(sd(random_accuracy_indv) * 100, 1), "%\n",
    "Random 95% range:", round(quantile(random_accuracy_indv, 0.025) * 100, 1),
    "% - ", round(quantile(random_accuracy_indv, 0.975) * 100, 1), "%\n")

##### == VARIABLE IMPORTANCE (INDIVIDUAL)
importance(rf_model_indv)
importance_results <- data.frame(Variable = rownames(importance(rf_model_indv)),
                                 importance(rf_model_indv))
importance_results

# Plot importance
varImpPlot(rf_model_indv, main = "Random Forest variable importance (Individual)")

##### ========================== REVIEW RF RESULTS ========================== #####
cat("Colony-level RF classification:", "\n",
    "Random Forest accuracy:", round(rf_accuracy * 100, 1), "%\n",
    "Random classification accuracy:", round(mean(random_accuracy) * 100, 1), "%\n",
    "Random classification accuracy SD:", round(sd(random_accuracy) * 100, 1), "%\n",
    "Random 95% range:", round(quantile(random_accuracy, 0.025) * 100, 1),
    "% - ", round(quantile(random_accuracy, 0.975) * 100, 1), "%\n", "\n",
    "Individual-level RF classification:", "\n",
    "Random Forest accuracy:", round(rf_accuracy_indv * 100, 1), "%\n",
    "Random classification accuracy:", round(mean(random_accuracy_indv) * 100, 1), "%\n",
    "Random classification accuracy SD:", round(sd(random_accuracy_indv) * 100, 1), "%\n",
    "Random 95% range:", round(quantile(random_accuracy_indv, 0.025) * 100, 1),
    "% - ", round(quantile(random_accuracy_indv, 0.975) * 100, 1), "%\n")

# If you want to rerun:
rm(rf_model)
rm(rf_model_indv)

