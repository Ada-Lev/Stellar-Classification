################################################################################
# LIBRARIES
################################################################################
library(nnet) # logistic regression with multiple classes
library(ggplot2)
library(gridExtra)
library(grid)
library(plotly)
library(dplyr)
library(forcats)
library(tidyverse)
################################################################################
# READING THE DATA
################################################################################
data <- read.csv("data/star_classification.csv", stringsAsFactors = TRUE)
################################################################################
# FEATURE SELECTION
################################################################################
# Let us visualise the distribution of attributes

# We need to select meaningful attributes for the analysis using domain 
# knowledge and exploratory analysis. 

# 1. obj_ID is the object identifier, the unique value that identifies the 
# object in the image catalog used by the CAS, which should not be used in
# further analysis because it is not a spectral characteristic

# 2. alpha and delta specify the right ascension and declination angles, 
# respectively. These angles specify the position of the object in the sky,
# making the use of the feature irrelevant for the model

# 3. u, r, g, i, redshift, and z are spectral characteristics. They should be 
# included
# into the model.

# 4. run_ID, rerun_ID, cam_col, fiber_ID, plate, MJD, and field_ID are features 
# related to the the
# measurement setup, which makes them irrelevant for stellar classification
# based on spectroscopic data.

# 5. spec_obj_ID feature is related to spectroscopic observations. We will use
# it in our model.


# Now, let us see the correlation between features
features <- data.frame(u = data$u, r = data$g, g = data$g, i = data$i, 
                       z = data$z, spec_obj_ID = data$spec_obj_ID, 
                       redshift = data$redshift)
target <- data$class

# Compute the correlation matrix
cor_mat <- cor(features)
print(cor_mat)

# We observe that r, g, u, and z are strongly correlated. 
features <- data.frame(r = data$r, g = data$g, u = data$u, z = data$z, 
                       spec_obj_ID = data$spec_obj_ID, 
                       redshift = data$redshift, i = data$i)

# Write that the model performance is maximum if all spectral features are 
# included. Show this #TODO
################################################################################
# DATA PREPROCESSING
################################################################################
features <- scale(features)
data_prep <- as.data.frame(cbind(target, features))
n <- dim(data_prep)[1]
# Data division
set.seed(12345)
tr_ind <- sample(1:n, floor(0.8*n))
tr <- data_prep[tr_ind, ]
te <- data_prep[-tr_ind, ]
# logistic regression with multinom
################################################################################
# DATA VISUALISATION
################################################################################

# Let us visualise the distribution of data points as spec_obj_ID versus red_shift
# colored by class

plot <- ggplot(data = te) + theme_bw() +
  geom_point(mapping = aes(x = redshift, y = spec_obj_ID, 
                           color = as.factor(target))) +
  scale_color_manual(values = c("orange", "#d50000", "#4527a0"), 
                     name = NULL,
                     labels = c("Galaxies", "Quasars", "Stars")) + 
  labs(title = "Clusters of galaxies, quasars, and stars\nin the redshift-spec_obj_ID coordinate system") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 10, face = "plain")) 

print(plot)

# We can observe three distinct clusters that should be easily identified 
# using logistic regression. We see that the red shift value that corresponds to 
# quasars is the highest because these are the most distant objects we can
# observe. Stars have the lowest red shift value. The value of the red shift of 
# galaxies is in between.
################################################################################
# LOGISTIC REGRESSION
################################################################################
logreg <- multinom(formula = target ~ ., data = tr)

pred_tr <- predict(logreg, tr, type = "class")

MCR_tr <- mean(pred_tr != tr$target)
print("MCR_tr is"); print(MCR_tr)

pred_te <- predict(logreg, te, type = "class")
MCR_te <- mean(pred_te != te$target)

print("MCR_te is"); print(MCR_te)

CM_te <- table(te$target, pred_te)
print("Confusion matrix for test data is"); print(CM_te)
# 1 is GALAXY, 2 is QUASAR, 3 is STAR
# We see that stars are easily classified. The confusion can arise between 
# the identification of galaxies and quasars. TODO
################################################################################
# CLASS DISTRIBUTION
################################################################################
# Visualise the distribution of classes
classes = c("GALAXY", "STAR", "QSO")
counts = numeric(3)
for (i in 1:3){
  counts[i] = length(which(data$class == classes[i]))
}

df_classes <- data.frame(classes = classes, counts = counts)

plot <-df_classes %>%
  mutate(classes = fct_reorder(classes, counts)) %>%
  ggplot(data = df_classes, mapping = aes(x = counts, y = reorder(classes, counts))) + 
  theme_bw() +
  geom_bar(stat = "identity", fill = "darkblue", alpha = .6, width = .4,
           color = "black") +
  labs(title = "Distribution of Classes", x = "", y = "") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        title = element_text(size = 12, face = 'bold')) +
  scale_x_continuous(expand = c(0, 10), limits = c(0, 70000)) +
  geom_text(
    aes(1000, y = classes, label = c("Galaxies", "Stars", "Quasars")),
    hjust = 0,
    nudge_x = 0,
    color = "white",
    family = "Calibri",
    size = 7
  ) +
  geom_text(
    aes(counts + 1000, y = classes, label = counts),
    hjust = 0,
    nudge_x = 0.3,
    colour = "black",
    family = "Calibri",
    size = 6
  )
print(plot)

 