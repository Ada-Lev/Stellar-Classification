################################################################################
# LIBRARIES
################################################################################
library(nnet) # logistic regression with multiple classes
library(ggplot2)
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
# included. Future scope
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
pred_te <- predict(logreg, te, type = "class")
################################################################################
# PERFORMANCE EVALUATION
################################################################################
CM_tr <- table(pred_tr, tr$target, deparse.level = 0)
CM_te <- table(pred_te, te$target, deparse.level = 0)
rownames(CM_te) <- c("GAL_pr", "QUA_pr", "STA_pr")
colnames(CM_te) <- c("GAL_tr", "QUA_tr", "STA_tr")

cat("The confusion matrix for train data is\n"); print(CM_tr)
cat("The confusion matrix for test data is\n"); print(CM_te)
# 1 is GALAXY, 2 is QUASAR, 3 is STAR
# We see that stars are easily classified. The confusion can arise between 
# the identification of galaxies and quasars. 
# For assessment of the model performance on distinguishing between stars
# and other classes, one can simply use the accuracy metric because of the low number
# of misclassification cases. That is, the number of true negatives (TNs) for
# stars versus other classes is 7. 

# Let us introduce this renewed stars-versus-rest CM
CM_tr_stars_vs_rest <- matrix(c(CM_tr[1, 1] + CM_tr[1, 2] + CM_tr[2, 1] + CM_tr[2, 2], 
                                CM_tr[2, 3] + CM_tr[1, 3],
                                CM_tr[3, 1] + CM_tr[3, 2],
                                CM_tr[3, 3]), byrow = TRUE, ncol = 2)

colnames(CM_tr_stars_vs_rest) <- c("rest_tr", "STA_tr")
rownames(CM_tr_stars_vs_rest) <- c("rest_pr", "STA_pr")

print("For train data")
print(CM_tr_stars_vs_rest)
GT <- factor(c("Rest", "Stars", "Rest", "Stars"))
predicted <- factor(c("Rest", "Rest", "Stars", "Stars"))
values <- c(CM_tr_stars_vs_rest)
df <- data.frame(GT, predicted, values)

# Visualise this CM
plot <- ggplot(data =  df, mapping = aes(x = predicted, y = GT)) +
  geom_tile(aes(fill = values), color = "white") +
  geom_text(data = subset(df, values > 30000),
            aes(label = values), 
            vjust = 1, color = "white", fontface = "bold", family = "arial") +
  geom_text(data = subset(df, values < 30000),
            aes(label = values), 
            vjust = 1, color = "black", fontface = "bold", family = "arial") +
  scale_fill_gradient(low = "#bbdefb", high = "#0d47a1") +
  theme_bw() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev) +
  labs(x = "Predicted", y = "Ground Truth", fill = "# cases",
       title = "Train CM for stars-versus-rest case") +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 10, face = "plain"))

print(plot)

CM_te_stars_vs_rest <- matrix(c(CM_te[1, 1] + CM_te[1, 2] + CM_te[2, 1] + CM_te[2, 2], 
                                CM_te[2, 3] + CM_te[1, 3],
                                CM_te[3, 1] + CM_te[3, 2],
                                CM_te[3, 3]), byrow = TRUE, ncol = 2)

colnames(CM_te_stars_vs_rest) <- c("rest_tr", "STA_tr")
rownames(CM_te_stars_vs_rest) <- c("rest_pr", "STA_pr")

print("For test data")
print(CM_te_stars_vs_rest)
GT <- factor(c("Rest", "Stars", "Rest", "Stars"))
predicted <- factor(c("Rest", "Rest", "Stars", "Stars"))
values <- c(CM_te_stars_vs_rest)
df <- data.frame(GT, predicted, values)

# Visualise this CM
plot <- ggplot(data =  df, mapping = aes(x = predicted, y = GT)) +
  geom_tile(aes(fill = values), color = "white") +
  geom_text(data = subset(df, values > 10000),
            aes(label = values), 
            vjust = 1, color = "white", fontface = "bold", family = "arial") +
  geom_text(data = subset(df, values < 10000),
            aes(label = values), 
            vjust = 1, color = "black", fontface = "bold", family = "arial") +
  scale_fill_gradient(low = "#bbdefb", high = "#0d47a1") +
  theme_bw() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev) +
  labs(x = "Predicted", y = "Ground Truth", fill = "# cases",
       title = "Test CM for stars-versus-rest case") +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 10, face = "plain"))

print(plot)

cat("Train accuracy for distinguishing between stars and other classes\n")
MCR_tr <- (CM_tr_stars_vs_rest[1, 1] + CM_tr_stars_vs_rest[2, 2]) / sum(CM_tr_stars_vs_rest) 
cat(round(MCR_tr*100, 1), "%", sep="")

cat("Test accuracy for distinguishing between stars and other classes\n")
MCR_te <- (CM_te_stars_vs_rest[1, 1] + CM_te_stars_vs_rest[2, 2]) / sum(CM_te_stars_vs_rest) 
cat(round(MCR_te*100, 1), "%", sep="")

# Now, let us consider distinguishing quasars from other data. The confusion
# matrix is as follows:
CM_tr_quasars_vs_rest <- matrix(c(CM_tr[1, 1] + CM_tr[1, 3] + CM_tr[3, 1] + CM_tr[3, 3], 
                                CM_tr[1, 2] + CM_tr[3, 2],
                                CM_tr[2, 1] + CM_tr[2, 3],
                                CM_tr[2, 2]), byrow = TRUE, ncol = 2)

colnames(CM_tr_quasars_vs_rest) <- c("rest_tr", "QUA_tr")
rownames(CM_tr_quasars_vs_rest) <- c("rest_pr", "QUA_pr")

print("For train data")
print(CM_tr_quasars_vs_rest)

GT <- factor(c("Rest", "Quasars", "Rest", "Quasars"))
predicted <- factor(c("Rest", "Rest", "Quasars", "Quasars"))
values <- c(CM_tr_quasars_vs_rest)
df <- data.frame(GT, predicted, values)

# Visualise this CM
plot <- ggplot(data =  df, mapping = aes(x = predicted, y = GT)) +
  geom_tile(aes(fill = values), color = "white") +
  geom_text(data = subset(df, values > 40000),
            aes(label = values), 
            vjust = 1, color = "white", fontface = "bold", family = "arial") +
  geom_text(data = subset(df, values < 40000),
            aes(label = values), 
            vjust = 1, color = "black", fontface = "bold", family = "arial") +
  scale_fill_gradient(low = "#bbdefb", high = "#0d47a1") +
  theme_bw() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev) +
  labs(x = "Predicted", y = "Ground Truth", fill = "# cases",
       title = "Train CM for quasars-versus-rest case") +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 10, face = "plain"))

print(plot)

CM_te_quasars_vs_rest <- matrix(c(CM_te[1, 1] + CM_te[1, 3] + CM_te[3, 1] + CM_te[3, 3], 
                                CM_te[1, 2] + CM_te[3, 2],
                                CM_te[2, 1] + CM_te[2, 3],
                                CM_te[2, 2]), byrow = TRUE, ncol = 2)

colnames(CM_te_quasars_vs_rest) <- c("rest_tr", "QUA_tr")
rownames(CM_te_quasars_vs_rest) <- c("rest_pr", "QUA_pr")

print("For test data")
print(CM_te_quasars_vs_rest)

GT <- factor(c("Rest", "Quasars", "Rest", "Quasars"))
predicted <- factor(c("Rest", "Rest", "Quasars", "Quasars"))
values <- c(CM_te_quasars_vs_rest)
df <- data.frame(GT, predicted, values)

# Visualise this CM
plot <- ggplot(data =  df, mapping = aes(x = predicted, y = GT)) +
  geom_tile(aes(fill = values), color = "white") +
  geom_text(data = subset(df, values > 6000),
            aes(label = values), 
            vjust = 1, color = "white", fontface = "bold", family = "arial") +
  geom_text(data = subset(df, values < 6000),
            aes(label = values), 
            vjust = 1, color = "black", fontface = "bold", family = "arial") +
  scale_fill_gradient(low = "#bbdefb", high = "#0d47a1") +
  theme_bw() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev) +
  labs(x = "Predicted", y = "Ground Truth", fill = "# cases",
       title = "Test CM for quasars-versus-rest case") +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 10, face = "plain"))

print(plot)

# Here, we observe the large number of misclassification cases. To assess the
# model performance, considering the data imbalance (low number of quasars wrt
# other classes, as shown in Figure X), we use the Phi coefficient, which is a
# recommended choice for binary classification in the case of imbalanced data

phi_coef <- function(M){
  numerator <- (M[1, 1]*M[2, 2] - M[1, 2]*M[2, 1])
  
  MCC <- numerator/sqrt(M[1, 1] + M[1, 2])
  MCC <- MCC/sqrt(M[2, 1] + M[2, 2])
  MCC <- MCC/sqrt(M[1, 1] + M[2, 1])
  MCC <- MCC/sqrt(M[1, 2] + M[2, 2])
  return(MCC)
}

cat("Train Phi coefficienct for distinguishing between quasars and other classes\n", 
    round(phi_coef(CM_tr_quasars_vs_rest)*100, 1), "%", sep = "")
cat("Test Phi coefficienct for distinguishing between quasars and other classes\n", 
    round(phi_coef(CM_te_quasars_vs_rest)*100, 1), "%", sep = "")

# Considering the considerable imbalance, the model performance is decent in
# distinguishing quasars from other classes

# Finally, we want to assess the model performance in distinguishing galaxies 
# from other classes. 

CM_tr_galaxies_vs_rest <- matrix(c(CM_tr[2, 2] + CM_tr[2, 3] + CM_tr[3, 2] + CM_tr[3, 3], 
                                  CM_tr[1, 2] + CM_tr[1, 3],
                                  CM_tr[2, 1] + CM_tr[3, 1],
                                  CM_tr[1, 1]), byrow = TRUE, ncol = 2)

colnames(CM_tr_galaxies_vs_rest) <- c("rest_tr", "GAL_tr")
rownames(CM_tr_galaxies_vs_rest) <- c("rest_pr", "GAL_pr")

print("For train data")
print(CM_tr_galaxies_vs_rest)
GT <- factor(c("Rest", "Galaxies", "Rest", "Galaxies"))
predicted <- factor(c("Rest", "Rest", "Galaxies", "Galaxies"))
values <- c(CM_tr_galaxies_vs_rest)
df <- data.frame(GT, predicted, values)

# Visualise this CM
plot <- ggplot(data =  df, mapping = aes(x = predicted, y = GT)) +
  geom_tile(aes(fill = values), color = "white") +
  geom_text(data = subset(df, values > 30000),
            aes(label = values), 
            vjust = 1, color = "white", fontface = "bold", family = "arial") +
  geom_text(data = subset(df, values < 30000),
            aes(label = values), 
            vjust = 1, color = "black", fontface = "bold", family = "arial") +
  scale_fill_gradient(low = "#bbdefb", high = "#0d47a1") +
  theme_bw() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev) +
  labs(x = "Predicted", y = "Ground Truth", fill = "# cases",
       title = "Train CM for galaxies-versus-rest case") +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 10, face = "plain"))

print(plot)

CM_te_galaxies_vs_rest <- matrix(c(CM_te[2, 2] + CM_te[2, 3] + CM_te[3, 2] + CM_te[3, 3], 
                                  CM_te[1, 2] + CM_te[1, 3],
                                  CM_te[2, 1] + CM_te[3, 1],
                                  CM_te[1, 1]), byrow = TRUE, ncol = 2)

colnames(CM_te_galaxies_vs_rest) <- c("rest_tr", "GAL_tr")
rownames(CM_te_galaxies_vs_rest) <- c("rest_pr", "GAL_pr")

print("For test data")
print(CM_te_galaxies_vs_rest)
GT <- factor(c("Rest", "Galaxies", "Rest", "Galaxies"))
predicted <- factor(c("Rest", "Rest", "Galaxies", "Galaxies"))
values <- c(CM_te_galaxies_vs_rest)
df <- data.frame(GT, predicted, values)

# Visualise this CM
plot <- ggplot(data =  df, mapping = aes(x = predicted, y = GT)) +
  geom_tile(aes(fill = values), color = "white") +
  geom_text(data = subset(df, values > 6000),
            aes(label = values), 
            vjust = 1, color = "white", fontface = "bold", family = "arial") +
  geom_text(data = subset(df, values < 6000),
            aes(label = values), 
            vjust = 1, color = "black", fontface = "bold", family = "arial") +
  scale_fill_gradient(low = "#bbdefb", high = "#0d47a1") +
  theme_bw() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev) +
  labs(x = "Predicted", y = "Ground Truth", fill = "# cases",
       title = "Test CM for galaxies-versus-rest case") +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 10, face = "plain"))

print(plot)
# For assessing the model performance for distinguishing galaxies from other classes,
# we can use the accuracy metric because the imbalance is not a problem here. 

cat("Train accuracy for distinguishing between galaxies and other classes\n")
MCR_tr <- (CM_tr_galaxies_vs_rest[1, 1] + CM_tr_galaxies_vs_rest[2, 2]) / sum(CM_tr_galaxies_vs_rest) 
cat(round(MCR_tr*100, 1), "%", sep="")

cat("Test accuracy for distinguishing between galaxies and other classes\n")
MCR_te <- (CM_te_galaxies_vs_rest[1, 1] + CM_te_galaxies_vs_rest[2, 2]) / sum(CM_te_galaxies_vs_rest) 
cat(round(MCR_te*100, 1), "%", sep="")
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

plot <- ggplot(data = df_classes, mapping = aes(x = counts, y = reorder(classes, counts))) + 
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

 