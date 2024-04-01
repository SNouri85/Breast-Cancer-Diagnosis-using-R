library("ggplot2")
library("corrgram")
library("car")
library("lattice")
library("ROCR")
library("plotly")
library("tree")
library(devtools)
library(tidyverse)
library(corrplot)
library(PerformanceAnalytics)
library(psych)
library(GGally)



dt <- read.csv("F:\\Project\\datasets_180_408_data.csv", header=TRUE, sep=",")
View(dt)
colnames(dt)
dim(dt)
table(dt$diagnosis)
library(dplyr)
glimpse(dt)
#----------------------------------------------------------------------
#check for missing value
sapply(dt, function(x) sum(is.na(x)))
library(dplyr)
dt <- dt %>% select(-c("id","X"))
ncol(dt)
#----------------------------------------------------------------------
summary(dt)
summary(dt$diagnosis)
table(dt$diagnosis)
#-----------------------------------------------------------------------
dt$y[dt$diagnosis=="M"] = 1
dt$y[dt$diagnosis=="B"] = 0
#---------------------------------------------------------------------
#Correlation chart for means
library(PerformanceAnalytics)
chart.Correlation(dt[, c(3:12)], histogram=TRUE, col="grey10", pch=1, main="Cancer Means")
library(psych)
pairs.panels(dt[,c(3:12)], method="pearson",
             hist.col = "#1fbbfa", density=TRUE, ellipses=TRUE, show.points = TRUE,
             pch=1, lm=TRUE, cex.cor=1, smoother=F, stars = T, main="Means")

#correlation chart for SE
chart.Correlation(dt[, c(13:22)], histogram=TRUE, col="grey10", pch=1, main="Cancer SEs")
pairs.panels(dt[,c(13:22)], method="pearson",
             hist.col = "#1fbbfa", density=TRUE, ellipses=TRUE, show.points = TRUE,
             pch=1, lm=TRUE, cex.cor=1, smoother=F, stars = T, main="SEs")


#correlation chart for worst
chart.Correlation(dt[, c(23:32)], histogram=TRUE, col="grey10", pch=1, main="Cancer SEs")
pairs.panels(dt[,c(23:32)], method="pearson",
             hist.col = "#1fbbfa", density=TRUE, ellipses=TRUE, show.points = TRUE,
             pch=1, lm=TRUE, cex.cor=1, smoother=F, stars = T, main="Worst")

#-----------------------------------------------------------------------
library(corrplot)
corMat <- cor(dt[,2:31])
corrplot(corMat, order = "hclust", tl.cex = 0.7, col = terrain.colors(200))
#dropping covariates with more than .9 corr
library(caret)
highlyCor <- colnames(dt)[findCorrelation(corMat, cutoff = 0.9, verbose = TRUE)]
print(highlyCor)

#10 columns are flagged for removal.
dt_cor <- dt[, which(!colnames(dt) %in% highlyCor)]
ncol(dt_cor)

#-----------------------------------------------------------------------
#split data train and test subsets (70/30)
sample_size = floor(0.7 * nrow(dt))

# set the seed to make your partition reproductible
set.seed(1729)
train_set = sample(seq_len(nrow(dt)), size = sample_size)
training = dt[train_set, ]
testing = dt[-train_set, ]
head(training)
nrow(training)
nrow(testing)

#---------------------------------------------------------------------
set.seed(1234)
df <- cbind(diagnosis=dt$y, dt_cor)
train_indx <- createDataPartition(df$diagnosis, p = 0.7, list = FALSE)

train_set <- df[train_indx,]
test_set <- df[-train_indx,]

nrow(train_set)
nrow(test_set)

#-----------------------------------------------------------------------
library("ggplot2")
ggplot(dt, aes(x = diagnosis)) +
  geom_bar(aes(fill = "blue")) +
  ggtitle("Distribution of diagnosis") +
  theme(legend.position="none")
#for training set
ggplot(training, aes(x = diagnosis)) + 
  geom_bar(aes(fill = "blue")) + 
  ggtitle("Distribution of diagnosis for the training subset") + 
  theme(legend.position="none")
#for testing set
ggplot(testing, aes(x = diagnosis)) + 
  geom_bar(aes(fill = "blue")) + 
  ggtitle("Distribution of diagnosis for the testing subset") + 
  theme(legend.position="none")
#----------------------------------------------------------------------
#fitting GLM on the training data
model = glm(diagnosis ~. ,family=binomial(link='logit'), control = list(maxit = 50),data=train_set)
print(summary(model))
print(anova(model, test="Chisq"))

#selecting variables based on backward elimination method
step(model, direction = "backward", trace=FALSE )

#final model for training using selected variables
model_final=glm(formula = diagnosis ~ smoothness_mean + concave.points_mean + 
      area_se + concavity_se + smoothness_worst + compactness_worst + 
      concavity_worst + symmetry_worst, family = binomial(link = "logit"), 
    data = train_set, control = list(maxit = 50))


#---------------------------------------------------------------------
#prdictions using final model for training set
prediction_training = predict(model_final,train_set, type = "response")
prediction_training = ifelse(prediction_training > 0.5, 1, 0)
error = mean(prediction_training != train_set$diagnosis)
print(paste('Model Accuracy',1-error))


#ROCcurve for training
p = predict(model_final, train_set, type="response")
library(ROCR)
pr = prediction(p, train_set$diagnosis)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
#AUC for training
auc = performance(pr, measure = "auc")
auc = auc@y.values[[1]]
print(paste("Model Accuracy", auc))
#----------------------------------------------------------------------
#prdictions using final model for testing set
prediction_testing = predict(model_final,test_set, type = "response")
prediction_testing = ifelse(prediction_testing > 0.5, 1, 0)
error = mean(prediction_testing != test_set$diagnosis)
print(paste('Model Accuracy',1-error))


#ROCcurve for testing
p = predict(model_final, test_set, type="response")
library(ROCR)
pr = prediction(p, test_set$diagnosis)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
#AUC for training
auc = performance(pr, measure = "auc")
auc = auc@y.values[[1]]
print(paste("Model Accuracy", auc))
#--------------------------------------------------------------------
#Using PCA
dt.pca <- prcomp(dt[, 2:31], center=TRUE, scale=TRUE)
plot(dt.pca, type="l", main='')
grid(nx = 10, ny = 14)
title(main = "Principal components ", sub = NULL, xlab = "Components")
box()

pca_var <- dt.pca$sdev^2
pve_df <- pca_var / sum(pca_var)
cum_pve <- cumsum(pve_df)
pve_table <- tibble(comp = seq(1:ncol(select(dt,-c("diagnosis","y")))), pve_df, cum_pve)
ggplot(pve_table, aes(x = comp, y = cum_pve)) + 
   geom_point(shape=22, fill="blue", color="blue", size=3) + 
  geom_abline(intercept = 0.95, color = "green", slope = 0)+xlab("Principal components") +
     ylab("cumulative variance") +
     ggtitle("Cumulative variance explained by principal components")

#selecting the first 10 PC for analysis
pca_df <- as.data.frame(dt.pca$x)
pca_dff <- pca_df[,1:10]
View(pca_dff)
#---------------------------------------------------------------------
set.seed(1234)
library(caret)
dfpca <- cbind(diagnosis=dt$y, pca_dff)
trainpca_indx <- createDataPartition(dfpca$diagnosis, p = 0.7, list = FALSE)

trainpca_set <- dfpca[trainpca_indx,]
testpca_set <- dfpca[-trainpca_indx,]

nrow(trainpca_set)
nrow(testpca_set)
#---------------------------------------------------------------------
#fitting GLM on the trainingpca data
modelpca = glm(diagnosis ~. ,family=binomial(link='logit'), control = list(maxit = 50),data=trainpca_set)
print(summary(modelpca))
print(anova(modelpca, test="Chisq"))

#---------------------------------------------------------------------
#prdictions using final model for trainingpca set
prediction_trainingpca = predict(modelpca,trainpca_set, type = "response")
prediction_trainingpca = ifelse(prediction_trainingpca > 0.5, 1, 0)
error = mean(prediction_trainingpca != trainpca_set$diagnosis)
print(paste('Modelpca Accuracy',1-error))

#ROCcurve for trainingpca
p = predict(modelpca, trainpca_set, type="response")
library(ROCR)
pr = prediction(p, trainpca_set$diagnosis)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
#AUC for training
auc = performance(pr, measure = "auc")
auc = auc@y.values[[1]]
print(paste("Area under the Roc curve", auc))
#---------------------------------------------------------------------
#fitting GLM on the testingpca data
modelpca = glm(diagnosis ~. ,family=binomial(link='logit'), control = list(maxit = 50),data=testpca_set)
print(summary(modelpca))
print(anova(modelpca, test="Chisq"))

#---------------------------------------------------------------------
#prdictions using final model for testingpca set
prediction_testingpca = predict(modelpca,testpca_set, type = "response")
prediction_testingpca = ifelse(prediction_testingpca > 0.5, 1, 0)
error = mean(prediction_testingpca != testpca_set$diagnosis)
print(paste('Modelpca Accuracy',1-error))

#ROCcurve for testingpca
p = predict(modelpca, testpca_set, type="response")
library(ROCR)
pr = prediction(p, trainpca_set$diagnosis)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
#AUC for training
auc = performance(pr, measure = "auc")
auc = auc@y.values[[1]]
print(paste("Area under the Roc curve", auc))
#--------------------------------------------------------------------
trainpca_set$diagnosis<-factor(trainpca_set$diagnosis,levels=c(0,1),labels=c("Benign","Malignant"))

df_control <- trainControl(method="cv",
                           number = 15,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary)

model_rf <- train(diagnosis ~., data = trainpca_set,
                     method = "rf", 
                     metric = 'ROC',
                     trControl = df_control)
#---------------------------------------------------------------------
#prediction on trainingpca_set
library(caret)
library(e1071)
prediction_rf <- predict(model_rf, trainpca_set)
cm_rf <- confusionMatrix(prediction_rf, trainpca_set$diagnosis, positive = "Malignant")
cm_rf

#prediction on testingpca_set
testpca_set$diagnosis<-factor(testpca_set$diagnosis,levels=c(0,1),labels=c("Benign","Malignant"))

prediction_rf <- predict(model_rf, testpca_set)
cm_rf <- confusionMatrix(prediction_rf, testpca_set$diagnosis, positive = "Malignant")
cm_rf
#--------------------------------------------------------------------
#logistic model
model_logreg <- train(diagnosis ~., data = trainpca_set, method = "glm", 
                         metric = "ROC", 
                         trControl = df_control)

prediction_logreg <- predict(model_logreg, testpca_set)
cm_logreg <- confusionMatrix(prediction_logreg, testpca_set$diagnosis, positive = "Malignant")
cm_logreg

#ROCcurve for testingpca

library(ROCR)
pr = prediction(prediction_rf, trainpca_set$diagnosis)
prf = performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
#AUC for training
auc = performance(pr, measure = "auc")
auc = auc@y.values[[1]]
print(paste("Area under the Roc curve", auc))
