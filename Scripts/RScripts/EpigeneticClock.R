#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#	Bedmethyl.R
#
#	Purpose 
#		This code was made for creating epigenetic clock - age prediction model
#   Using machine learning model (2023, Anastasiadi and Preffer... )
#
# Author comment:
#	Lisa Koole
#	lisa.koole@maastrichtuniversity.nl
#
#
#-----------------------------------------------------------------------------------------------------#
#							Main settings or load
#-----------------------------------------------------------------------------------------------------#
# Get all general settings 

s_ROOT_dir <<- "G:/.shortcut-targets-by-id/1-U3PYANHthlhUCklfcHOQA8_gPsFjKjA/Lisa Koole/Data (raw) and analysis/Oxford Nanopore Epi-AD/"
source(paste0(s_ROOT_dir,"Scripts/.Main/Settings.R"))


#-----------------------------------------------------------------------------------------------------#
#							        Libraries needed
#-----------------------------------------------------------------------------------------------------#
library(tidymodels)
library(readr)


#-----------------------------------------------------------------------------------------------------#
#					            	Load data
#-----------------------------------------------------------------------------------------------------#

load(paste0(s_OUT_dir,"Bedmethyl/methyl_all.rdata"))
load(paste0(s_ROOT_dir,s_out_folder,"Pheno/Pheno.Rdata")) 


#-----------------------------------------------------------------------------------------------------#
#					            Switching columns and rows 
#-----------------------------------------------------------------------------------------------------#

methyl_flipped <- data.frame(t(methyl_all))
methyl_flipped <- cbind(methyl_flipped, "age" = pheno$Age)
rm(methyl_all)

#-----------------------------------------------------------------------------------------------------#
#					            Splitting dataset
#-----------------------------------------------------------------------------------------------------#

set.seed(123)
splits <- initial_split(methyl_flipped, strata = age)

age_other <- training(splits)
age_test <- testing(splits)

# Training set proportions by age class
age_other %>%
  count(age) %>%
  mutate(prop = n/sum(n))

# Test set proportions by age class
age_test %>%
  count(age) %>%
  mutate(prop = n/sum(n))


#-----------------------------------------------------------------------------------------------------#
#					    Exclusion CpGs with zero variance across ages
#-----------------------------------------------------------------------------------------------------#

library(caret)
library(dplyr)

## Detect features and visualize them
nzv.cpg <- nearZeroVar(age_other, saveMetrics= TRUE,
                       names=TRUE, freqCut = 85/15, uniqueCut = 50)
boxplot(nzv.cpg$percentUnique)
boxplot(nzv.cpg$freqRatio)

## Detect features, exclude them and save the object
nzv.cpg.list <- nearZeroVar(age_other, freqCut = 85/15,
                            uniqueCut = 50) filteredDescr <- age_other[, -nzv.cpg.list]
dim(filteredDescr)

#-----------------------------------------------------------------------------------------------------#
#					    Exclusion highly correlated variables
#-----------------------------------------------------------------------------------------------------#

highlyCorDescr <- findCorrelation(filteredDescr, cutoff = 0.8)
filteredDescr.cor <- filteredDescr[,-highlyCorDescr]

#-----------------------------------------------------------------------------------------------------#
#					   Transformation, mean= 0 and sd = 1
#-----------------------------------------------------------------------------------------------------#

preProcValues <- preProcess(filteredDescr.cor, method = c
                            (“center”, “scale”))
trainTransformed <- predict(preProcValues, filteredDescr.cor)

#-----------------------------------------------------------------------------------------------------#
#					   Imputation of missing values
#-----------------------------------------------------------------------------------------------------#

library(mice)

init = mice(meth.age.df, maxit=0)
meth = init$method
predM = init$predictorMatrix
colnames(meth.age.df)
predM[, c("age")]=0
meth[c("age")]=""
set.seed(100)
imputed = mice(meth.age.df, method=meth,predictorMatrix=predM, m=5)

#-----------------------------------------------------------------------------------------------------#
#					   Prediction Model tuning and evaluation
#-----------------------------------------------------------------------------------------------------#

# repeated cross-validation
fitControl <- trainControl(method = ‘repeatedcv ’,
                           number=10, repeats=10)

# Define range of lambda to be tested
lambda <- 10^seq(-3, 3, length = 100)

# run penalized regressions
set.seed(123)
elastic_model <- train(age ~., data = trainTransformed, method =
                         “glmnet”, trControl = fitControl, tuneLength = 10)


models_compare <- resamples(list(R=ridge_model,
                                 LM=lasso_model, EM=elastic_model, EM05=elastic_model.05))
summary(models_compare)

# sum of CpGs 
sum(coef(elastic_model$finalModel, elastic_model$bestTune
         $lambda)!=0)


#-----------------------------------------------------------------------------------------------------#
#					   Prediction Model evaluation in training dataset
#-----------------------------------------------------------------------------------------------------#

predicted.age <- predict.train(elastic_model)
postResample(pred = predicted.age, trainTransformed$age)
cor.test(predicted.age, trainTransformed$age)


#-----------------------------------------------------------------------------------------------------#
#					   Prediction Model evaluation in testing dataset
#-----------------------------------------------------------------------------------------------------#

predict.enet.test <- predict(elastic_model, testTransformed)
postResample(pred = predict.enet.test, testTransformed$age)
cor.test(predict.enet.test, testTransformed$age)
     

#-----------------------------------------------------------------------------------------------------#
#					  Build final modeel
#-----------------------------------------------------------------------------------------------------#

inalmodelCtrl <- trainControl(method = “none”)
set.seed(123)
final <- train(age ~., data = trainTransformed, method =
                 "glmnet", trControl=finalmodelCtrl, tuneGrid = expand.grid(alpha
                                                                            = bestalpha, lambda = bestlambda))
predicted.final.train <- predict(final, trainTransformed)
cor.test(predicted.final.train, trainTransformed$age)        
         