
#########################################################

#AFOLA KOSSI MAWOUENA SAMUEL
# July 21, 2023
# Prediction Cardiovascular disease

source('HLtest.r')
source('val.prob.ci.dec08.r')

source('dca.r')

# Load the required package
#library(survival)
#?coxph

source("val.prob.ci.dec08.r")

smart <- readRDS(file = "SMARTs_P2-3.rds")

smart <- na.omit(smart) 

###################################################################################
#							 #
# CARDIOVASCULAR DISEASE     #
#############################
valother  <- smart[smart$outcome==1 | smart$outcome==0,]



##########################################Training data##################
set.seed(150283)  # Set a seed for reproducibility

# Define the proportion for each set (adjust as needed)
train_prop <- 0.7  # Proportion for the training set
val_prop <- 0.15  # Proportion for the validation set
dev_prop <- 0.15  # Proportion for the development set

# Calculate the number of observations for each set
n <- nrow(smart)
train_size <- round(train_prop * n)
val_size <- round(val_prop * n)
dev_size <- n - train_size - val_size

# Randomly shuffle the indices of the dataset
indices <- sample(1:n)

# Create the train, validation, and development sets based on the shuffled indices
train <- smart[indices[1:train_size], ]
val <- smart[indices[(train_size + 1):(train_size + val_size)], ]
full <- smart[indices[(train_size + val_size + 1):n], ]




valother  <- smart[smart$outcome==1 | smart$outcome==0,]




##########################################Training data##################
set.seed(150283)  # Set a seed for reproducibility

# Define the proportion for each set (adjust as needed)
train_prop <- 0.7  # Proportion for the training set
val_prop <- 0.15  # Proportion for the validation set
dev_prop <- 0.15  # Proportion for the development set

# Calculate the number of observations for each set
n <- nrow(smart)
train_size <- round(train_prop * n)
val_size <- round(val_prop * n)
dev_size <- n - train_size - val_size

# Randomly shuffle the indices of the dataset
indices <- sample(1:n)

# Create the train, validation, and development sets based on the shuffled indices
train <- smart[indices[1:train_size], ]
val <- smart[indices[(train_size + 1):(train_size + val_size)], ]
full <- smart[indices[(train_size + val_size + 1):n], ]



##############################################################################

library(ggplot2)
library(plyr)
library(MASS)
library(reshape2)
install.packages("lmtest")
library(lmtest)
install.packages("rms")
library(rms)
library(pROC)
install.packages("doParallel")
library(dplyr)
library(doParallel)
library(prodlim)
library(survival)
library(Hmisc)
library(MASS)

library(Hmisc)
library(MASS)
library(foreign)




full <- glm(CARDIAC ~ AGE + DIABETES + CEREBRAL + BMIO + SBP,
            data = train, family="binomial" )



# Convert full$y to numeric

full$y <- as.numeric(full$y)

full 

B     <- mean((full$y) * (1-plogis(full$linear.predictors))^2 + 
                (1-full$y) * plogis(full$linear.predictors)^2)
B
Bmax  <- mean(full$y) * (1-mean(full$y))
Bmax
Bscaled <- 1 - B/Bmax
Bscaled

# Compare to Pearson R2
cor(x=plogis(full$linear.predictors), y=full$y)^2

#######################################


# Conversion

val$CARDIAC <- as.numeric(val$CARDIAC)

###################################### Model validation

lp.val <- predict(object = full, newdata = val)

# Brier max
B     <- mean((val$CARDIAC) * (1-plogis(lp.val))^2 + 
                (1-val$CARDIAC) * plogis(lp.val)^2)
B
Bmax  <- mean(val$CARDIAC) * (1-mean(val$CARDIAC))
Bmax
Bscaled <- 1 - B/Bmax
Bscaled
cor(x=plogis(lp.val), y=val$CARDIAC)^2





library(Hmisc)
library(MASS)
library(foreign)




full <- glm(CARDIAC ~ AGE + DIABETES + CEREBRAL + BMIO + SBP,
            data = train, family="binomial" )
########################################################################



# Convert full$y to numeric

full$y <- as.numeric(full$y)

full 

B     <- mean((full$y) * (1-plogis(full$linear.predictors))^2 + 
                (1-full$y) * plogis(full$linear.predictors)^2)
B
Bmax  <- mean(full$y) * (1-mean(full$y))
Bmax
Bscaled <- 1 - B/Bmax
Bscaled

# Compare to Pearson R2
cor(x=plogis(full$linear.predictors), y=full$y)^2

#######################################


# Conversion

val$CARDIAC <- as.numeric(val$CARDIAC)

###################################### Model validation

lp.val <- predict(object = full, newdata = val)

# Brier max
B     <- mean((val$CARDIAC) * (1-plogis(lp.val))^2 + 
                (1-val$CARDIAC) * plogis(lp.val)^2)
B
Bmax  <- mean(val$CARDIAC) * (1-mean(val$CARDIAC))
Bmax
Bscaled <- 1 - B/Bmax
Bscaled
cor(x=plogis(lp.val), y=val$CARDIAC)^2
########################################################





# Conversion

val$CARDIAC <- ifelse(val$CARDIAC == 1, 0, 1)




#####################################

## H-L tests
hl.ext2(p=plogis(full$linear.predictor),y=full$y,g=10,df=8)
hl.ext2(p=plogis(lp.val),y=val$CARDIAC,g=10,df=9)

## CI around c stat
cstatNo <- rcorr.cens(full$linear.predictors, full$y) 
cat(cstatNo[1], "[", cstatNo[1]-1.96/2*cstatNo[3], " - ", cstatNo[1]+1.96/2*cstatNo[3],"]")
cstatvAL <- rcorr.cens(lp.val, val$CARDIAC) 
cat(cstatvAL[1], "[", cstatvAL[1]-1.96/2*cstatvAL[3], " - ", cstatvAL[1]+1.96/2*cstatvAL[3],"]")

#############################
## Validation plots
par(mfrow = c(1,2), pty='s', mar=c(4.2,4,4,1),cex=0.95, font =1, col=1 )
# Apparent
val.prob.ci(logit=full$linear.predictor,y=full$y, pl=T,smooth=T,logistic.cal=F, g=10,
            xlab="Predicted risk CARDIAC",
            ylab="Observed masses with CARDIAC",riskdist='predicted',
            d1lab="CARDIAC", d0lab="Necrosis", dist.label=-0.95, cutoff=.2)
title("Development, n=2468")

# External
val.prob.ci(logit=lp.val,y=val$CARDIAC, pl=T,smooth=T,logistic.cal=F, g=10,
            xlab="Predicted risk CARDIAC",
            ylab="Observed masses with tumor at validation",riskdist='predicted',
            d1lab="Tumor", d0lab="Necrosis", dist.label=-0.95, cutoff=.2)
title("Validation, n=529")


#############################
## Boxplots

par(mfrow = c(1,2), pty='s', mar=c(4.2,4,4,1),cex=0.95, font =1, col=1 )
boxplot(plogis(full$linear.predictors)~full$y,
        ylab="Predicted risk without LDH", xlab="Tumor",ylim=c(0,1))
boxplot(c(mean(plogis(full$linear.predictors[full$y==0])), mean(plogis(full$linear.predictors[full$y==1])))~ c(0,1),add=T,
        boxlty=0, staplelty=0, medlty=0, medlwd=0, medpch=15)
title(paste("Development: Slope=", 
            round(mean(plogis(full$linear.predictors[full$y==1])) - 
                    mean(plogis(full$linear.predictors[full$y==0])),2),sep=""))

boxplot(plogis(lp.val)~val$CARDIAC,
        ylab="Predicted risk without CARDIAC", xlab="CARDIAC at validation",ylim=c(0,1))
boxplot(c(mean(plogis(lp.val[val$CARDIAC==0])), mean(plogis(lp.val[val$CARDIAC==1])))~ c(0,1),add=T,
        boxlty=0, staplelty=0, medlty=0, medlwd=0, medpch=15)
title(paste("Validation: Slope=", 
            round(mean(plogis(lp.val[val$CARDIAC==1])) - 
                    mean(plogis(lp.val[val$CARDIAC==0])),2),sep=""))



###########################


# External validation
plot(x=dcaVal$threshold, y=dcaVal[,1], type='l', lty=1, lwd=2, las=1, ylim=c(-.05,.8),
     ylab="Net Benefit at validation", xlab="Threshold risk for resection of cardiac (%)", cex.lab=1.2)
lines(x=dcaVal$threshold, y=dcaVal[,2], lty=3, lwd=1)
lines(x=dcaVal$threshold, y=dcaVal[,3], lty=1, lwd=1)

arrows(20,max(dcaVal[,2],na.rm=T)-.05,20,0)
text(x=20,y=max(dcaVal[,2],na.rm=T), "Threshold")
text(x=80,y=.03, "Treat none")
text(x=60,y=.2, "Treat\nin all")
text(x=75,y=.5, "Resection\nif tumor risk\n> threshold")
title("Validation, n=529")





########################
# Internal validation  #
val.full  <- validate(full, B=200)
val.full
val.full[1,1:5]/2 + .5  # index.corrected c stat 0.811

# Bootstrap discrimination slope and DCA
nrowB	<- nrow(train)  # nrow from development set
B <- 20            # 200 bootstraps
matB <- matrix(NA,nrow=B,ncol=3) # Matrix for results
dimnames(matB) <- list(c(1:B), Cs(Slopeapp, Slopetest, optimism ))
matDCA  <- matrix(0,nrow=99,ncol=3) # Matrix for results of DCA
dimnames(matDCA) <- list(c(1:99), Cs(NBorig, NBval, NBoptimism ))

# Start loop
for (i in 1:B) {
  if (i%%10==0) cat("Start Bootstrap sample nr", i, "\n")
  Brows <- sample(nrowB,replace=T)
  
  # Bsample is bootstrap sample from development set
  Bsample	<- train[Brows,]
  # devfull <- lrm(CARDIAC ~ AGE+DIABETES+CEREBRAL+BMIO+SBP, data=Bsample,linear.predictors=T, x=T, y=T)
  devfull <- glm(CARDIAC ~ AGE+DIABETES+CEREBRAL+BMIO+SBP,family="binomial" )
  
  
  
  matB[i,1] <- mean(plogis(devfull$linear.predictors[devfull$y==1])) - 
    mean(plogis(devfull$linear.predictors[devfull$y==0]))
  lp  <- full$x %*% devfull$coef[2:length(devfull$coef)] + devfull$coef[1] # lp with coefs from bootstrap
  matB[i,2] <- mean(plogis(lp[full$y==1])) - mean(plogis(lp[full$y==0]))  # Testing on original sample
  
  dcaorig   <- dca(yvar=devfull$y, xmatrix=plogis(devfull$linear.predictor), prob="Y") # bootstrap sample
  dcaval    <- dca(yvar=full$y, xmatrix=plogis(devfull$linear.predictor), prob="Y") # Testing on original sample
  matDCA[,1]  <- (i-1)/i * matDCA[,1] + 1/i * (dcaorig[,1] - dcaorig[,2])  # NB orig
  matDCA[,2]  <- (i-1)/i * matDCA[,2] + 1/i * (dcaval[,1]  - dcaval[,2])   # NB validated
  
} 

### End example performance of prediction models ###
#################################################################################################

#############################################################################################
##				##
##				##
###########################  Future Work #######################################################






smart <- readRDS(file= "SMARTs_P1.rds")
class(smart) # This tells us the class/data type of the object smart (hopefully returning data.frame)
smart <- na.omit(smart)
dim(smart)
sapply(smart,class) # This tells us the class of each variable in the smart data, eg, categorical='factor' or continuous='numeric'


##############################################


# Albumin
smart$albumin <- as.factor(smart$albumin)
levels(smart$albumin)
smart$albumin <- revalue(smart$albumin, c("1"="No", "2"="Low", "3"="High"))
levels(smart$albumin)


#############################################


head(smart)

#########################################
table(is.na(smart$SYSTH)==FALSE | is.na(smart$SYSTBP)==FALSE)


# Load required libraries
library(caret)

# Set seed for reproducibility
set.seed(123)

# Split the data into training and testing sets
train_indices <- sample(1:nrow(smart), 0.7*nrow(smart))
train_data <- smart[train_indices, ]
test_data <- smart[-train_indices, ]

# Define the predictor variables
predictors <- c("TEVENT", "EVENT", "SEX", "AGE", "DIABETES", "CEREBRAL", "AAA", "PERIPH", "STENOSIS", "SYSTBP", "DIASTBP", "SYSTH", "DIASTH", "LENGTHO", "WEIGHTO", "BMIO", "CHOLO", "albumin", "SMOKING", "packyrs", "alcohol")

# Create a formula for the model
formula <- as.formula(paste("CARDIAC ~", paste(predictors, collapse = "+")))

# Train the model using random forest algorithm
model <- train(formula, data = train_data, method = "rf")

# Make predictions on the test data
predictions <- predict(model, newdata = test_data)

# Evaluate the model performance
confusionMatrix(predictions, test_data$CARDIAC)


# Evaluate the model performance
confusionMatrix(predictions, test_data$CARDIAC)


#Create the train , val, dev
heart.split <- initial_split(smart)
heart.train <- training(heart.split)
heart.test <- testing(heart.split)


# Devellpment 
heart.full <- glm(CARDIAC~., data = heart.train, family = "binomial")
summary(heart.full)


# set engine
heart_model <- logistic_reg() %>%
  set_engine("glm")

# create recipe
heart_recipe <- recipe(CARDIAC ~., data = heart.train) %>%
  step_zv(all_predictors())

# build work flow
heart_wflow <- workflow() %>%
  add_model(heart_model) %>%
  add_recipe(heart_recipe)

# fit training data through the work flow 
heart_fit <- heart_wflow %>%
  fit(data = heart.train)
tidy(heart_fit)
######################################################



# Fit the logistic regression model
logit_model <- glm(outcome ~ ., data = heart.train, family = "binomial")

# Predict probabilities for the test data
probabilities <- predict(logit_model, newdata = heart.test, type = "response")

# Create a roc object
roc_obj <- roc(response = heart.test$outcome, predictor = probabilities)

# Plot the roc curve
plot(roc_obj, main = "ROC Curve", xlab = "Specificity",col = "blue", ylab = "Sensitivity")

# Compute the area under the roc curve
roc_area <- auc(roc_obj)
roc_area


########################################

set.seed(470)
folds <- vfold_cv(heart.train, v=5)

heart_fit_rs <- heart_wflow %>%
  fit_resamples(folds)

metrics <- data.frame(collect_metrics(heart_fit_rs, summarize = FALSE))

metrics <- metrics %>%
  select(-.config)
colnames(metrics) <- c("Fold", "Metric", "Estimator", "Estimate")
metrics

##################################################
heart_disease_pred <- predict(heart_fit, new_data = heart.test) %>%
  bind_cols(heart.test %>% select(CARDIAC))

test_accuracy <- accuracy(heart_disease_pred, truth = CARDIAC, estimate = .pred_class)
test_specificity <- spec(heart_disease_pred, truth = CARDIAC, estimate = .pred_class)
test_sensitivity <- sens(heart_disease_pred, truth = CARDIAC, estimate = .pred_class)

test.values <- data.frame(test_accuracy$.estimate, test_sensitivity$.estimate, test_specificity$.estimate)
colnames(test.values) <- c("Test set Accuracy", "Test set Sensitivity", "Test set Specificity")
test.values



################################################################

# Compute the sample size
sample_size <- nrow(heart.train)
sample_size



# Set the seed for reproducibility
set.seed(470)

# Perform 5-fold cross-validation
folds <- vfold_cv(heart.train, v = 5)

# Fit the workflow using resampling
heart_fit_rs <- heart_wflow %>%
  fit_resamples(folds)

# Collect the metrics from the resampled fits
metrics <- data.frame(collect_metrics(heart_fit_rs, summarize = FALSE))

# Remove the .config column from the metrics dataframe
metrics <- metrics %>%
  select(-.config)

# Rename the columns of the metrics dataframe
colnames(metrics) <- c("Fold", "Metric", "Estimator", "Estimate")

# Print the metrics dataframe
metrics

# Predict the outcome for the test set using the fitted model
heart_disease_pred <- predict(heart_fit, new_data = heart.test) %>%
  bind_cols(heart.test %>% select(CARDIAC))

# Calculate the accuracy, sensitivity, and specificity for the test set predictions
test_accuracy <- accuracy(heart_disease_pred, truth = CARDIAC, estimate = .pred_class)
test_specificity <- spec(heart_disease_pred, truth = CARDIAC, estimate = .pred_class)
test_sensitivity <- sens(heart_disease_pred, truth = CARDIAC, estimate = .pred_class)

# Create a dataframe to store the test set evaluation metrics
test.values <- data.frame(test_accuracy$.estimate, test_sensitivity$.estimate, test_specificity$.estimate)

# Rename the columns of the test.values dataframe
colnames(test.values) <- c("Test set Accuracy", "Test set Sensitivity", "Test set Specificity")

# Print the test.values dataframe
test.values

# Calculate the sample size of the training data
sample_size <- nrow(heart.train)

# Print the sample size
sample_size


###################################

# Load the required library
library(caret)

# Convert outcome and EVENT variables to factors with the same levels
heart.test$outcome <- factor(heart.test$outcome)
heart.test$EVENT <- factor(heart.test$EVENT)

# Create a confusion matrix
confusion_matrix <- confusionMatrix(data = heart.test$outcome, reference = heart.test$EVENT)

# Print the confusion matrix
confusion_matrix
