#Part One
require(glmnet)
require(dplyr)
USPRC <- readr::read_csv("G:/UnitedStatesandPuertoRicoCancerStatistics.csv", col_names = TRUE)
head(USPRC)

USPRC <- USPRC %>%
  #remove Not Applicable
  filter(!Population %in% "Not Applicable") %>%
  mutate(Population = as.double(Population))
  #replace space
colnames(USPRC) <- gsub(" ", "_", colnames(USPRC)); colnames(USPRC)

USPRC$Leading_Cancer_Sites <- USPRC$Leading_Cancer_Sites %>% 
  gsub(" ", "_", .)
USPRC$States <- USPRC$States %>% 
  gsub(" ", "_", .)
USPRC$Race <- USPRC$Race %>% 
  gsub(" ", "_", .)

# as.factor(), 
Leading_Cancer_Sites <- as.factor(USPRC$Leading_Cancer_Sites)
States <- as.factor(USPRC$States)
Race <- as.factor(USPRC$Race)
Year <- as.factor(USPRC$Year)
Sex <- as.factor(USPRC$Sex)

# Crude_Rate
Crude_Rate <- as.double(USPRC$Crude_Rate)

factorVars <- model.matrix(Crude_Rate ~ Leading_Cancer_Sites + States + Race + Year + Sex)[, -1]
Count <- USPRC$Count
Population <- USPRC$Population

X <- data.frame(Count, Population, factorVars) %>% 
  as.matrix()

nr <- nrow(X); nr

#sample train
trainingSize <- ceiling(nr/2) 
#set seed
set.seed(201812)
train <- sample(1:nr, trainingSize)
#generate train set
X_train <- X[train,]
Crude_Rate_train <- Crude_Rate[train]
glimpse(X_train)
glimpse(Crude_Rate_train)

#Build regression,alpha = 1,type lasson
USPRCmodel <- glmnet(x = X_train, 
                     y = Crude_Rate_train, 
                     alpha = 1, 
                     family = "gaussian")

plot(USPRCmodel, las = 1)

####
set.seed(201812)
# k=10
cv_USPRCmodel <- cv.glmnet(X_train, 
                           Crude_Rate_train, 
                           alpha = 1,
                           family = "gaussian")
# log(lambda),MSE
plot(cv_USPRCmodel)

###
cbind(lambda = round(cv_USPRCmodel$lambda, 2), 
      cvm = signif(cv_USPRCmodel$cvm, 7), 
      cvsd = signif(cv_USPRCmodel$cvsd, 5))

#cvm
cvm <- cv_USPRCmodel$cvm
#Find mix
sub1 <- which.min(cvm); sub1

# cvsd
# cvmTest
cvmTest <- cvm[sub1] + cv_USPRCmodel$cvsd[sub1]; cvmTest

# find all cvm < cvmTest, take out the first one
sub2 <- which(cvm < cvmTest)[1]; sub2

# lambda value reflect on sub2
cv_USPRCmodel$lambda[sub2]

# lambda.1se -- largest value of lambda such that 
# error is within 1 standard error of the minimum.
cv_USPRCmodel$lambda.1se 

#min lambda
lambda.min <- cv_USPRCmodel$lambda.min; lambda.min

lambda.1se <- cv_USPRCmodel$lambda.1se; lambda.1se


#usd lambda.min predict 
lasso.predBest <- predict(USPRCmodel, s = lambda.min, newx = X_test)
#compare by lambda.1se
lasso.pred1se <- predict(USPRCmodel, s = lambda.1se, newx = X_test)
# compare MSE
mean((lasso.predBest - Crude_test)^2)

mean((lasso.pred1se - Crude_test)^2)

# view lambda.min
lasso.coef <- predict(USPRCmodel, type = "coefficients", s = lambda.min)
# lasso.coef to matrix
coefMat <- as.matrix(lasso.coef)
# remove zero
coefMat <- cbind(coefMat, coefMat != 0)
colnames(coefMat) <- c("Coefficent", "Keep")
round(coefMat[coefMat[, 2] == 1, ], 3)

# view lambda.1se
lasso.coef.1se <- predict(USPRCmodel, type = "coefficients", s = lambda.1se)
coefMat.1se <- as.matrix(lasso.coef.1se)
coefOther <- round(coefMat.1se[coefMat.1se != 0, ], 3); coefOther


#############################################################################

#Part Two
require(micromap)
USPRC <- readr::read_csv("G:/UnitedStatesandPuertoRicoCancerStatistics.csv", col_names = TRUE)
USPRC <- USPRC %>%
  filter(!Population %in% "Not Applicable") %>%
  mutate(Population = as.double(Population))
USPRC$`Leading Cancer Sites` %>%
  unique()

Ustates <- USPRC$States %>%
  unique()

# df_cancer
df_cancer <- as.data.frame(matrix(NA, nrow = length(Ustates), ncol = 3))

colnames(df_cancer) <- c("states", "Brain_and_Other_Nervous_System", "Stomach")

for (i in 1:length(Ustates)) {
  stateName <- Ustates[i]
  df_cancer[i, "states"] <- stateName
  BrainCancer <- USPRC %>%
    filter(
      States %in% stateName,
      `Leading Cancer Sites` %in% "Brain and Other Nervous System"
    )
  df_cancer[i, "Brain_and_Other_Nervous_System"] <-
    (sum(BrainCancer$Count) / sum(BrainCancer$Population)) %>%
    `*`(100000) %>%
    round(2)
  StomachCancer <- USPRC %>%
    filter(
      States %in% stateName,
      `Leading Cancer Sites` %in% "Stomach"
    )
  df_cancer[i, "Stomach"] <-
    (sum(StomachCancer$Count) / sum(StomachCancer$Population)) %>%
    `*`(100000) %>%
    round(2)
}

df_cancer


print(USPRC %>%
        filter(
          States %in% "District of Columbia",
          `Leading Cancer Sites` %in% "Brain and Other Nervous System"
        ))


df_cancer[51, "Brain_and_Other_Nervous_System"] <- 0


df_cancer[51, "states"] <- "Washington D.C."

stateRegion <- readr::read_csv("G:/stateRegion.csv", col_names = TRUE)

df_cancer$StateAb <- NA
df_cancer <- df_cancer[match(sort(stateRegion$state), df_cancer$states), ]
df_cancer$states <- as.factor(df_cancer$states)
df_cancer$StateAb <- stateRegion$StateAb
df_cancer <- df_cancer[, c(1, 4, 2, 3)]
head(df_cancer)

statePolys <- readr::read_csv("g:/statePolys.csv", col_names = TRUE) %>% 
  as.data.frame()
head(statePolys)








