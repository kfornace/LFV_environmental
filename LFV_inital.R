library(lme4)
library(plyr)
library(corrplot)

## Read data
setwd()
df <- read.csv("Lassa_merged.csv")
ev <- read.csv("All_LFV_environmental.csv")
ev[5:131] <- scale(ev[5:131])
## Merge by village
df <- merge(df, ev, by.x="village", by.y="vill.name")

## Set outcome variable 
df$result <- factor(df$result)
df$lfv <- ifelse(df$result=="NEG", 0, 1)

## Village as group
df$village <- factor(df$village)
df$ageCat <- factor(findInterval(df$age, c(15, 26, 37, 48, 59)))


###################################################### Macenta + Gueckedou ##########################################################################

## Subset data 
df2 <- subset(df, study>1)

## Null model 
m0 <- glmer(lfv ~ + (1 | village), data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m0)

## Including age as continuous variable
m1 <- glmer(lfv ~ age + (1 | village), data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m1)
anova(m0,m1)

## Including age as a categorical value
m2 <- glmer(lfv ~ ageCat + (1 | village), data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m2)
anova(m1,m2)

## Including gender - exclude from all data model 
m3 <- glmer(lfv ~ ageCat + sex + (1 | village), data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m3)
anova(m2,m3)

## Correlation for all environmental variables
ev <- df2[15:140]
corrplot(cor(ev, method="spearman"), method="color")

## Run model for all environmental variables
vars <- colnames(ev)
for(i in 1:length(vars)){
  
  # Set formula
  fmla <- paste0("lfv ~ + ", vars[i], " + (1 | village)")
  
  # Run model 
  m1 <- glmer(fmla, data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
              nAGQ = 10)
  
  # Extract results
  se <- sqrt(diag(vcov(m1)))
  cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
  cof <- data.frame(exp(cof))
  a <- anova(m0, m1)
  pv <- a$`Pr(>Chisq)`[2]
  cof$pvalue <- pv
  cof$aic <- AIC(m1)
  cof$bic <- BIC(m1)
  
  # Name variables and remove intercept
  cof$variable <- row.names(cof)
  row.names(cof) <- NULL
  cof <- subset(cof, variable !="(Intercept)")
  
  # Add to main dataset
  # If merged dataset does not exist, create dataset
  if(!exists("dataset")){
    dataset <- cof
  }
  
  # If merged dataset exists, append data
  if(exists("dataset")){
    dataset <- rbind(dataset, cof)
    rm(cof)
  }
}

## Write results
write.csv(dataset, "LFV_initial_Macenta_Guek.csv")
rm(dataset)

###################################################### Coastal ##########################################################################

## Subset data 
df2 <- subset(df, study<2)

## Null model 
m0 <- glmer(lfv ~ + (1 | village), data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m0)

## Including age as continuous variable
m1 <- glmer(lfv ~ age + (1 | village), data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m1)
anova(m0,m1)

## Including age as a categorical value
m2 <- glmer(lfv ~ ageCat + (1 | village), data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m2)
anova(m1,m2)

## Including gender - exclude from all data model 
m3 <- glmer(lfv ~ age + sex + (1 | village), data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m3)
anova(m2,m3)

m0 <- m2

## Correlation for all environmental variables
ev <- df2[15:140]
corrplot(cor(ev, method="spearman"), method="color")

## Run model for all environmental variables
vars <- colnames(ev)
for(i in 1:length(vars)){
  
  # Set formula
  fmla <- paste0("lfv ~ + ageCat +", vars[i], " + (1 | village)")
  
  # Run model 
  m1 <- glmer(fmla, data = df2, family = binomial, control = glmerControl(optimizer = "bobyqa"),
              nAGQ = 10)
  
  # Extract results
  se <- sqrt(diag(vcov(m1)))
  cof <- cbind(Est = fixef(m1), LL = fixef(m1) - 1.96 * se, UL = fixef(m1) + 1.96 *se)
  cof <- data.frame(exp(cof))
  a <- anova(m0, m1)
  pv <- a$`Pr(>Chisq)`[2]
  cof$pvalue <- pv
  cof$aic <- AIC(m1)
  cof$bic <- BIC(m1)
  
  # Name variables and remove intercept
  cof$variable <- row.names(cof)
  row.names(cof) <- NULL
  cof <- subset(cof, variable !="(Intercept)")
  
  # Add to main dataset
  # If merged dataset does not exist, create dataset
  if(!exists("dataset")){
    dataset <- cof
  }
  
  # If merged dataset exists, append data
  if(exists("dataset")){
    dataset <- rbind(dataset, cof)
    rm(cof)
  }
}

## Write results
dataset <- subset(dataset, variable != "ageCat1")
dataset <- subset(dataset, variable != "ageCat2")
dataset <- subset(dataset, variable != "ageCat3")
dataset <- subset(dataset, variable != "ageCat4")
dataset <- subset(dataset, variable != "ageCat5")

write.csv(dataset, "LFV_initial_Coastal.csv")
rm(dataset)