library(stringr)
library(caret)

df <- read.csv("Lassa_merged.csv")
ev <- read.csv("All_LFV_environmental.csv")
ev[5:131] <- scale(ev[5:131])

## Merge by village
df <- merge(df, ev, by.x="village", by.y="vill.name")
df$ageCat <- factor(findInterval(df$age, c(15, 26, 37, 48, 59)))

age <- df %>% select(ageCat, age)

## Set outcome variable 
df$result <- factor(df$result)
df$lfv <- ifelse(df$result=="NEG", 0, 1)

## Village as group
df$village <- factor(df$village)

macenda <- subset(df, study > 1)
coastal <- subset(df, study < 2)



###################################################### Macenta + Gueckedou ##########################################################################

mac_set <- read.csv("LFV_initial_Macenta_Guek.csv")
mac_set_sig <- subset(mac_set, pvalue<0.2)
mac_sig_vars <- mac_set_sig$variable

macenda_sig <- macenda %>% select("lfv", "village", one_of(mac_sig_vars))


## Check for correlation between selected variables
correlationMatrix_mac <- cor(macenda_sig[mac_sig_vars])

## summarize the correlation matrix
print(correlationMatrix_mac)

## Find attributes that are highly corrected (ideally >0.75)
hc_mac <- findCorrelation(correlationMatrix_mac, cutoff=0.75)

## Print indexes of highly correlated attributes
print(hc_mac)

## Remove highly correlated variables
hc_mac = sort(hc_mac)
macenda2 = macenda[mac_sig_vars]
macenda2 = macenda2[,-c(hc_mac)]
head(macenda2)
vars_mac <- names(macenda2)

mac_set_nocorr_sig <- mac_set_sig[mac_set_sig$variable %in% vars_mac,]

mac_set_nocorr_sig %>% arrange(aic) 


#model selection forested


vars_mac<- c("lfv", "village", vars_mac)

mac_vars <- names(macenda_sig)[(names(macenda_sig) %in% vars_mac)]

var <- paste(names(mac_vars[3:15]), collapse="+")

frmla <- as.formula(paste("lfv ~ +" , var, " + (1 | village)"))

mac_vars <- macenda_sig[, mac_vars]

mixed_model <- glmer(frmla, data = mac_vars,  family = binomial, control = glmerControl(optimizer = "bobyqa"),
                     nAGQ = 10)


options(na.action = "na.fail") # Required for dredge to run
mac_dredge <- dredge(mixed_model, beta = "none", evaluate = T, rank = AICc)

## best model

top_model <- get.models(mac_dredge, subset = 1)[[1]]
summary(top_model)

## averaged model

summary(model.avg(mac_dredge, subset = delta <= 2))
avg_model_mac<- model.avg(mac_dredge, subset = delta <= 2)



write.csv(rtab_av_mac, "mac_OR.csv")

###################################################### Coastal ##########################################################################

bg_set <- read.csv("LFV_initial_coastal.csv")
bg_set_sig <- subset(bg_set, pvalue<0.2)
bg_sig_vars <- bg_set_sig$variable

coastal_sig <- coastal %>% select("lfv", "ageCat", "village", one_of(bg_sig_vars))


## Check for correlation between selected variables
correlationMatrix_bg <- cor(coastal_sig[bg_sig_vars])

## summarize the correlation matrix
print(correlationMatrix_bg)

## Find attributes that are highly corrected (ideally >0.75)
hc_bg <- findCorrelation(correlationMatrix_bg, cutoff=0.75)

## Print indexes of highly correlated attributes
print(hc_bg)

## Remove highly correlated variables
hc_bg = sort(hc_bg)
coastal2 = coastal[bg_sig_vars]
coastal2 = coastal2[,-c(hc_bg)]
head(coastal2)
vars_bg <- names(coastal2)

bg_set_nocorr_sig <- bg_set_sig[bg_set_sig$variable %in% vars_bg,]

bg_set_nocorr_sig %>% arrange(aic) 

### model selection coastal

vars_bg <- c("lfv", "village", "ageCat", vars_bg)

bg_vars <- coastal_sig[, vars_bg]

var <- paste(names(bg_vars[3:15]), collapse="+")

frmla <- as.formula(paste("lfv ~ " , var, " + (1 | village)"))

bg_vars <- coastal_sig[, bg_vars]

mixed_model <- glmer(frmla, data = bg_vars,  family = binomial, control = glmerControl(optimizer = "bobyqa"),
                     nAGQ = 10)


options(na.action = "na.fail") # Required for dredge to run
bc_dredge <- dredge(mixed_model, beta = "none", evaluate = T, rank = AICc)


## best model

top_model_bc <- get.models(bc_dredge, subset = 1)[[1]]
summary(top_model_bc)

cc2 <- confint(top_model_bc,parm = "beta_", method="Wald")
ctab2 <- cbind(est=fixef(top_model_bc),cc2)
rtab_bg2 <- exp(ctab2)
rtab_bg2

## averaged model

summary(model.avg(bc_dredge, subset = delta <= 2))
avg_model_bc <- model.avg(bc_dredge, subset = delta <= 2)
summary(avg_model_bc)


write.csv(rtab_av_bc, "bg_OR.csv")