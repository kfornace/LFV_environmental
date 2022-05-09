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

#null model

m0 <- glmer(lfv ~ + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
               nAGQ = 10)
summary(m0)

#of5000fd

m1 <- glmer(lfv ~ + of5000fd + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m1)
anova(m0, m1)

#of20000pa
m2 <- glmer(lfv ~ + of5000fd + of20000pa + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m2)
anova(m1, m2)

#shr5000pa

m3 <- glmer(lfv ~ + of5000fd + of20000pa + shr5000pa + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m3)
anova(m2, m3)

#bt20000
m4 <- glmer(lfv ~ + of5000fd + of20000pa + bt20000 +(1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m4)
anova(m2, m4)

#bt5000sh
m5 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m5)
anova(m2,m5)

#vg5000pa
m6 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + vg5000pa + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m6)
anova(m5, m6)

#bt20000sh ### already have a bt sh w lower aic ## look anyway?
m7 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + bt20000sh + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m7)
anova(m5,m7)
#bt5000pa
m8 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + bt5000pa + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m8)
anova(m5,m8)

#shr2000pa
m9 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + shr2000pa + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m9)
anova(m5,m9)

#cf10000pa

m10 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + shr2000pa + cf10000pa + (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m10)
anova(m9, m10)

#vg500pa
m11 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + shr2000pa + cf10000pa + vg500pa + 
               (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m11)
anova(m10, m11)

#of2000fd ## already have an of fd lower aic ## look anyway
m12 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + shr2000pa + cf10000pa + of2000fd +
               (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m12)
anova(m10,m12)

#vg20000pa
m13 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + shr2000pa + cf10000pa + vg20000pa +
               (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m13)
anova(m10, m13)

#cf500pa
m14 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + shr2000pa + cf10000pa + vg20000pa + cf500pa +
               (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m14)
anova(m13, m14)


#bt2000sh
m15 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + shr2000pa + cf10000pa + vg20000pa + bt2000sh +
               (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m15)
anova(m13,m15)

#shr2000

m16 <- glmer(lfv ~ + of5000fd + of20000pa + bt5000sh + shr2000pa + cf10000pa + vg20000pa + shr2000 +
               (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)

m17 <- glmer(lfv ~ + of5000fd + bt5000sh + shr2000pa + cf10000pa + vg20000pa + shr2000 +
               (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)

summary(m16)
anova(m13, m16)

msigman <- m16 <- glmer(lfv ~ + of5000fd  + cf10000pa + vg20000pa +
                          (1 | village), data = macenda_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
                        nAGQ = 10)

summary(m16)
summary(msigman)
cc <- confint(m16,parm = "beta_", method="Wald")
ctab <- cbind(est=fixef(m16),cc)
rtab_mac <- exp(ctab)
rtab_mac

write.csv(rtab, "mac_OR.csv")

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


m0 <- glmer(lfv ~ + ageCat +
                     (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
                   nAGQ = 10)
summary(m0)

#bt2000pa

m1 <- glmer(lfv ~ + ageCat + bt2000pa +
                  (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
                nAGQ = 10)
summary(m1)
anova(m0,m1)

#cf500pa


m2 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa+
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m2)
anova(m1,m2)


#shr5000fd
m3 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + shr5000fd +
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m3)
anova(m2,m3)

#cf2000
m4 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa  + cf2000 +
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m4)
anova(m2, m4)

#of10000
m5 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + of10000 +
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m5)
anova(m2, m5)

#bt10000
m6 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + of10000 + bt10000 +
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m6)
anova(m5, m6)

#vg1000sh
m7 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa  + of10000 + vg1000sh+
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m7)
anova(m5, m7)

#bt5000sh
m8 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa  + of10000 + bt5000sh +
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)

summary(m8)
anova(m5, m8)


#vg500pa

m9 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + of10000 + bt5000sh + vg500pa +
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m9)
anova(m8, m9)

#vg1000fd
m10 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + of10000 + bt5000sh + vg1000fd +
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m10)
anova(m8, m10)

#of500sh

m11 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + of10000 + bt5000sh + of500sh +
              (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
summary(m11)
anova(m8, m11)


#cf5000pa ### already an cf pa
m12 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + of10000 + bt5000sh + cf5000pa +
               (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m12)
anova(m8, m12)

#vg20000pa

m13 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa  + of10000 + bt5000sh + vg20000pa +
               (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m13)
anova(m8, m13)

#vg5000fd ### already an fd
m14 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa  + of10000 + bt5000sh + vg20000pa + vg5000fd +
               (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m14)
anova(m13, m14)

#vg10000
m15 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + of10000 + bt5000sh + vg20000pa + vg10000 +
               (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m15)
anova(m13, m15)

#bt10000pa

m16 <- glmer(lfv ~ + ageCat + bt2000pa + cf500pa + of10000 + bt5000sh + vg20000pa + bt10000pa +
               (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
summary(m16)
anova(m13, m16)

msig <- glmer(lfv ~ + ageCat + bt2000pa + of10000  +
               (1 | village), data = coastal_sig, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)


summary(msig)
summary(m13)

cc <- confint(m13,parm = "beta_", method="Wald")
ctab <- cbind(est=fixef(m13),cc)
rtab_bg <- exp(ctab)
rtab_bg

write.csv(rtab, "bg_OR.csv")

