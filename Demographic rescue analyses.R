#==========================================================================================================#
# Script created by Catherine Searle, contact at searlester@gmail.com, searlec@purdue.edu
# Script created in R version 4.3.1
# This Script: Analysis of demographic rescue data
# Usage Notes: Run one line at a time
#==========================================================================================================#

# Set working directory, import packages, load files, initialize global variables, source functions

#======================================================================================================================#

# load data:
rawdat <- read.table("Demographic_rescue_5.txt", header=TRUE, sep="\t", na.strings="", dec=".", strip.white=TRUE) 
dat<- rawdat[which(rawdat$beaker.id!="B39" & rawdat$beaker.id!="B36"),] # excluding 2 beakers that were accidentally exposed 
# infected only beakers:
dat.inf = dat[dat$metsch.treat == 'INF', ]

library(glmmTMB)
library(bbmle)
library(lme4)

#####################################################
# POPULATION DENSITY
# main model

densmixmodel<-glmer.nb(dat$tot.num~
                 dat$metsch.treat.id+dat$sup.treat.id#+dat$week 
                #+dat$metsch.treat.id:dat$sup.treat.id 
                #+dat$block
                +(1|dat$week)
                ) 

summary(densmixmodel)
dispersion_glmer(densmixmodel) 

densmixmodelB<-glmer.nb(dat$tot.num~
                         dat$metsch.treat.id+dat$sup.treat.id+dat$week
                       +dat$metsch.treat.id:dat$sup.treat.id 
                       +(1|dat$beaker.id)
) 

summary(densmixmodelB)
dispersion_glmer(densmixmodelB) 


# Pairwise comparisons 
densmixmodel2<-glmer.nb(dat$tot.num~
                       #dat$metsch.treat.id # metsch removed
                       +dat$sup.treat.id
                       +dat$week 
                       +dat$metsch.treat.id:dat$sup.treat.id 
                       +dat$block
                       +(1|dat$beaker.id)) 

densmixmodel3<-glmer.nb(dat$tot.num~
                          dat$metsch.treat.id
                        #+dat$sup.treat.id # supplemental treatment removed
                        +dat$week 
                        +dat$metsch.treat.id:dat$sup.treat.id 
                        +dat$block
                        +(1|dat$beaker.id)) 

densmixmodel4<-glmer.nb(dat$tot.num~
                          dat$metsch.treat.id
                        +dat$sup.treat.id
                        +dat$week 
                        #+dat$metsch.treat.id:dat$sup.treat.id  # interaction removed 
                        +dat$block
                        +(1|dat$beaker.id)) 

AIC(densmixmodel, densmixmodel2, densmixmodel3, densmixmodel4, densmixmodelB)

densmixmodel5<-glmer.nb(dat$tot.num~
                        dat$metsch.treat.id
                        +dat$sup.treat.id
                        +dat$week 
                        #+dat$metsch.treat.id:dat$sup.treat.id  # interaction removed 
                        #+dat$block
                        +(1|dat$beaker.id),
                        data=dat) 

AIC(densmixmodel5) 

densmixmodel5.1<-glmer.nb(dat$tot.num~
                         # dat$metsch.treat.id
                        +dat$sup.treat.id
                        +dat$week 
                        +(1|dat$beaker.id),
                        data=dat) 
densmixmodel5.2<-glmer.nb(dat$tot.num~
                            dat$metsch.treat.id
                          #+dat$sup.treat.id
                          +dat$week 
                          +(1|dat$beaker.id),
                          data=dat) 
densmixmodel5.3<-glmer.nb(dat$tot.num~
                            dat$metsch.treat.id
                          +dat$sup.treat.id
                          #+dat$week 
                          +(1|dat$beaker.id),
                          data=dat) 

anova(densmixmodel5,densmixmodel5.1) # metsch
anova(densmixmodel5,densmixmodel5.2) # supplement
anova(densmixmodel5,densmixmodel5.3) # week


###########################
# Each week separately
# week 1
dat.1 = dat[dat$week == '1', ]
model.1<-glm(dat.1$tot.num~
             dat.1$metsch.treat.id+dat.1$sup.treat.id
             ,family=poisson(link="log"))
drop1(model.1,~., test="Chisq") 

# week 2
dat.2 = dat[dat$week == '2', ]
model.2<-glm(dat.2$tot.num~
              dat.2$metsch.treat.id+dat.2$sup.treat.id
             ,family=poisson(link="log"))
drop1(model.2,~., test="Chisq") 

# week 3
dat.3 = dat[dat$week == '3', ]
model.3<-glm(dat.3$tot.num~
              dat.3$metsch.treat.id+dat.3$sup.treat.id
             ,family=poisson(link="log"))
drop1(model.3,~., test="Chisq") 

# week 4
dat.4 = dat[dat$week == '4', ]
model.4<-glm(dat.4$tot.num~
              dat.4$metsch.treat.id+dat.4$sup.treat.id
             ,family=poisson(link="log"))
drop1(model.4,~., test="Chisq") 

# week 5
dat.5 = dat[dat$week == '5', ]
model.5<-glm(dat.5$tot.num~
              dat.5$metsch.treat.id+dat.5$sup.treat.id
             ,family=poisson(link="log"))
drop1(model.5,~., test="Chisq") 

# week 6
dat.6 = dat[dat$week == '6', ]
model.6<-glm(dat.6$tot.num~
              dat.6$metsch.treat.id+dat.6$sup.treat.id,
             family=poisson(link="log"))
drop1(model.6,~., test="Chisq") 

# week 7
dat.7 = dat[dat$week == '7', ]
model.7<-glm(dat.7$tot.num~
              dat.7$metsch.treat.id+dat.7$sup.treat.id,
             family=poisson(link="log"))
drop1(model.7,~., test="Chisq") 

# week 8
dat.8 = dat[dat$week == '8', ]
model.8<-glm(dat.8$tot.num~
              dat.8$metsch.treat.id+dat.8$sup.treat.id,
             family=poisson(link="log"))
drop1(model.8,~., test="Chisq") 


###############
# NUMBER INFECTED 

infdens<-glmer(dat.inf$tot.num.inf~dat.inf$sup.treat.id+dat.inf$week # main effects
                +dat.inf$block
                +(1|dat.inf$beaker.id), 
                family=poisson (link = "log"))
summary(infdens) 
dispersion_glmer(infdens) 
drop1(infdens,~., test="Chisq") 

# model selection: 
infdens1<-glmer(dat.inf$tot.num.inf~#dat.inf$sup.treat.id 
               +dat.inf$week
               +dat.inf$block
               +(1|dat.inf$beaker.id), 
               family=poisson (link = "log"))
infdens2<-glmer(dat.inf$tot.num.inf~dat.inf$sup.treat.id
               #+dat.inf$week # drop week
               +dat.inf$block
               +(1|dat.inf$beaker.id), 
               family=poisson (link = "log"))
infdens3<-glmer(dat.inf$tot.num.inf~dat.inf$sup.treat.id
               +dat.inf$week
               #+dat.inf$block # drop block 
               +(1|dat.inf$beaker.id), 
               family=poisson (link = "log"))
AIC(infdens1,infdens2,infdens3) 

# pairwise comparisons for supplementation: 
dat.inf2A<- dat.inf[which(dat.inf$sup.treat.id!="0"),] # dropping control
dat.inf2B<- dat.inf[which(dat.inf$sup.treat.id!="1"),] # dropping low
dat.inf2C<- dat.inf[which(dat.inf$sup.treat.id!="2"),] # dropping high

infdens2A<-glmer(dat.inf2A$tot.num.inf~dat.inf2A$sup.treat.id+dat.inf2A$week # low versus high
               +dat.inf2A$block
               +(1|dat.inf2A$beaker.id), 
               family=poisson (link = "log"))
drop1(infdens2A,~., test="Chisq") 

infdens2B<-glmer(dat.inf2B$tot.num.inf~dat.inf2B$sup.treat.id+dat.inf2B$week # control vs. high
                 +dat.inf2B$block
                 +(1|dat.inf2B$beaker.id), 
                 family=poisson (link = "log"))
drop1(infdens2B,~., test="Chisq")

infdens2C<-glmer(dat.inf2C$tot.num.inf~dat.inf2C$sup.treat.id+dat.inf2C$week # control vs. low
                 +dat.inf2C$block
                 +(1|dat.inf2C$beaker.id), 
                 family=poisson (link = "log"))
drop1(infdens2C,~., test="Chisq")


###############
# INFECTION PREVALENCE
infprev2<-glmer(dat.inf$prop.inf
               ~dat.inf$sup.treat.id+dat.inf$week
               +dat.inf$sup.treat.id:dat.inf$week
               +dat.inf$block
               +(1|dat.inf$beaker.id), 
               family=binomial ) 
summary(infprev2) 
dispersion_glmer(infprev2)
drop1(infprev2,~., test="Chisq") 

###############
# PROPORTION MALES 

maleprev3<-glmmTMB(dat$prop.male
                 ~dat$sup.treat.id+dat$metsch.treat.id#+dat$week 
                 #+dat$sup.treat.id:dat$metsch.treat.id
                 #+dat$block
                 +(1|dat$beaker.id),
                 data=dat,
                 family=nbinom1) 
summary(maleprev3) 
anova(maleprev3)
drop1(maleprev3,~., test="Chisq")


###############
# PROPORTION JUVENILES
propjuv<-glmer(dat$prop.juv
                ~dat$sup.treat.id+dat$metsch.treat.id+dat$week 
                +dat$sup.treat.id:dat$metsch.treat.id
                +dat$block
                +(1|dat$beaker.id), 
                family=binomial ) 
summary(propjuv) 
dispersion_glmer(propjuv) 
drop1(propjuv,~., test="Chisq") 


