#####WOOD DECOMPOSITION ACROSS A Distrubance gradient in tropical montane forest in Mengsong for 3 years southwest China#####
###Clean the R's brain
#rm(list=ls())

#setwd("")

###Load libraries in need
library(ggplot2)
library(lme4)
library(car)
library(arm)
library(MuMIn)
library(nlme)
library(rms)
library(AICcmodavg)
library(LMERConvenienceFunctions)
library(plyr)
library(dplyr)
library(tidyr)

##Import data 
Mengsong_field <- read.csv("Mengsong_field_data_FEM.csv")


###Import Chemsitry data

ini_chemistry<-read.csv(file="Mengsong_samples_chemical_analysis.csv", header=T, sep=",")


###Wood logs chemical content 
###Piping
chemistry_logs<-ini_chemistry%>%
  group_by(Species,Tissue)%>%
  summarize(mean_C=mean(Carb),
            mean_N=mean(Nitro),
            mean_P=mean(Phosp),
            mean_K=mean(Pota),
            mean_sugar=mean(Sugar),
            mean_Lig=mean(Lignin),
            mean_Cell=mean(Cellu),
            mean_Hemi=mean(Hemicellu),
            mean_Fiber=mean(Fiber_content),
            mean_Tannin=mean(Tannin),
            mean_C_N_ratio=mean(Carb_Nitro_ratio),
            mean_Lignin_N_ratio=mean(Lignin_Nitro_ratio),
            mean_Lignin_C_ratio=mean(Lignin_Carb_ratio),
  )

chemistry_logs
###Table 1

###Wood logs characteristics

###Compute mean length of logs
###Piping
logs_initial_characteristics<-Mengsong_field%>%
  group_by(Species)%>%
  summarize(mean_length=mean(L),
            sd_length=sd(L),
            se_length=sd(L)/sqrt(n()),
            mean_WSG_initial=mean(WSG_ini.y),
            sd_WSG_initial=sd(WSG_ini.y),
            se_WSG_initial=sd(WSG_ini.y)/sqrt(n()),
            mean_Log_weight_initial=mean(Log_mass_ini),
            sd_Log_weight_initial=sd(Log_mass_ini),
            se_Log_weight_initial=sd(Log_mass_ini)/sqrt(n()),
            mean_bark_thick_initial=mean(Bark_thickness_T1),
            sd_Log_bark_thick_initial=sd(Bark_thickness_T1),
            se_Log_bark_thick_initial=sd(Bark_thickness_T1)/sqrt(n())##n() counts rows in each group
  )
logs_initial_characteristics

logs_initial_characteristics<-as.data.frame(logs_initial_characteristics)
logs_initial_characteristics


######WOOD SPECIFIC GRAVITY LOSS ANALYSIS

####Modeling of wood specific gravity loss
max(Mengsong_field$Per_WSG_loss[Mengsong_field$Number_days>1100 & Mengsong_field$Species=="Lit_cub"])
min(Mengsong_field$Per_WSG_loss[Mengsong_field$Number_days>1100 & Mengsong_field$Species=="Lit_cub"])

max(Mengsong_field$Per_WSG_loss[Mengsong_field$Number_days>1100])

min(Mengsong_field$Per_WSG_loss[Mengsong_field$Number_days>1100])

par(mfrow=c(1,1))
plot(Mengsong_field$Number_days)

summary(Mengsong_field)
Mengsong_field$PLOT_S_PLOT<-paste(Mengsong_field$PLOT,Mengsong_field$S_PLOT, sep="_")
m1 <- lmer(log(Per_WSG_loss) ~  Number_days+Species+For_type+Pos_soil+Species:Pos_soil+
             (1|PLOT/S_PLOT),
           control=lmerControl(optCtrl=list(maxfun=20000)),data = Mengsong_field,na.action=na.omit)

#Diagnostics
resids <- resid(m1,type='pearson')
#Checking linearity
par(mfrow=c(1,1))
plot(resids~fitted(m1))
lines(lowess(resids~fitted(m1)),col='red') #slight curvature

#checking homoscedasticity
plot(resids~log(Mengsong_field$Per_WSG_loss))## no homoscedacity perhaps should consider lme
plot(sqrt(abs(resids))~fitted(m1))
lines(lowess(sqrt(abs(resids))~fitted(m1)),col='red') # slight curvature at higher and lower values

#checking normality at all levels
qqPlot(resids) # Over dispersion but a bit terrible
qqPlot(ranef(m1)$'PLOT'$'(Intercept)') #ok

##prediction
preddat <- expand.grid(Number_days = c(91,182,365,548,730,1095),
                       Species = c('Cast_mek','Lit_cub'),
                       For_type=c("OC","CC", "OL"),
                       Pos_soil=c('up','down'))

dim(preddat)

m1.pred <- predict(m1,newdata=preddat,re.form=~0)
#now for confindence intervals
m1.form <- formula(m1,fixed.only=T)
m1.form <- update(m1.form, NULL~.)
mod.mat <- model.matrix(m1.form,preddat)
vcv <- vcov(m1)
semod <- sqrt(diag(mod.mat%*%vcv%*%t(mod.mat)))

uc_all <- m1.pred + semod*1.96
lc_all <- m1.pred - semod*1.96

#convert to original scale
m1.pred <- exp(m1.pred)
uc.all <- exp(uc_all)
lc.all <- exp(lc_all)

##results
dcmp.results_Dossa <- list()
dcmp.results_Dossa$confint <- confint(m1,method = 'boot',nsim=999,parallel='snow',ncpus=4)
dcmp.results_Dossa$r.squared <- r.squaredGLMM(m1)
dcmp.results_Dossa$model.predictions <- cbind(preddat,m1.pred,uc.all,lc.all)


#######Using lme because of heteroscedacity 


m2 <- lme(fixed=log(Per_WSG_loss) ~  poly(Number_days,2)+Species+For_type+Pos_soil+poly(Number_days,2):Pos_soil+poly(Number_days,2):Species+Species:Pos_soil+poly(Number_days,2):For_type+For_type:Pos_soil,random=~1|PLOT_S_PLOT,data = Mengsong_field,na.action=na.omit)

m2 <- lme(fixed=log(Per_WSG_loss) ~  poly(Number_days,2)+Species+For_type+Pos_soil+poly(Number_days,2):Species+Species:Pos_soil,random=~1|PLOT/S_PLOT,data = Mengsong_field,na.action=na.omit)
summary(m2)
AIC(m2)

plot(m2, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance has a decreasing trend. 


##Dealing with the variance trend
## model the variance as a power function of the fitted values

m2a <- update(m2, weights=varPower(form=~fitted(.)))
plot(m2a, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## it is improved a lot

AIC(m2a)
## model the variance as a exponential function of the fitted values
m2b <- update(m2, weights=varExp(form=~fitted(.)))
plot(m2b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## much similarto previous model m2a
AIC(m2b)

## let's  make the variance be different 
m2c <- update(m2, weights=varIdent(form=~1|Species*Pos_soil))
plot(m2c, abs(resid(., type='pearson'))~fitted(.)|Pos_soil, type=c('p', 'r'), 
     abline=0) 
plot(m2c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## Became worst than previous
AIC(m2c)

AIC(m2, m2a,m2b,m2c)

summary(m2b)
anova(m2b, type = "marginal")
##Table

r.squaredLR(m2b)

intervals(m2b)
###Table 

#### Consider interaction

m3 <- lme(fixed=log(Per_WSG_loss) ~  poly(Number_days,2)+Species+For_type+Pos_soil+as.factor(Termi.assum)+poly(Number_days,2):Pos_soil+poly(Number_days,2):Species+Species:Pos_soil+poly(Number_days,2):For_type+For_type:Pos_soil+
            Species:as.factor(Termi.assum)+poly(Number_days,2):as.factor(Termi.assum)+Pos_soil:as.factor(Termi.assum)+For_type:as.factor(Termi.assum),random=~1|PLOT_S_PLOT,data = Mengsong_field,na.action=na.omit)

summary(m3)
AIC(m3)

plot(m3, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance has a decreasing trend. 


m4 <- lme(fixed=log(Per_WSG_loss) ~  poly(Number_days,2)+Species+For_type+Pos_soil+as.factor(Termi.assum)+poly(Number_days,2):Species+Species:Pos_soil+For_type:as.factor(Termi.assum),random=~1|PLOT/S_PLOT,data = Mengsong_field,na.action=na.omit)
summary(m4)
AIC(m4)

plot(m4, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance has a decreasing trend. 


##Dealing with the variance trend
## model the variance as a power function of the fitted values

m4a <- update(m4, weights=varPower(form=~fitted(.)))
summary(m4a)
plot(m4a, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## it is improved a lot

AIC(m4a)
## model the variance as a exponential function of the fitted values
m4b <- update(m4, weights=varExp(form=~fitted(.)))
summary(m4b)
plot(m4b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## much similar to previous model m2a
AIC(m4b)

## let's  make the variance be different 
m4c <- update(m4, weights=varIdent(form=~1|Species*Pos_soil))
summary(m4c)
plot(m4c, abs(resid(., type='pearson'))~fitted(.)|Pos_soil, type=c('p', 'r'), 
     abline=0) 
plot(m4c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## Became worst than previous
AIC(m4c)

AIC(m4, m4a,m4b,m4c)

summary(m4b)
anova(m4b, type = "marginal")


r.squaredLR(m4b)

intervals(m4b)
####Adding termite presence absence does not improve the model. just very little effect of interaction and at one level with forest type (open land)
###Thus we chose the model m2b that does not include termites

###Best model is m2b 

##prediction
preddat <- expand.grid(Species = c('Cast_mek','Lit_cub'),
                       Pos_soil=c('up','down'),
                       PLOT=unique(Mengsong_field$PLOT),
                       S_PLOT=unique(Mengsong_field$S_PLOT),
                       Number_days = c(91,182,365,548,730,1095))

dim(preddat)


#write.csv(preddat,file="preddat.csv")
###We added a new column about Forest type accordingly 
preddat <- read.csv("preddat.csv")
head(preddat)
dim(preddat)
str(preddat)
preddat$PLOT
summary(m2b)

###
preddat$pred <- predict(m2b, newdata=preddat, level=0)

head(preddat)
###Prediction interval 
## [-2] drops response from formula

Designmat <- model.matrix(formula(m2b)[-2], preddat)
predvar <- diag(Designmat %*% vcov(m2b) %*% t(Designmat)) 
preddat$SE <- sqrt(predvar) 
preddat$SE2 <- sqrt(predvar+m2b$sigma^2) ## sigma= residual 

preddat$SE_uc<-preddat$pred+1.96*preddat$SE
preddat$SE_lc<-preddat$pred-1.96*preddat$SE
head(preddat)
preddat$Sp_Pos = paste(preddat$Species, preddat$Pos_soil, sep="_")
preddat$Sp_Pos <- as.factor(preddat$Sp_Pos)
levels(preddat$Sp_Pos)
preddat$For_type_Sp_Pos = paste(preddat$For_type, preddat$Sp_Pos, sep="_")
preddat$For_type_Sp_Pos <- as.factor(preddat$For_type_Sp_Pos)
levels(preddat$For_type_Sp_Pos)

preddat$For_type_Sp_Pos = paste(preddat$For_type, preddat$Sp_Pos, sep="_")
preddat$For_type_Sp_Pos <- as.factor(preddat$For_type_Sp_Pos)
levels(preddat$For_type_Sp_Pos)

#convert to original scale
preddat_orig<-preddat[,1:6]
preddat_orig$pred <- exp(preddat$pred)
preddat_orig$SE_uc<- exp(preddat$SE_uc)
preddat_orig$SE_lc<- exp(preddat$SE_lc)
preddat_orig$Sp_Pos<-preddat$Sp_Pos
preddat_orig$For_type_Sp_Pos<-preddat$For_type_Sp_Pos
preddat_orig$Number_mo<-as.factor(preddat$Number_days)
levels(preddat_orig$Number_mo)
levels(preddat_orig$Number_mo)<-c("3 months", "6 months", "12 months", "18 months", "24 months", "36 months")
levels (preddat_orig$Number_mo)
levels(preddat_orig$For_type)
levels(preddat_orig$For_type) <- c("MATURE", "REGENERATING", "OPEN")


head(preddat_orig)

library(ggplot2)
library(gridExtra)
library(reshape)
library(grid)

str(preddat_orig)

####From here on we only consider plotting harvest times 6, 12, 18, 24, 36 months.
####But for presentting good visual we need to have a blank 30 months so the interval between sampling is 6 months
#### we use the information from 3 months, but everything concerning WSG is kept blank (NA)
preddat_orig$Number_mo<-(ifelse(preddat_orig$Number_mo=="3 months", "30 months", as.character(preddat_orig$Number_mo)))
preddat_orig$Number_mo<-as.factor(preddat_orig$Number_mo)

### Let create NA data for 30 months as we did not sample that time. 
preddat_orig$pred<-(ifelse(preddat_orig$Number_mo=="30 months",NA, preddat_orig$pred))
preddat_orig$SE_uc<-(ifelse(preddat_orig$Number_mo=="30 months",NA, preddat_orig$SE_uc))
preddat_orig$SE_lc<-(ifelse(preddat_orig$Number_mo=="30 months",NA, preddat_orig$SE_lc))

summary(preddat_orig)

##change the levels of factor Number_mo order
preddat_orig$Number_mo<-factor(preddat_orig$Number_mo, levels=c("6 months","12 months","18 months","24 months","30 months","36 months"))
levels(preddat_orig$Number_mo)


r.squaredLR(m2b)

######July 1st 2019

###Plotting WSG against time using an extra 30 mo as blank for helping visualization
###and using open and fill circles for respecitively up and down
preddat_orig
dim(preddat_orig)

##graphing
##Relevel Species as Litsea cubeba and Castanopsis mekongensis
levels(preddat_orig$Species)
levels(preddat_orig$Species) <- c("Castanopsis mekongensis","Litsea cubeba")
###Make sure Litsea is taken as reference level
preddat_orig$Species<-relevel(preddat_orig$Species, ref="Litsea cubeba")
levels(preddat_orig$Species)

levels(preddat_orig$For_type)
levels(preddat_orig$For_type) <- c("Mature forest", "Regenerating forest", "Open land")

levels(preddat_orig$Pos_soil)
###Make sure up is taken as reference level for position
preddat_orig$Pos_soil<-relevel(preddat_orig$Pos_soil, ref="up")
levels(preddat_orig$Pos_soil)



Lit_cast_up_down_no_3mo<-preddat_orig

Lit_cast_up_down_no_3mo_graph<-ggplot(data=Lit_cast_up_down_no_3mo, aes(y=pred, x=Number_mo, colour= For_type,shape=Species, fill=Pos_soil))+
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  geom_rect(data=NULL,aes(xmin=3.5,xmax=4.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  geom_rect(data=NULL,aes(xmin=5.5,xmax=6.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  scale_shape_manual(values=c(24,21)) +
  scale_fill_manual(values=c(NA, "black"),guide=guide_legend(override.aes=list(shape=21)))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  geom_errorbar(aes(ymax = Lit_cast_up_down_no_3mo$SE_uc, ymin = Lit_cast_up_down_no_3mo$SE_lc), width=0, size = 0.8, position=position_dodge(1))+
  scale_y_continuous(limits=c(10,55))+
  scale_x_discrete(limits=levels(Lit_cast_up_down_no_3mo$Number_mo))+
  xlab(" Incubation duration") +
  ylab("Percentage wood specific gravity loss") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        #legend.position = c("bottom"),
        legend.position = c(0.1,0.7),
        legend.text=element_text(size=12,face="italic"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

Lit_cast_up_down_no_3mo_graph


###Figure 2 WSG dynamics 
tiff(filename="Figure 2 with shapes and symboles for up down with 30 mo 2.tiff", res = 400, width=5700, height=3000, compression = "lzw")
Lit_cast_up_down_no_3mo_graph
dev.off()


####Modeling with environmental factors (soil and vegetation)
##Merge Environmental factors with Mengsong_field_data_full_analysis1_Rhett
## input data
Envi_factors_Mengsong <- read.csv("Environmental factors_Mengsong.csv")
head(Envi_factors_Mengsong)
Mengsong_field
head(Mengsong_field)
dim(Envi_factors_Mengsong)
dim(Mengsong_field)
Mengsong_field_envir<-merge(Mengsong_field,Envi_factors_Mengsong, by.x="PLOT_S_PLOT", by.y= "PLOT_SBPLOT", all=TRUE)
dim(Mengsong_field_envir)
head(Mengsong_field_envir)
str(Mengsong_field_envir)

Mengsong_field_envir$PLOT<-Mengsong_field_envir$PLOT.x


##Start with soil
##remove Forest type and replaceby soil PC1 &2 and put SOIL_PC1, SOIL_PC2 and their respective interactions.
m2 <- lme(fixed=log(Per_WSG_loss) ~  poly(Number_days,2)+Species+For_type+Pos_soil+poly(Number_days,2):Species+Species:Pos_soil,random=~1|PLOT/S_PLOT,data = Mengsong_field,na.action=na.omit)


soil_m2<- lme(fixed=log(Per_WSG_loss) ~  poly(Number_days,2)+Species+SOIL_PC1+SOIL_PC2+Pos_soil+poly(Number_days,2):Species+Species:Pos_soil+poly(Number_days,2):SOIL_PC1+poly(Number_days,2):SOIL_PC2,random=~1|PLOT/S_PLOT,data = Mengsong_field_envir,na.action=na.omit)

##remove insignificant variables interaction soil PC and Pos_soil
summary(soil_m2)

plot(soil_m2, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance has a decreasing trend. 


##Dealing with the variance trend
## model the variance as a power function of the fitted values

soil_m2a <- update(soil_m2, weights=varPower(form=~fitted(.)))
plot(soil_m2a, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## it is improved a lot

AIC(soil_m2a)

## model the variance as a exponential function of the fitted values
soil_m2b <- update(soil_m2, weights=varExp(form=~fitted(.)))
plot(soil_m2b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## much similarto previous model m2a
AIC(soil_m2b)

## let's  make the variance be different 
soil_m2c <- update(soil_m2, weights=varIdent(form=~1|Species*Pos_soil))
plot(soil_m2c, abs(resid(., type='pearson'))~fitted(.)|Pos_soil, type=c('p', 'r'), 
     abline=0) 
plot(soil_m2c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## Became worst than previous
AIC(soil_m2c)

AIC(soil_m2, soil_m2a,soil_m2b,soil_m2c)

summary(soil_m2b)
anova(soil_m2b, type = "marginal")

r.squaredGLMM(m2b)###
r.squaredLR(soil_m2b)

intervals(soil_m2b)
##Remove Poly_days:PC2
soil_m3<-lme(fixed=log(Per_WSG_loss) ~  poly(Number_days,2)+Species+SOIL_PC1+SOIL_PC2+Pos_soil+poly(Number_days,2):Species+Species:Pos_soil+poly(Number_days,2):SOIL_PC1,
             random=~1|PLOT/S_PLOT,data = Mengsong_field_envir,na.action=na.omit,weights=varExp(form=~fitted(.)))
#str(Mengsong_field_envir)
summary(soil_m3)
anova(soil_m3)
AIC(soil_m2b,soil_m3)

##Now modeleing for vegetation

veg_m2<- lme(fixed=log(Per_WSG_loss) ~  poly(Number_days,2)+Species+veg_MDS1+veg_MDS2+Pos_soil+poly(Number_days,2):Species+Species:Pos_soil+poly(Number_days,2):veg_MDS1+poly(Number_days,2):veg_MDS2 +Species:veg_MDS1+Species:veg_MDS2+Pos_soil:veg_MDS1+Pos_soil:veg_MDS2,random=~1|PLOT/S_PLOT,data = Mengsong_field_envir,na.action=na.omit)
summary(veg_m2)
plot(veg_m2, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance has a decreasing trend. 

#Remove insignificant interaction terms Species:veg_MDS1 and Species:veg_MDS2, Pos_soil:veg_MDS1+Pos_soil:veg_MDS2
veg_m3<-update(veg_m2, ~.-Species:veg_MDS1-Species:veg_MDS2-Pos_soil:veg_MDS1-Pos_soil:veg_MDS2)

plot(veg_m3, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance has a decreasing trend. 
summary(veg_m3)

##Dealing with the variance trend
## model the variance as a power function of the fitted values

veg_m3a <- update(veg_m3, weights=varPower(form=~fitted(.)))
plot(veg_m3a, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## it is improved a lot

AIC(veg_m3a)

## model the variance as a exponential function of the fitted values
veg_m3b <- update(veg_m3, weights=varExp(form=~fitted(.)))
plot(veg_m3b, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## much similarto previous model m2a
AIC(veg_m3b)

## let's  make the variance be different 
veg_m3c <- update(veg_m3, weights=varIdent(form=~1|Species*Pos_soil))
plot(veg_m3c, abs(resid(., type='pearson'))~fitted(.)|Pos_soil, type=c('p', 'r'), 
     abline=0) 
plot(veg_m3c, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## Became worst than previous
AIC(veg_m3c)

AIC(veg_m3, veg_m3a,veg_m3b,veg_m3c)

summary(veg_m3b)
###Table
anova(veg_m3b, type = "marginal")
###Table

r.squaredLR(veg_m3b)

intervals(veg_m3b)
### Table

#results
#for species model
decomp_field_Dossa.results <- list()
decomp_field_Dossa.results$model.predictions<-preddat_orig
decomp_field_Dossa.results$r.squared<-r.squaredLR(m2b)
decomp_field_Dossa.results$confint_m2b<-intervals(m2b)
decomp_field_Dossa.results$anova_m2b<-anova(m2b, type="marginal")
decomp_field_Dossa.results$summary_m2b<-summary(m2b)


#for soil model
decomp_field_Dossa.results$soil$r.squared<- r.squaredLR(soil_m3) 
decomp_field_Dossa.results$soil$confint_soil_m3<-intervals(soil_m3)
decomp_field_Dossa.results$soil$anova_soil_m3<-anova(soil_m3, type="marginal")
decomp_field_Dossa.results$soil$summary_soil_m3<-summary(soil_m3)

#for veg model
decomp_field_Dossa.results$veg$r.squared<- r.squaredLR(veg_m3b) 
decomp_field_Dossa.results$veg$confint_veg_m2b<-intervals(veg_m3b)
decomp_field_Dossa.results$veg$anova_veg_m2b<-anova(veg_m3b, type="marginal")
decomp_field_Dossa.results$veg$summary_veg_m2b<-summary(veg_m3b)


save(decomp_field_Dossa.results, file="decomp_field_Dossa.results.Rdata")

load("decomp_field_Dossa.results.Rdata")

str(decomp_field_Dossa.results)
decomp_field_Dossa.results$model.predictions

decomp_field_Dossa.results$anova_m2b
decomp_field_Dossa.results$confint_m2b
decomp_field_Dossa.results$r.squared
decomp_field_Dossa.results$summary_m2b

decomp_field_Dossa.results$soil$anova
decomp_field_Dossa.results$soil$confint
decomp_field_Dossa.results$soil$r.squared
decomp_field_Dossa.results$soil$summary_soil_m3

decomp_field_Dossa.results$veg$anova_veg_m2b
decomp_field_Dossa.results$veg$confint
decomp_field_Dossa.results$veg$r.squared
decomp_field_Dossa.results$veg$summary_veg_m2b



###########MASS LOSS AND DECAY CONSTANT RATE k ANALYSIS

#####Compute decay constant rate k for individual logs using the formula from Oberle et al 2018 
Mengsong_field$k_ind<--365.25*(log(Mengsong_field$m_harvest/Mengsong_field$m_initial)/Mengsong_field$t)


str(Mengsong_field)
#####NEED TO ONLY PUT k_ind for only 36 mo because only at that time we measured the remaining mass#####

Mengsong_field$k_ind<-ifelse(Mengsong_field$Coll_No.x=="Coll_6",Mengsong_field$k_ind,NA)

summary(Mengsong_field$k_ind)
hist(Mengsong_field$k_ind)
which(Mengsong_field$k_ind>1)



head(Mengsong_field)
summary(Mengsong_field)



###Survival analysis
###This analysis asks the question 
####...whether the species differ in 
####...terms of their susceptibility to termite attack



#######Infested logs
###### 
Mengsong_field_down<-Mengsong_field[Mengsong_field$Pos_soil=="down",]
head(Mengsong_field_down)


Mengsong_field_down$Termi.assum

###Because for litsea at coll5 and 6 the plot #119 data is missing we need to exclude that 
###plot info from the data set. Also for plot #75 the information is missing for coll 6 thus
###we need to exclude it as well from other collections otherwise this will biase the proportion 
###calculation of termites presence abscence

####Use Surv() with an Ho being no difference in
survdiff(Surv(t,Termi.assum) ~ Species, data = Mengsong_field_down, rho=1)


temp_s_ass <- Surv(Mengsong_field_down$t,Mengsong_field_down$Termi.assum)
temp_fit_ass <- coxph(temp_s_ass~1)
temp_gr_ass <- survfit(temp_fit_ass)
plot(temp_gr_ass)
print(temp_gr_ass)

plot(survfit(data=Mengsong_field_down,Surv(t,Termi.assum)~Species), lty=c(2,4), xlab="Incubation duration (Days)")

surv.model2.logi_ass<-survreg(data=Mengsong_field_down,Surv(t,Termi.assum)~Species+For_type, dist="logistic")
summary(surv.model2.logi_ass)
str(Mengsong_field_down)


#####Table 5
##Put time t in month by dividing it by 30
surv.model2.logi_ass<-survreg(data=Mengsong_field_down,Surv((t/30),Termi.assum)~Species+For_type, dist="logistic")
summary(surv.model2.logi_ass)



##Check for iinteraction
surv.model.logi_ass3<-survreg(data=Mengsong_field_down,Surv(t/30,Termi.assum)~Species*For_type, dist="logistic")
summary(surv.model.logi_ass3)

##None of the interaction were significant so there were not included


#Interpretation
#From these results it shows none of the wood species was more susceptible to termite infestation... 
##... and open canopy is also more susceptible to termite infestation than mature forest. 





###Change the reference to CC to check difference between CC and OL
levels(Mengsong_field_down$For_type)
Mengsong_field_down$For_type<-relevel(Mengsong_field_down$For_type, ref="OC")
levels(Mengsong_field_down$For_type)

surv.model2.logi_ass2<-survreg(data=Mengsong_field_down,Surv((t/30),Termi.assum)~Species+For_type, dist="logistic")
summary(surv.model2.logi_ass2)

###Tables S

###Let's do plotting of number of infested logs vs non infested per species and across habitats
str(Mengsong_field_down)


###Let's build the contignecy table
####Get reference back to "CC"
levels(Mengsong_field_down$For_type)
Mengsong_field_down$For_type<-relevel(Mengsong_field_down$For_type, ref="CC")
levels(Mengsong_field_down$For_type)

###Load library(dplyr)
library(dplyr)


####Input the data hat gives count of log infested and non-infested by termite
termite<-read.csv(file="Termite_infestation_per_plot.csv", header=T, sep=",")
termite$term_prop<-(termite$n_presence/(termite$n_absence+termite$n_presence))*100
head(termite)
summary(termite)
termite$term_prop


###subset the data from 6 mo
termi_conti_pres_abs_6mo2<-termite[termite$Coll_No.x.x!="Coll_1",]
##Because we did not to a presence absence of termites systematically for all logs at time mo
###We need to set the recordwe have at that time to NA
#termi_conti_pres_abs_6mo2$term_prop<-ifelse(is.na(termi_conti_pres_abs_6mo2$term_prop)==T,0,termi_conti_pres_abs_6mo2$term_prop)
#termi_conti_pres_abs_6mo2$term_prop<-ifelse(termi_conti_pres_abs_6mo2$Coll_No.x.x=="Coll_2",0,termi_conti_pres_abs_6mo2$term_prop)

ggplot(data=termi_conti_pres_abs_6mo2, aes(x=Coll_No.x.x,y=term_prop,col=For_type.x, shape=Species.x.x))+
  geom_point(position=position_dodge(1),stat="identity")



####Reduce this table to mean and SD
termi_conti3<-termi_conti_pres_abs_6mo2%>%
  group_by(Species.x.x,For_type.x,Coll_No.x.x)%>%
  summarise(mean_term_prop=mean(term_prop,na.rm = T),
            sd_term_prop=sd(term_prop,na.rm = T))#,
termi_conti3

#tiff(filename="Figure with shapes and symboles for up down with 30 mo 2.tiff", res = 400, width=5700, height=3000, compression = "lzw")
ggplot(data=termi_conti3, aes(y=mean_term_prop, x=as.numeric(Coll_No.x.x), colour= For_type.x,shape=Species.x.x))+
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  geom_rect(data=NULL,aes(xmin=3.5,xmax=4.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  geom_rect(data=NULL,aes(xmin=5.5,xmax=6.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  geom_point(position=position_dodge(1),stat="identity", size=6, aes(y=mean_term_prop, x=Coll_No.x.x))+ 
  geom_errorbar(aes(ymax = termi_conti3$mean_term_prop+termi_conti3$sd_term_prop, ymin = termi_conti3$mean_term_prop-termi_conti3$sd_term_prop), width=0, size = 0.8, position=position_dodge(1))+
  scale_x_discrete(limits=levels(termi_conti3$Coll_No.x.x))+
  xlab(" Incubation duration") +
  ylab("Percentage of infested logs") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        legend.position = c(0.5,0.7),
        legend.text=element_text(size=12,face="italic"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 


###Figures for publications

####Make the interactive plot
dim(termi_conti3)
levels(termi_conti3$Species.x.x)

##Rename Species levels to "Castanopsis mekongensis", "Litsea cubeba"
levels(termi_conti3$Species.x.x)<-c("Castanopsis mekongensis", "Litsea cubeba")
levels(termi_conti3$Species.x.x)
##Rename Coll_No.x levels to "3 mo","6 mo","12 mo","18 mo","24 mo","36 mo"
levels(termi_conti3$Coll_No.x.x)
levels(termi_conti3$Coll_No.x.x)<-c("3 mo","6 mo","12 mo","18 mo","24 mo","36 mo")
###Reducing the levels of collection to only the one that are represented in the dataset from 6 mo
termi_conti3[] <- lapply(termi_conti3, function(x) if(is.factor(x)) factor(x) else x)
levels(termi_conti3$Coll_No.x.x)
##Renaming habitats category
levels(termi_conti3$For_type.x)
levels(termi_conti3$For_type.x)<-c("Mature forest", "Regenerating forest", "Open land")
levels(termi_conti3$For_type.x)

###Because we don't have systematic check for presence absence of termite at 6 mo collection 
###We need to leave it blank not to bias the results
termi_conti3$term_prop<-ifelse(termi_conti3$Coll_No.x.x=="6 mo",NA,termi_conti3$mean_term_prop)
termi_conti3$sd_term_prop<-ifelse(termi_conti3$Coll_No.x.x=="6 mo",NA,termi_conti3$sd_term_prop)

Percent_logs_infested<-ggplot(data=termi_conti3, aes(y=mean_term_prop, x=as.numeric(Coll_No.x.x), colour= For_type.x,shape=Species.x.x))+
  geom_rect(data=NULL,aes(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  geom_rect(data=NULL,aes(xmin=3.5,xmax=4.5,ymin=-Inf,ymax=Inf),
            fill="lightgray", colour = NA, alpha = 0.1) +
  geom_point(position=position_dodge(1),stat="identity", size=6, aes(y=term_prop, x=Coll_No.x.x))+ 
  geom_errorbar(aes(ymax = termi_conti3$mean_term_prop+termi_conti3$sd_term_prop, ymin = termi_conti3$mean_term_prop-termi_conti3$sd_term_prop), width=0, size = 0.8, position=position_dodge(1))+
  scale_x_discrete(limits=levels(termi_conti3$Coll_No.x.x))+
  xlab(" Incubation duration") +
  ylab("Percentage of infested logs") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        legend.text=element_text(size=12,face="italic"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

Percent_logs_infested

###Figure 5

tiff(filename="Figure 5 final error bar of percentage of log infested by termites per species per forest type.tiff", res = 600, width=6000, height=3400, compression = "lzw")
Percent_logs_infested
dev.off()


#####February 1st 2020

######

head(Mengsong_field[Mengsong_field$Coll_No.x=="Coll_6",])
str(Mengsong_field[Mengsong_field$Coll_No.x=="Coll_6",])
summary(Mengsong_field[Mengsong_field$Coll_No.x=="Coll_6",])


###Subset data collection at the end Coll_6
mass.loss.36mo<-Mengsong_field[Mengsong_field$Coll_No.x=="Coll_6",]



###Reduced by only one position down
mass.loss.36mo.down<-mass.loss.36mo[mass.loss.36mo$Pos_soil=="down",]
dim(mass.loss.36mo)
dim(mass.loss.36mo.down)

summary(mass.loss.36mo.down)
hist(mass.loss.36mo.down$ML_percent)
hist(log(mass.loss.36mo.down$ML_percent))
hist(sqrt(mass.loss.36mo.down$ML_percent))

###Better to use the non transformed data

levels(mass.loss.36mo.down$For_type)


###Modeling mass loss
####Use lme

ggplot(data=Mengsong_field, aes(x=Coll_No.x, y=Per_WSG_loss, group=For_type)) +
  geom_point() + geom_smooth(method='lm') +
  facet_wrap(~Species+Pos_soil)

ggplot(data=Mengsong_field[Mengsong_field$Coll_No.y=="Coll_6",], aes(x=Coll_No.x, y=ML_percent, color=For_type, shape=Species)) +
  geom_point(position=position_dodge(1)) + geom_smooth(method='lm')



###Graph  
str(mass.loss.36mo.down)
##Rename For_type levels to "Mature forest", "Regenerating Forest", "Open land"
levels(mass.loss.36mo.down$For_type)
levels(mass.loss.36mo.down$For_type)<-c("Mature forest", "Regenerating Forest", "Open land")
levels(mass.loss.36mo.down$For_type)

levels(mass.loss.36mo.down$Species)
##Rename Species levels 

levels(mass.loss.36mo.down$Species)<-c("Castanopsis mekongensis", "Litsea cubeba")
levels(mass.loss.36mo.down$Species)


####Mass loss modeling
#####Species and termites interaction, Forest type and termites
mass.loss.lme1<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species*as.factor(Termi.assum)+For_type+as.factor(Termi.assum)*For_type, random=~1|PLOT,na.action=na.omit)
summary(mass.loss.lme1)

plot(mass.loss.lme1, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance trend not bad.

plot(mass.loss.lme1)

####
AIC(mass.loss.lme1)

##Dealing with the variance trend

mass.loss.lme2 <- update(mass.loss.lme1, weights=varExp(form=~fitted(.)))
summary(mass.loss.lme2)
plot(mass.loss.lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## variance slight increase


## let's  make the variance be different for each species
mass.loss.lme3 <- update(mass.loss.lme1, weights=varIdent(form=~PLOT|Species))
summary(mass.loss.lme3)
plot(mass.loss.lme3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## the trend Worsens again

## let's  make the variance be different for each species
mass.loss.lme4 <- update(mass.loss.lme1, weights=varIdent(form=~1|Species))
summary(mass.loss.lme4)
plot(mass.loss.lme4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## the trend Worsens again

## let's  make the variance
mass.loss.lme5 <- update(mass.loss.lme1, weights=varPower(form=~fitted(.)))
summary(mass.loss.lme5)
plot(mass.loss.lme5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## the trend has a bit increasing trend 

####Species and termites Species and forest type interactions
mass.loss.lme6<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species*as.factor(Termi.assum)+For_type+Species*For_type, random=~1|PLOT, na.action=na.omit)
summary(mass.loss.lme6)
plot(mass.loss.lme6, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## variance very good
AIC(mass.loss.lme6)
r.squaredGLMM(mass.loss.lme6)

mass.loss.lme7<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species*as.factor(Termi.assum)+For_type, random=~1|PLOT,na.action=na.omit)
summary(mass.loss.lme7)
plot(mass.loss.lme7, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## variance very good


AIC(mass.loss.lme7)

####

summary(mass.loss.lme1)
summary(mass.loss.lme2)
summary(mass.loss.lme3)
summary(mass.loss.lme4)
summary(mass.loss.lme5)
summary(mass.loss.lme6)
AIC(mass.loss.lme1,mass.loss.lme2, mass.loss.lme3,mass.loss.lme4,mass.loss.lme5, mass.loss.lme6)

anova(mass.loss.lme1,mass.loss.lme2)
anova(mass.loss.lme2,mass.loss.lme3)
anova(mass.loss.lme3,mass.loss.lme4)
anova(mass.loss.lme4,mass.loss.lme5)


##Best model mass.loss.lme 6
summary(mass.loss.lme6)
r.squaredGLMM(mass.loss.lme6) ###r^2 m== 25.12 %, r^2c ==44.95%



###Let's change the reference level
levels(mass.loss.36mo.down$For_type)
mass.loss.36mo.down$For_type<-relevel(mass.loss.36mo.down$For_type, ref="Regenerating Forest")

mass.loss.lme6_2<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species*as.factor(Termi.assum)+For_type+Species*For_type, random=~1|PLOT, na.action=na.omit)
summary(mass.loss.lme6_2)
plot(mass.loss.lme6_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## variance very good
AIC(mass.loss.lme6_2)
r.squaredGLMM(mass.loss.lme6_2)


###Let's change back the reference level
levels(mass.loss.36mo.down$For_type)
mass.loss.36mo.down$For_type<-relevel(mass.loss.36mo.down$For_type, ref="Mature forest")



###Let's proceed in simplification of the model
###Remove Termite habitat interaction

mass.loss.lme7<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species*as.factor(Termi.assum)+For_type, random=~1|PLOT,na.action=na.omit)
summary(mass.loss.lme7)
plot(mass.loss.lme7, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## variance very good

#AIC(mass.loss.lme6,mass.loss.lme7)

##anova with type="marginal" gives the same result with summary
anova(mass.loss.lme7, type="marginal")
summary(mass.loss.lme7)

r.squaredGLMM(mass.loss.lme7)

#### mass model Without termites

mass.loss.lme8<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species+For_type, random=~1|PLOT,na.action=na.omit)
summary(mass.loss.lme8)
plot(mass.loss.lme8, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0) ## variance slight decrease trend

plot(mass.loss.lme8, resid(., type='p')~fitted(.)|Species*For_type*as.factor(Termi.assum), 
     layout=c(3,4)) ## 

#AIC(mass.loss.lme6,mass.loss.lme7)

##anova with type="marginal" gives the same result with summary
anova(mass.loss.lme8, type="marginal")
summary(mass.loss.lme8)

r.squaredGLMM(mass.loss.lme8)###r^2 m== 23.05 %, r^2c ==49.11%


###Mass loss model with interaction Species with habitat releveled
levels(mass.loss.36mo.down$For_type)
mass.loss.36mo.down$For_type<-relevel(mass.loss.36mo.down$For_type, ref="Regenerating Forest")

mass.loss.lme3<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species*as.factor(Termi.assum)+Species*For_type,random=~1|PLOT,na.action=na.omit)
summary(mass.loss.lme3)#
anova(mass.loss.lme3)
plot(mass.loss.lme3, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance has a slight decreasing trend
par(mfrow=c(2,2))
plot(mass.loss.lme3)
r.squaredGLMM(mass.loss.lme3)

##Back the default baseline for habitat
levels(mass.loss.36mo.down$For_type)
mass.loss.36mo.down$For_type<-relevel(mass.loss.36mo.down$For_type, ref="Mature forest")

####
###Mass loss model with interaction Species with termite, termite with habitat
mass.loss.lme4<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species*as.factor(Termi.assum)+as.factor(Termi.assum)*For_type,random=~1|PLOT,na.action=na.omit)
summary(mass.loss.lme4)
anova(mass.loss.lme4)
plot(mass.loss.lme4, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance good


###Mass loss model with three way interaction Species with termite and with habitat
mass.loss.lme5<-lme(data=mass.loss.36mo.down, fixed=ML_percent~Species*as.factor(Termi.assum)*For_type,random=~1|PLOT,na.action=na.omit)
summary(mass.loss.lme5)
anova(mass.loss.lme5)
plot(mass.loss.lme5, abs(resid(.))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance has
r.squaredGLMM(mass.loss.lme5)

#####
ggplot(data=mass.loss.36mo.down, aes(x=For_type,  y=ML_percent, shape=Species, color=as.factor(Termi.assum)))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  xlab(" Land cover type") +
  ylab("Percentage mass loss")+ 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        legend.position = c(0.5,0.7),
        legend.text=element_text(size=12,face="italic"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 


OL_logs<-mass.loss.36mo.down[mass.loss.36mo.down$For_type=="OL",]
OL_logs_litsea<-OL_logs[OL_logs$Species=="Litsea cubeba",]

###Piping
mass.loss.36mo.down$Termite_A_P<-as.factor(mass.loss.36mo.down$Termi.assum)
sum_mass_loss<-mass.loss.36mo.down%>%
  group_by(Species, For_type, Termite_A_P)%>%
    summarize(mean_ML_percent=mean(ML_percent),
              se_ML_percent=sd(ML_percent)/sqrt(n())##n() counts rows in each group
              )

sum_mass_loss
str(sum_mass_loss)

##renaming Termite_A_P levels from 0 and 1 to Absence and Presence
levels(sum_mass_loss$Termite_A_P)
levels(sum_mass_loss$Termite_A_P)<-c("Absence","Presence")
levels(sum_mass_loss$Termite_A_P)

##renaming For_type levels from "CC" "OC" and "OL" to "Mature forest" "Regenerating forest" "Open land"
levels(sum_mass_loss$For_type)
levels(sum_mass_loss$For_type)<-c("Mature forest", "Regenerating forest", "Open land")
levels(sum_mass_loss$For_type)


mass_loss_with_regard_termite<-ggplot(data=sum_mass_loss,aes(x=For_type,  y=mean_ML_percent, shape=Species, color=Termite_A_P))+
geom_point(position=position_dodge(0.5),stat="identity", size=6)+ 
geom_errorbar(data=sum_mass_loss, aes(ymax =mean_ML_percent+se_ML_percent, ymin = mean_ML_percent-se_ML_percent), width=0, size = 0.8, position=position_dodge(0.5))+
  xlab("") + # Land cover type
  ylab("Percentage mass loss")+ 
  labs(colour="Termite status")+
  labs(shape="Species")+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_line(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        legend.text=element_text(size=12,face="italic"),
        legend.title = element_text(size=15),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 


###Figure S5 Mass loss along distrubance gradient categories
tiff(filename="Figure S5 mass loss with regard to termite status.tiff", res = 600, width=6000, height=3400, compression = "lzw")
mass_loss_with_regard_termite
dev.off()



#####DECAY CONSTANT RATE Analysis
######Individual log decay analysis
str(mass.loss.36mo.down)
summary(mass.loss.36mo.down)
hist(mass.loss.36mo.down$k_ind)
###Need transformation because this is at the moment skewed on the right

hist(log(mass.loss.36mo.down$k_ind))
#weights=varIdent(form=~species|bark_treatment.x)



##Back the default baseline for habitat
levels(mass.loss.36mo.down$For_type)
mass.loss.36mo.down$For_type<-relevel(mass.loss.36mo.down$For_type, ref="Mature forest")

ggplot(data=mass.loss.36mo.down, aes(x=For_type,  y=k_ind, shape=Species, color=as.factor(Termi.assum)))+
  geom_point(position=position_dodge(1),stat="identity", size=6)+ 
  xlab("") +
  ylab(bquote("Decay rate constant (k) based on mass loss ( "* year^-1*")"))+ 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        legend.position = c(0.5,0.7),
        legend.text=element_text(size=12,face="italic"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 


OL_logs<-mass.loss.36mo.down[mass.loss.36mo.down$For_type=="OL",]
OL_logs_litsea<-OL_logs[OL_logs$Species=="Litsea cubeba",]

###Piping
sum_decay<-mass.loss.36mo.down%>%
  group_by(Species, For_type, Termite_A_P)%>%
  summarize(mean_k_ind=mean(k_ind),
            se_k_ind=sd(k_ind)/sqrt(n())##n() counts rows in each group
  )

sum_decay
str(sum_decay)

##renaming Termite_A_P levels from 0 and 1 to Absence and Presence
levels(sum_decay$Termite_A_P)
levels(sum_decay$Termite_A_P)<-c("Absence","Presence")
levels(sum_decay$Termite_A_P)

##renaming For_type levels from "CC" "OC" and "OL" to "Mature forest" "Regenerating forest" "Open land"
levels(sum_decay$For_type)
levels(sum_decay$For_type)<-c("Mature forest", "Regenerating forest", "Open land")
levels(sum_decay$For_type)


sum_decay_with_regard_termite<-ggplot(data=sum_decay,aes(x=For_type,  y=mean_k_ind, shape=Species, color=Termite_A_P))+
  geom_point(position=position_dodge(0.5),stat="identity", size=6)+ 
  geom_errorbar(data=sum_decay, aes(ymax =mean_k_ind+se_k_ind, ymin = mean_k_ind-se_k_ind), width=0, size = 0.8, position=position_dodge(0.5))+
  xlab("") + # Land cover type
  ylab(bquote("Decay rate constant (k) based on mass loss ( "* year^-1*")"))+  
  labs(colour="Termite status")+
  labs(shape="Species")+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_line(),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        legend.position=c(0,1), 
        legend.justification=c(0,1),
        #legend.position = c("bottom"),
        #legend.position = c(0.5,0.7),
        legend.text=element_text(size=12,face="italic"),
        legend.title = element_text(size=15),
        legend.background = element_blank(),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 
sum_decay_with_regard_termite

###Figure 3 Decay constant rate 
tiff(filename="Figure 3 decay with regard to termite status 2.tiff", res = 600, width=6000, height=3400, compression = "lzw")
sum_decay_with_regard_termite
dev.off()


#####There might be spatial autocorrelation, our samples here have true replicate at PLOT level, so a random factor here is PLOT
###Full model three way interaction
Decay.lme0<-lme(data=mass.loss.36mo.down, fixed=log(k_ind)~Species*For_type*as.factor(Termi.assum), random=~1|PLOT)

summary(Decay.lme0)##
anova(Decay.lme0)
plot(Decay.lme0, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Decreasing trend


Decay.lme0_1<-update(Decay.lme0,weights=varExp(form=~fitted(.)))
summary(Decay.lme0_1)##
anova(Decay.lme0_1)
plot(Decay.lme0_1, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly increasing trend


Decay.lme0_2<-update(Decay.lme0,weights=varPower(form=~fitted(.)))
summary(Decay.lme0_2)##
anova(Decay.lme0_2)
plot(Decay.lme0_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly better

Decay.lme0_3<-update(Decay.lme0,weights=varIdent(form=~1|Species*as.factor(Termi.assum)))
summary(Decay.lme0_3)##
anova(Decay.lme0_3)
plot(Decay.lme0_3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## worsened trend

## let's  make the variance be different for each species
#weights=varIdent(form=~PLOT|Species)
#weights=varExp(form=~fitted(.))

AIC(Decay.lme0,Decay.lme0_1,Decay.lme0_2,Decay.lme0_3)

####Let's use Decay.lme0_1

####let's do simplification
Decay.lme0_1<-lme(data=mass.loss.36mo.down, fixed=log(k_ind)~Species*as.factor(Termi.assum)*For_type, random=~1|PLOT)
summary(Decay.lme0_1)##
anova(Decay.lme0_1)
plot(Decay.lme0_1, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly decreasing trend
r.squaredGLMM(Decay.lme0_1)
###Remove three way interaction

Decay.lme0_4<-update(Decay.lme0_1, .~.-Species:as.factor(Termi.assum):For_type)
summary(Decay.lme0_4)##
anova(Decay.lme0_4)
plot(Decay.lme0_4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly decreasing trend

anova(Decay.lme0_1,Decay.lme0_4)
AIC(Decay.lme0_1,Decay.lme0_4)
r.squaredGLMM(Decay.lme0_4)

####accounting for variance trend by allowing different variance per species per forest
Decay.lme0_4_2<-update(Decay.lme0_1, .~.-Species:as.factor(Termi.assum):For_type,weights=varIdent(form=~1|For_type*Species))
summary(Decay.lme0_4_2)##
anova(Decay.lme0_4_2)
plot(Decay.lme0_4_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly decreasing trend

anova(Decay.lme0_1,Decay.lme0_4_2)
AIC(Decay.lme0_1,Decay.lme0_4_2)
r.squaredGLMM(Decay.lme0_4_2)


###Remove two way interaction forest:termite

Decay.lme0_5<-update(Decay.lme0_4, .~.-as.factor(Termi.assum):For_type,weights=varIdent(form=~1|For_type*Species))
summary(Decay.lme0_5)##
anova(Decay.lme0_5)
plot(Decay.lme0_5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly decreasing trend

anova(Decay.lme0_4,Decay.lme0_5)
AIC(Decay.lme0_4,Decay.lme0_5)
r.squaredGLMM(Decay.lme0_5)


###Remove two way interaction species:forest

Decay.lme0_6<-update(Decay.lme0_5, .~.-Species:For_type,weights=varIdent(form=~1|For_type))
summary(Decay.lme0_6)##
anova(Decay.lme0_6)
plot(Decay.lme0_6, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly decreasing trend

anova(Decay.lme0_5,Decay.lme0_6)
AIC(Decay.lme0_5,Decay.lme0_6)
r.squaredGLMM(Decay.lme0_6)


###Allow variance to change across forest types and species
Decay.lme0_7<-update(Decay.lme0_6, .~.,weights=varIdent(form=~1|For_type*Species))
summary(Decay.lme0_7)##
anova(Decay.lme0_7)
plot(Decay.lme0_7, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly decreasing trend

anova(Decay.lme0_6,Decay.lme0_7)
AIC(Decay.lme0_6,Decay.lme0_7)
r.squaredGLMM(Decay.lme0_7)

###Remove two way interaction species:termites
Decay.lme0_8<-update(Decay.lme0_6, .~.-Species:as.factor(Termi.assum),weights=varIdent(form=~1|For_type*Species))
summary(Decay.lme0_8)##
anova(Decay.lme0_8)
plot(Decay.lme0_8, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##  decreasing trend

anova(Decay.lme0_7,Decay.lme0_8)
AIC(Decay.lme0_7,Decay.lme0_8)
r.squaredGLMM(Decay.lme0_8)


####BEST MODEL Decay.lme0_6

summary(Decay.lme0_6)##
anova(Decay.lme0_6, type="marginal")
r.squaredGLMM(Decay.lme0_6)#### r^2marg=25.50% r^2cond==46.66%

####Change reference 
levels(mass.loss.36mo.down$For_type)
mass.loss.36mo.down$For_type<-relevel(mass.loss.36mo.down$For_type, ref="Regenerating Forest")
levels(mass.loss.36mo.down$For_type)

Decay.lme0_6_2<-update(Decay.lme0_5, .~.-Species:For_type,weights=varIdent(form=~1|For_type))
summary(Decay.lme0_6_2)##
anova(Decay.lme0_6_2)
plot(Decay.lme0_6_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly decreasing trend

r.squaredGLMM(Decay.lme0_6_2)

####Change back to default reference 
levels(mass.loss.36mo.down$For_type)
mass.loss.36mo.down$For_type<-relevel(mass.loss.36mo.down$For_type, ref="Mature forest")
levels(mass.loss.36mo.down$For_type)



#####Decay model without termite

Decay.lme5<-lme(data=mass.loss.36mo.down, fixed=log(k_ind)~Species+For_type, random=~1|PLOT)
summary(Decay.lme5)
anova(Decay.lme5)
par(mfrow=c(2,2))
plot(Decay.lme5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## decreasing trend
plot(Decay.lme5)

### Fix the variance decresing trend by allowing different variance for forest type


Decay.lme6<-lme(data=mass.loss.36mo.down, fixed=log(k_ind)~Species+For_type, random=~1|PLOT,weights=varIdent(form=~1|For_type))
summary(Decay.lme6)
anova(Decay.lme6)
par(mfrow=c(2,2))
plot(Decay.lme6, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## decreasing trend

r.squaredGLMM(Decay.lme6)####Without termite the decay model explains r^2marg=23.25%, r^2cond=49.43%




####Decay  and environmental and soil parameters

########USE LME

####Decay  and environmental and soil parameters

### Merge environmental, plant diversity factors with Mengsong field data 
####Import envrionmental factors

Mengsong_envi<- read.csv(file="Mengsong_environmental_factors.csv", header=T, sep=",")

str(Mengsong_envi)
summary(is.na(Mengsong_envi$FOR_TYPE))
Mengsong_envi$PLOT_SBPLOT

###Create plot_subplot for merging later
mass.loss.36mo.down
str(mass.loss.36mo.down)
mass.loss.36mo.down$plot_sbplot<-as.factor(paste(mass.loss.36mo.down$PLOT,mass.loss.36mo.down$S_PLOT, sep="_"))

mass.loss.36mo.down$plot_sbplot
levels(mass.loss.36mo.down$plot_sbplot)

####Merge now
mengsong_k_envi<-merge(mass.loss.36mo.down, Mengsong_envi, by.x ="plot_sbplot"  , by.y = "PLOT_SBPLOT", all.x=T)
head(mengsong_k_envi)
dim(mengsong_k_envi)
str(mengsong_k_envi)

mengsong_k_envi$Termi.assum
mengsong_k_envi$Termite_A_P

####Model k constant with 
####Consider adding vegetation
lme_k_veg<-lme(data=mengsong_k_envi,log(k_ind)~Species+veg_MDS1+veg_MDS2, random=~1|PLOT.x)
summary(lme_k_veg)
anova(lme_k_veg, type="marginal")
AIC(lme_k_veg)
plot(lme_k_veg, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## decreasing trend

###Simplify by removing Veg_MDS2

lme_k_veg2<-lme(data=mengsong_k_envi,log(k_ind)~Species+veg_MDS1, random=~1|PLOT.x)
summary(lme_k_veg2)##With forest type replaced by vegetation composition r^2 marg=24.27, r^2cond=45.95%
anova(lme_k_veg2, type="marginal")
AIC(lme_k_veg2)
plot(lme_k_veg2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## decreasing trend

r.squaredGLMM(lme_k_veg2)



###Improving model by allowing different variance per vegetaion type
lme_k_veg2_2<-lme(data=mengsong_k_envi,log(k_ind)~Species+veg_MDS1, random=~1|PLOT.x,weights=varIdent(form=~1|veg_MDS1*veg_MDS2))
summary(lme_k_veg2_2)##With forest type replaced by vegetation composition r^2 marg=36.69, r^2cond=77.32%
anova(lme_k_veg2_2, type="marginal")
AIC(lme_k_veg2_2)
plot(lme_k_veg2_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## increasing trend

r.squaredGLMM(lme_k_veg2_2)

#

###Improving model by allowing different variance per forest type
lme_k_veg2_2<-lme(data=mengsong_k_envi,log(k_ind)~Species+veg_MDS1, random=~1|PLOT.x,weights=varExp(form=~fitted(.)))
summary(lme_k_veg2_2)##With forest type replaced by vegetation composition r^2 marg=39.40, r^2cond=77.21%
anova(lme_k_veg2_2, type="marginal")
AIC(lme_k_veg2_2)
plot(lme_k_veg2_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## decreasing trend

r.squaredGLMM(lme_k_veg2_2)


####Consider adding vegetation and including termites and interactions with bth species and vegetation
lme_k_veg3<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+veg_MDS1*as.factor(Termi.assum)+veg_MDS2*as.factor(Termi.assum), random=~1|PLOT.x)
summary(lme_k_veg3) 
anova(lme_k_veg3, type="marginal")
AIC(lme_k_veg3)
plot(lme_k_veg3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##decreasing trend

####Simplify by removing the vegetation and termite interaction
lme_k_veg4<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+veg_MDS1+veg_MDS2, random=~1|PLOT.x)
summary(lme_k_veg4) ##
anova(lme_k_veg4, type="marginal")
AIC(lme_k_veg4)
plot(lme_k_veg4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##decreasing trend

### Simplify by removing veg_MDS2
lme_k_veg5<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+veg_MDS1, random=~1|PLOT.x)
summary(lme_k_veg5) ##
anova(lme_k_veg5, type="marginal")
AIC(lme_k_veg5)
plot(lme_k_veg5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)### Decreasing trend



###Simplify by removing species and termites interaction 

lme_k_veg6<-lme(data=mengsong_k_envi,log(k_ind)~Species+as.factor(Termi.assum)+veg_MDS1, random=~1|PLOT.x)
summary(lme_k_veg6)##With forest type replaced by vegetation composition
anova(lme_k_veg6, type="marginal")
AIC(lme_k_veg6)
plot(lme_k_veg6, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)### Decreasing trend

###Simplify by removing termites

lme_k_veg7<-lme(data=mengsong_k_envi,log(k_ind)~Species+veg_MDS1, random=~1|PLOT.x,weights=varIdent(form=~1|For_type))
summary(lme_k_veg7)##With forest type replaced by vegetation composition
anova(lme_k_veg7, type="marginal")
AIC(lme_k_veg7)
plot(lme_k_veg7, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Decreasing trend
r.squaredGLMM(lme_k_veg7)

####We chose lme_k_veg2_2
##Table S11
summary(lme_k_veg2_2)##With forest type replaced by vegetation composition r^2 marg=39.40, r^2cond=77.21%
anova(lme_k_veg2_2, type="marginal")
AIC(lme_k_veg2_2)
r.squaredGLMM(lme_k_veg2_2)



###mass loss and vegetation
dim(mengsong_k_envi)
hist(mengsong_k_envi$ML_percent)

lme_mass_veg<-lme(data=mengsong_k_envi,ML_percent~Species+veg_MDS1+veg_MDS2, random=~1|PLOT.x)
summary(lme_mass_veg)
anova(lme_mass_veg, type="marginal")
AIC(lme_mass_veg)
plot(lme_mass_veg, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## decreasing trend

###Simplify by removing Veg_MDS2

lme_mass_veg2<-lme(data=mengsong_k_envi,ML_percent~Species+veg_MDS1, random=~1|PLOT.x)
summary(lme_mass_veg2)##With forest type replaced by vegetation composition r^2 marg=24.27, r^2cond=45.95%
anova(lme_mass_veg2, type="marginal")
AIC(lme_mass_veg2)
plot(lme_mass_veg2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly decreasing trend

r.squaredGLMM(lme_mass_veg2)



###Improving model by allowing different variance per vegetaion type
lme_mass_veg2_2<-lme(data=mengsong_k_envi,ML_percent~Species+veg_MDS1, random=~1|PLOT.x,weights=varIdent(form=~1|veg_MDS1*veg_MDS2))
summary(lme_mass_veg2_2)
anova(lme_mass_veg2_2, type="marginal")
AIC(lme_mass_veg2_2)
plot(lme_mass_veg2_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## increasing trend

r.squaredGLMM(lme_mass_veg2_2)


###Improving model by allowing different variance per forest type
lme_mass_veg2_3<-lme(data=mengsong_k_envi,ML_percent~Species+veg_MDS1, random=~1|PLOT.x,weights=varExp(form=~fitted(.)))

###Table S12
summary(lme_mass_veg2_3)
anova(lme_mass_veg2_3, type="marginal")
AIC(lme_mass_veg2_3)
plot(lme_mass_veg2_3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## good

r.squaredGLMM(lme_mass_veg2_3)##r^2marg=16.52 % r^cond=31.13 %

###Without termites BEST MODEL lme_mass_veg2_3
summary(lme_mass_veg2_3)
anova(lme_mass_veg2_3, type="marginal")
r.squaredGLMM(lme_mass_veg2_3)##r^2marg=16.52 % r^cond=31.13 %

####Consider adding vegetation and including termites and interactions with both species and vegetation
lme_mass_veg3<-lme(data=mengsong_k_envi,ML_percent~Species*as.factor(Termi.assum)+veg_MDS1*as.factor(Termi.assum)+veg_MDS2*as.factor(Termi.assum), random=~1|PLOT.x)
summary(lme_mass_veg3) 
anova(lme_mass_veg3, type="marginal")
AIC(lme_mass_veg3)
plot(lme_mass_veg3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##slightly ecreasing

####Simplify by removing the vegetation and termite interaction
lme_mass_veg4<-lme(data=mengsong_k_envi,ML_percent~Species*as.factor(Termi.assum)+veg_MDS1+veg_MDS2, random=~1|PLOT.x)
summary(lme_mass_veg4) ##
anova(lme_mass_veg4, type="marginal")
AIC(lme_mass_veg4)
plot(lme_mass_veg4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##slightly decreasing trend

### Simplify by removing veg_MDS2
lme_mass_veg5<-lme(data=mengsong_k_envi,ML_percent~Species*as.factor(Termi.assum)+veg_MDS1, random=~1|PLOT.x)
summary(lme_mass_veg5) ##
anova(lme_mass_veg5, type="marginal")
AIC(lme_mass_veg5)
plot(lme_mass_veg5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)### slightly Decreasing trend

###Simplify by removing species and termites interaction 

lme_mass_veg6<-lme(data=mengsong_k_envi,ML_percent~Species+as.factor(Termi.assum)+veg_MDS1, random=~1|PLOT.x)
summary(lme_mass_veg6)##With forest type replaced by vegetation composition
anova(lme_mass_veg6, type="marginal")
AIC(lme_mass_veg6)
plot(lme_mass_veg6, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)

###Simplify by removing termites

lme_mass_veg7<-lme(data=mengsong_k_envi,ML_percent~Species+veg_MDS1, random=~1|PLOT.x)
summary(lme_mass_veg7)##With forest type replaced by vegetation composition
anova(lme_mass_veg7, type="marginal")
AIC(lme_mass_veg7)

plot(lme_mass_veg7, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)
r.squaredGLMM(lme_mass_veg7)

####We chose lme_mass_veg2_3

summary(lme_mass_veg2_3)##With forest type replaced by vegetation composition r^2 marg=16.52, r^2cond=31.13%
anova(lme_mass_veg2_3, type="marginal")
AIC(lme_mass_veg2_3)
r.squaredGLMM(lme_mass_veg2_3)





####Consider adding soil/topography
lme_k_soil<-lme(data=mengsong_k_envi,log(k_ind)~Species+SOIL_PC1+SOIL_PC2, random=~1|PLOT.x, na.action = na.omit,weights=varIdent(form=~1|For_type))
summary(lme_k_soil)
anova(lme_k_soil, type="marginal")
AIC(lme_k_soil)
plot(lme_k_soil, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###Decreasing singtrend
r.squaredGLMM(lme_k_soil)


###Simplify by removing SOIL_PC2

lme_k_soil2<-lme(data=mengsong_k_envi,log(k_ind)~Species+SOIL_PC1,random=~1|PLOT.x, na.action = na.omit,weights=varIdent(form=~1|For_type))
summary(lme_k_soil2)##With forest type replaced by vegetation composition
anova(lme_k_soil2, type="marginal")
AIC(lme_k_soil2)
plot(lme_k_soil2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Dcreasing trend
r.squaredGLMM(lme_k_soil)

###Simplify by removing SOIL_PC1

lme_k_soil3<-lme(data=mengsong_k_envi,log(k_ind)~Species,random=~1|PLOT.x, na.action = na.omit,weights=varIdent(form=~1|For_type))
summary(lme_k_soil3) 
anova(lme_k_soil3, type="marginal")
AIC(lme_k_soil3)
plot(lme_k_soil3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)# Decreasing trend
r.squaredGLMM(lme_k_soil3)


###Without termites abs/pres, Soil / topography does not explain k values variance

####Consider adding soil/topograpy and including termites and interactions with both species and soil/topography
lme_k_soil4<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+SOIL_PC1*as.factor(Termi.assum)+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x, na.action=na.omit)
summary(lme_k_soil4) ##
anova(lme_k_soil4, type="marginal")
AIC(lme_k_soil4)
plot(lme_k_soil4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###Decreasing trend
r.squaredGLMM(lme_k_soil4)

####Simplify by removing the SOIL_PC1 and termite interaction
lme_k_soil5<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+SOIL_PC1+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x,na.action=na.omit)
summary(lme_k_soil5) ##
anova(lme_k_soil5, type="marginal")
AIC(lme_k_soil5)
plot(lme_k_soil5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Decreasing trend
r.squaredGLMM(lme_k_soil5)


####Simplify by removing the SOIL_PC1
lme_k_soil6<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x,na.action=na.omit)
summary(lme_k_soil6) ##
anova(lme_k_soil6, type="marginal")
AIC(lme_k_soil6)
plot(lme_k_soil6, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###Decreasing trend
r.squaredGLMM(lme_k_soil6)


#### Deal with decreasing trend in variance by allowing k to have different variance per forest type
lme_k_soil7<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x,na.action=na.omit, weights=varIdent(form=~1|For_type))
summary(lme_k_soil7) ##
anova(lme_k_soil7, type="marginal")
AIC(lme_k_soil7)
plot(lme_k_soil7, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Slightly decreasing
anova(lme_k_soil7,lme_k_soil6)
r.squaredGLMM(lme_k_soil7)

###Remove interaction between 
#### Deal with decreasing trend in variance by allowing k to have different variance per forest type
lme_k_soil7_2<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+SOIL_PC2, random=~1|PLOT.x,na.action=na.omit, weights=varIdent(form=~1|For_type))
summary(lme_k_soil7_2) ##
anova(lme_k_soil7_2, type="marginal")
AIC(lme_k_soil7_2)
plot(lme_k_soil7_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)#Slightly decreasing trend
anova(lme_k_soil7_2,lme_k_soil7)
r.squaredGLMM(lme_k_soil7_2)


#### Deal with decreasing trend in variance by allowing k to have different variance per forest type
lme_k_soil8<-lme(data=mengsong_k_envi,log(k_ind)~Species*as.factor(Termi.assum)+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x,na.action=na.omit, weights=varIdent(form=~1|For_type))
summary(lme_k_soil8) ##
anova(lme_k_soil8, type="marginal")
AIC(lme_k_soil8)
plot(lme_k_soil8, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###decreasing trend
anova(lme_k_soil8,lme_k_soil7)
r.squaredGLMM(lme_k_soil8)


###BEST MODEL None that uses SOIL PC



#####Mass loss and envrionmental
####Soil/topography

####Consider adding soil/topography
lme_mass_soil<-lme(data=mengsong_k_envi,ML_percent~Species+SOIL_PC1+SOIL_PC2, random=~1|PLOT.x,na.action=na.omit)
summary(lme_mass_soil)
anova(lme_mass_soil, type="marginal")
AIC(lme_mass_soil)
plot(lme_mass_soil, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##slightly decreasing

r.squaredGLMM(lme_mass_soil)

###Simplify by removing SOIL_PC2

lme_mass_soil2<-lme(data=mengsong_k_envi,ML_percent~Species+SOIL_PC1, random=~1|PLOT.x,na.action=na.omit)
summary(lme_mass_soil2)##With forest type replaced by vegetation composition
anova(lme_mass_soil2, type="marginal")
AIC(lme_mass_soil2)
plot(lme_mass_soil2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###slightly decreasing

r.squaredGLMM(lme_mass_soil2)


###Simplify by removing SOIL_PC1

lme_mass_soil3<-lme(data=mengsong_k_envi,ML_percent~Species, random=~1|PLOT.x,na.action=na.omit)
summary(lme_mass_soil3)##With forest type replaced by vegetation composition
anova(lme_mass_soil3, type="marginal")
AIC(lme_mass_soil3)
plot(lme_mass_soil3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##slightly decreasing

r.squaredGLMM(lme_mass_soil3)###


####Consider adding soil/topograpy and including termites and interactions with both species and soil/topography
lme_mass_soil4<-lme(data=mengsong_k_envi,ML_percent~Species*as.factor(Termi.assum)+SOIL_PC1*as.factor(Termi.assum)+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x,na.action=na.omit)
summary(lme_mass_soil4) ##

anova(lme_mass_soil4, type="marginal")
AIC(lme_mass_soil4)
plot(lme_mass_soil4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##slightly increasing



####Simplify by removing the SOIL_PC1 and termite interaction
lme_mass_soil5<-lme(data=mengsong_k_envi,ML_percent~Species*as.factor(Termi.assum)+SOIL_PC1+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x,na.action = na.omit)
summary(lme_mass_soil5) ##
anova(lme_mass_soil5, type="marginal")
AIC(lme_mass_soil5)
plot(lme_mass_soil5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##not bad

####Simplify by removing the SOIL_PC1
lme_mass_soil6<-lme(data=mengsong_k_envi,ML_percent~Species*as.factor(Termi.assum)+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x, na.action=na.omit)
summary(lme_mass_soil6) 

anova(lme_mass_soil6, type="marginal")
AIC(lme_mass_soil6)
plot(lme_mass_soil6, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)####slightly increasing

####Simplify by removing species termite interaction
lme_mass_soil7<-lme(data=mengsong_k_envi,ML_percent~Species+as.factor(Termi.assum)+SOIL_PC2*as.factor(Termi.assum), random=~1|PLOT.x, na.action=na.omit)
summary(lme_mass_soil7) 

anova(lme_mass_soil7, type="marginal")
AIC(lme_mass_soil7)
plot(lme_mass_soil7, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)
r.squaredGLMM(lme_mass_soil7)###r^2 marg=20.25%, r^2cond=45.00%



####Simplify by removing species termite interaction
lme_mass_soil8<-lme(data=mengsong_k_envi,ML_percent~Species, random=~1|PLOT.x, na.action=na.omit)
summary(lme_mass_soil8) 

anova(lme_mass_soil8, type="marginal")
AIC(lme_mass_soil8)
plot(lme_mass_soil8, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##slightly decreasing
r.squaredGLMM(lme_mass_soil8)###r^2 marg=19.17%, r^2cond=49.08%


###BEST MODEL none 

################################################


##Now compute the mean and SE for k for each forest 

##Compute the means
##Piping method
library(dplyr)
library(tidyr)
#summary(log_k_envi)
##Using piping method in R piping is done with %>%
k_ind_Means<-mengsong_k_envi%>%
  group_by(For_type,Species)%>%
  summarise(mean_k_ind=mean(k_ind),
            sd_k_ind=sd(k_ind))
k_ind_Means
k_ind_Means<-as.data.frame(k_ind_Means)
k_ind_Means

#####Litsea to Castonopsis ratio
####Mature forest  0.3962090/0.2410064=1.64
####Regenerating forest 0.353066/ 0.2512412=1.4
####Open land  0.5466251/0.3513492=1.55


#####Forest type ratio
######Open to mature for Litsea 0.5466251/0.3962090=1.38
######Open to mature for Casstanopsis 0.3513492/0.2410064=1.46


####Considering termite status
##Using piping method in R piping is done with %>%
k_ind_termites_Means<-mengsong_k_envi%>%
  group_by(For_type,Species, Termite_A_P)%>%
  summarise(mean_k_ind=mean(k_ind),
            sd_k_ind=sd(k_ind))
k_ind_termites_Means
k_ind_termites_Means<-as.data.frame(k_ind_termites_Means)
k_ind_termites_Means

#####Litsea to Castonopsis ratio
#####WITHOUT TERMITES 
####Mature forest  0.3730903/0.1920045=1.94
####Regenerating forest 0.3424894/0.2069781=1.65
####Open land  0.5834800/0.2768242=2.11

###Mean= 1.94, mean(1.94,1.65,2.11)

#####WITH TERMITES 
####Mature forest  0.4270338/0.2818414=1.52
####Regenerating forest 0.3871475/0.3293527=1.18
####Open land  0.4483452/0.3811592=1.18

###Mean=1.52 , mean(1.52,1.18,1.18)

#####Forest type ratio
#####WITHOUT TERMITES 
######Open to mature for Litsea 0.5834800/0.3730903=1.56
######Open to mature for Casstanopsis 0.2768242/0.1920045=1.44


#####WITH TERMITES 
######Open to mature for Litsea 0.4483452/0.4270338=1.04
######Open to mature for Casstanopsis 0.3811592/0.2818414=1.35


str(mengsong_k_envi)

levels(mengsong_k_envi$Termite_A_P)
levels(mengsong_k_envi$Termite_A_P)<-c("Absence","Presence")
levels(mengsong_k_envi$Termite_A_P)

######January 28th 2020
####Considering plotting down and up for wsg loss being a function of mass loss

head(Mengsong_field)
dim(Mengsong_field)

###Create the termite presence absence binary data 0=absence, 1=for presence
str(Mengsong_field)


####Subset the 36 mo data
Mengsong_field_36mo<-Mengsong_field[Mengsong_field$Coll_No.x=="Coll_6",]
dim(Mengsong_field_36mo)


###Litsea 36 mo
litsea_36mo<-Mengsong_field_36mo[Mengsong_field_36mo$Species_full=="Litsea_cubeba",]
dim(litsea_36mo)

######
###Modeling mass as function wsg up and down
hist(litsea_36mo$ML_percent)
summary(litsea_36mo$ML_percent)
hist(sqrt(litsea_36mo$ML_percent))
hist(log(litsea_36mo$ML_percent))
hist((litsea_36mo$ML_percent)^(1/3))


#######Litsea
#######Subset data Litsea without termites
str(litsea_36mo)
levels(litsea_36mo$Termi.assum)
levels(litsea_36mo$Termi.assum)<-c("0","1")
levels(litsea_36mo$Termi.assum)

litsea_36mo_no_termites<-litsea_36mo[litsea_36mo$Termi.assum==0,]
hist(litsea_36mo_no_termites$ML_percent)
hist(log(litsea_36mo_no_termites$ML_percent))
hist((litsea_36mo_no_termites$ML_percent)^(1/3))
hist(sqrt(litsea_36mo_no_termites$ML_percent))
dim(litsea_36mo_no_termites)
###massloss
L_mass_no_termi_lme1<-lme(data=litsea_36mo_no_termites, fixed=ML_percent~Per_WSG_loss+Pos_soil, random=~1|PLOT)
summary(L_mass_no_termi_lme1)##
anova(L_mass_no_termi_lme1, type="marginal")
AIC(L_mass_no_termi_lme1)
plot(L_mass_no_termi_lme1, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##decreasing trend
r.squaredGLMM(L_mass_no_termi_lme1)###r^2 marg=%, r^2cond=%

###Deal with variance trend
L_mass_no_termi_lme1_2<-lme(data=litsea_36mo_no_termites, fixed=ML_percent~Per_WSG_loss+Pos_soil, random=~1|PLOT, weights=varIdent(form=~1|For_type))
summary(L_mass_no_termi_lme1_2)##
anova(L_mass_no_termi_lme1_2, type="marginal")
AIC(L_mass_no_termi_lme1_2)
plot(L_mass_no_termi_lme1_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)#decreasing trend
r.squaredGLMM(L_mass_no_termi_lme1_2)###r^2 marg=%, r^2cond=%


###Deal with variance trend
L_mass_no_termi_lme1_2<-lme(data=litsea_36mo_no_termites, fixed=ML_percent~Per_WSG_loss+Pos_soil, random=~1|PLOT, weights=varPower(form=~fitted(.)))
summary(L_mass_no_termi_lme1_2)##
anova(L_mass_no_termi_lme1_2, type="marginal")
AIC(L_mass_no_termi_lme1_2)
plot(L_mass_no_termi_lme1_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##good
r.squaredGLMM(L_mass_no_termi_lme1_2)###r^2 marg=%, r^2cond=%



###Interaction of WSG and Position core
L_mass_no_termi_lme2<-lme(data=litsea_36mo_no_termites,fixed=ML_percent~Per_WSG_loss*Pos_soil,random=~1|PLOT)
summary(L_mass_no_termi_lme2)
anova(L_mass_no_termi_lme2, type="marginal")
AIC(L_mass_no_termi_lme2)
plot(L_mass_no_termi_lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)
r.squaredGLMM(L_mass_no_termi_lme2)###r^2 marg=%, r^2cond=%

###Deal with variance trend
L_mass_no_termi_lme2_2<-lme(data=litsea_36mo_no_termites, fixed=ML_percent~Per_WSG_loss*Pos_soil, random=~1|PLOT, weights=varPower(form=~fitted(.)))
summary(L_mass_no_termi_lme1_2)##
anova(L_mass_no_termi_lme2_2, type="marginal")
AIC(L_mass_no_termi_lme2_2)
plot(L_mass_no_termi_lme2_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)
r.squaredGLMM(L_mass_no_termi_lme2_2)

###BEST Model no transformation L_mass_no_termi_lme1_2
summary(L_mass_no_termi_lme1_2)##
anova(L_mass_no_termi_lme1_2, type="marginal")
AIC(L_mass_no_termi_lme1_2)

r.squaredGLMM(L_mass_no_termi_lme1_2)###r^2 marg=%, r^2cond=%


#####Model for k

hist(litsea_36mo_no_termites$k_ind)
hist(log(litsea_36mo_no_termites$k_ind))
hist((litsea_36mo_no_termites$k_ind)^(1/3))

L_k_no_termi_lme1<-lme(data=litsea_36mo_no_termites,fixed=(k_ind)^(1/3)~Per_WSG_loss+Pos_soil,random=~1|PLOT)
summary(L_k_no_termi_lme1)##adj r^2 ==00 %

anova(L_k_no_termi_lme1, type="marginal")
AIC(L_k_no_termi_lme1)
plot(L_k_no_termi_lme1, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###decreasing trend
r.squaredGLMM(L_k_no_termi_lme1)###r^2 marg=%, r^2cond=%


###Dealing with trend

L_k_no_termi_lme1_2<-lme(data=litsea_36mo_no_termites,fixed=(k_ind)^(1/3)~Per_WSG_loss+Pos_soil,random=~1|PLOT, weights=varPower(form=~fitted(.)))
summary(L_k_no_termi_lme1_2)##

anova(L_k_no_termi_lme1_2, type="marginal")
AIC(L_k_no_termi_lme1_2)
plot(L_k_no_termi_lme1_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###better
r.squaredGLMM(L_k_no_termi_lme1_2)###r^2 marg=%, r^2cond=%



###Interaction of WSG and Position core
L_k_no_termi_lme2<-lme(data=litsea_36mo_no_termites, fixed=(k_ind)^(1/3)~Per_WSG_loss*Pos_soil, random=~1|PLOT)
summary(L_k_no_termi_lme2)##adj r^2 
anova(L_k_no_termi_lme2, type="marginal")
AIC(L_k_no_termi_lme2)
plot(L_k_no_termi_lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###Better
r.squaredGLMM(L_k_no_termi_lme2)###r^2 marg=%, r^2cond=%

###Dealing with trend
L_k_no_termi_lme2_2<-lme(data=litsea_36mo_no_termites, fixed=(k_ind)^(1/3)~Per_WSG_loss*Pos_soil, random=~1|PLOT, weights=varPower(form=~fitted(.)))
summary(L_k_no_termi_lme2_2)##adj r^2 
anova(L_k_no_termi_lme2_2, type="marginal")
AIC(L_k_no_termi_lme2_2)
plot(L_k_no_termi_lme2_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###Better
r.squaredGLMM(L_k_no_termi_lme2_2)###r^2 marg=%, r^2cond=%


###BEST MODEL cubic transformation L_k_no_termi_lme1_2
summary(L_k_no_termi_lme1_2)##
anova(L_k_no_termi_lme1_2, type="marginal")
r.squaredGLMM(L_k_no_termi_lme1_2)###r^2 marg=%, r^2cond=%


#####Subset data Litsea with termites
#####Model for mass loss
litsea_36mo_termites<-litsea_36mo[litsea_36mo$Termi.assum==1,]
hist(litsea_36mo_termites$ML_percent)
hist(log(litsea_36mo_termites$ML_percent))
hist(sqrt(litsea_36mo_termites$ML_percent))
hist((litsea_36mo_termites$ML_percent)^(1/3))

L_mass_termi_lme1<-lme(data=litsea_36mo_termites, fixed=(ML_percent)^(1/3)~Per_WSG_loss+Pos_soil, random=~1|PLOT,na.action=na.omit)
summary(L_mass_termi_lme1)##
anova(L_mass_termi_lme1, type="marginal")
AIC(L_mass_termi_lme1)
plot(L_mass_termi_lme1, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###decreasing trend
r.squaredGLMM(L_mass_termi_lme1)###r^2 marg=%, r^2cond=%


### Dealing with variance
L_mass_termi_lme1_2<-lme(data=litsea_36mo_termites, fixed=(ML_percent)^(1/3)~Per_WSG_loss+Pos_soil, random=~1|PLOT,na.action=na.omit,weights=varPower(form=~fitted(.)))
summary(L_mass_termi_lme1_2)##
anova(L_mass_termi_lme1_2, type="marginal")
AIC(L_mass_termi_lme1_2)
plot(L_mass_termi_lme1_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###decreasing trend
r.squaredGLMM(L_mass_termi_lme1_2)###r^2 marg=%, r^2cond=%


### Dealing with variance
L_mass_termi_lme1_3<-lme(data=litsea_36mo_termites, fixed=(ML_percent)^(1/3)~Per_WSG_loss+Pos_soil, random=~1|PLOT,na.action=na.omit,weights=varExp(form=~fitted(.)))
summary(L_mass_termi_lme1_3)##
anova(L_mass_termi_lme1_3, type="marginal")
AIC(L_mass_termi_lme1_3)
plot(L_mass_termi_lme1_3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)###
r.squaredGLMM(L_mass_termi_lme1_3)###r^2 marg=%, r^2cond=%


## Dealing with variance
L_mass_termi_lme1_4<-lme(data=litsea_36mo_termites, fixed=(ML_percent)^(1/3)~Per_WSG_loss+Pos_soil, random=~1|PLOT,na.action=na.omit,weights=varIdent(form=~1|For_type))
summary(L_mass_termi_lme1_4)##
anova(L_mass_termi_lme1_4, type="marginal")
AIC(L_mass_termi_lme1_4)
plot(L_mass_termi_lme1_4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##worsen trend
r.squaredGLMM(L_mass_termi_lme1_4)###r^2 marg=%, r^2cond=%


###BEST MODEL ###  L_mass_termi_lme1_2
summary(L_mass_termi_lme1_2)##
anova(L_mass_termi_lme1_2, type="marginal")
r.squaredGLMM(L_mass_termi_lme1_2)###r^2 marg=%, r^2cond=%


###Interaction of WSG and Position core
L_mass_termi_lme2<-lme(data=litsea_36mo_termites, fixed=(ML_percent)^(1/3)~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(L_mass_termi_lme2)

anova(L_mass_termi_lme2, type="marginal")
AIC(L_mass_termi_lme2)
plot(L_mass_termi_lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Increasing trend
r.squaredGLMM(L_mass_termi_lme2)###r^2 marg=%, r^2cond=%


###Dealing with variance trend
L_mass_termi_lme2_2<-lme(data=litsea_36mo_termites, fixed=(ML_percent)^(1/3)~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varPower(form=~fitted(.)))
summary(L_mass_termi_lme2_2)

anova(L_mass_termi_lme2_2, type="marginal")
AIC(L_mass_termi_lme2_2)
plot(L_mass_termi_lme2_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Better
r.squaredGLMM(L_mass_termi_lme2_2)###r^2 marg=%, r^2cond=%



###Dealing with variance trend and simplify model by removing interaction
L_mass_termi_lme2_3<-lme(data=litsea_36mo_termites, fixed=(ML_percent)^(1/3)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varPower(form=~fitted(.)))
summary(L_mass_termi_lme2_3)

anova(L_mass_termi_lme2_3, type="marginal")
AIC(L_mass_termi_lme2_3)
plot(L_mass_termi_lme2_3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Better
r.squaredGLMM(L_mass_termi_lme2_3)###r^2 marg=%, r^2cond=%

###Best model cubic transformation L_mass_termi_lme2_3

summary(L_mass_termi_lme2_3)
anova(L_mass_termi_lme2_3, type="marginal")
r.squaredGLMM(L_mass_termi_lme2_3)###r^2 marg=%, r^2cond=%



#####Model for k

hist(litsea_36mo_termites$k_ind)
hist(log(litsea_36mo_termites$k_ind))
hist((litsea_36mo_termites$k_ind)^(1/3))

L_k_termi_lme1<-lme(data=litsea_36mo_termites, fixed=log(k_ind)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(L_k_termi_lme1)##adj r^2 
anova(L_k_termi_lme1, type="marginal")
AIC(L_k_termi_lme1)
plot(L_k_termi_lme1, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## increasing trend

###Dealing with trend

L_k_termi_lme1_2<-lme(data=litsea_36mo_termites, fixed=log(k_ind)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varPower(form=~fitted(.)))
summary(L_k_termi_lme1_2)##adj r^2 
anova(L_k_termi_lme1_2, type="marginal")
AIC(L_k_termi_lme1_2)
plot(L_k_termi_lme1_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Increasing trend

###Dealing with trend

L_k_termi_lme1_3<-lme(data=litsea_36mo_termites, fixed=log(k_ind)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varExp(form=~fitted(.)))
summary(L_k_termi_lme1_3)##adj r^2 
anova(L_k_termi_lme1_3, type="marginal")
AIC(L_k_termi_lme1_3)
plot(L_k_termi_lme1_3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Better
r.squaredGLMM(L_k_termi_lme1_3)

####BEST MODEL L_k_termi_lme1_3
summary(L_k_termi_lme1_3)##adj r^2 
anova(L_k_termi_lme1_3, type="marginal")
r.squaredGLMM(L_k_termi_lme1_3)


####With interaction
L_k_termi_lme2<-lme(data=litsea_36mo_termites, fixed=log(k_ind)~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(L_k_termi_lme2)##adj r^2
anova(L_k_termi_lme2, type="marginal")
AIC(L_k_termi_lme2)
plot(L_k_termi_lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## variance increasing 


###Dealing with trend
L_k_termi_lme2_2<-lme(data=litsea_36mo_termites, fixed=log(k_ind)~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varPower(form=~fitted(.)))
summary(L_k_termi_lme2_2)##adj r^2
anova(L_k_termi_lme2_2, type="marginal")
AIC(L_k_termi_lme2_2)
plot(L_k_termi_lme2_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slightly increasing trend


###Dealing with trend
L_k_termi_lme2_3<-lme(data=litsea_36mo_termites, fixed=log(k_ind)~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varExp(form=~fitted(.)))
summary(L_k_termi_lme2_3)##adj r^2
anova(L_k_termi_lme2_3, type="marginal")
AIC(L_k_termi_lme2_3)
plot(L_k_termi_lme2_3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Better 

###BEST MODEL is L_k_termi_lme1_3
summary(L_k_termi_lme1_3)
r.squaredGLMM(L_k_termi_lme1_3)### r^2 marg=0.65 % r^2cond=71.33%


##Graph
litsea_36mo$Termi.assum<-as.factor(litsea_36mo$Termi.assum)
levels(litsea_36mo$Termi.assum)
levels(litsea_36mo$Termi.assum)<-c("Absence","Presence")
levels(litsea_36mo$Termi.assum)

litsea_mass_WSG_loss_graph<-ggplot(data=litsea_36mo,aes(x=Per_WSG_loss, y=ML_percent, color=Termi.assum, shape= Pos_soil))+
  geom_point(data=litsea_36mo, stat="identity",size=3)+
  coord_fixed(ratio = 1)+
  geom_abline(aes(intercept=0, slope=1), lty=2)+
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100))+
  xlab(" Percent wood specific gravity loss") +
  ylab("Percent mass loss")+ 
  labs(color="Termites status")+
  labs(shape="Position of wood core")+
  labs(title="A-Litsea cubeba")+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_line(colour = "black"),
        axis.ticks.y=element_line(colour = "black"),
        axis.text = element_text(colour = "black", size =25),
        axis.title = element_text(colour = "black", size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        legend.position = c("left"),
        legend.text=element_text(size=25),
        legend.title = element_text(size=25),
        legend.background = element_blank(),
        plot.title = element_text(size=25, face="italic"),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

litsea_mass_WSG_loss_graph


#########Catsanopsis data
###Castanopsis 36mo
castanopsis_36mo<-Mengsong_field_36mo[Mengsong_field_36mo$Species_full=="Castanopsis_mekongensis",]
dim(castanopsis_36mo)

###Modeling k_ind as function wsg up and down
hist(castanopsis_36mo$k_ind)
hist(log(castanopsis_36mo$k_ind))



#####Castanopsis

#######Subset data castanopsis without termites
#####Model for k
cast_36mo_no_termites<-castanopsis_36mo[castanopsis_36mo$Termi.assum==0,]

hist(cast_36mo_no_termites$k_ind)
hist(log(cast_36mo_no_termites$k_ind))
#hist((casta_36mo_no_termites$k_ind)^(1/3))


#####Subset data with termites
#####Model for mass loss
cast_36mo_termites<-castanopsis_36mo[castanopsis_36mo$Termi.assum==1,]
hist(cast_36mo_termites$ML_percent)
hist(log(cast_36mo_termites$ML_percent))
hist(sqrt(cast_36mo_termites$ML_percent))
hist((cast_36mo_termites$ML_percent)^(1/3))


#######Subset data castanopsis without termites
#####Model for k
str(castanopsis_36mo)
levels(castanopsis_36mo$Termi.assum)
levels(castanopsis_36mo$Termi.assum)<-c("0","1")
levels(castanopsis_36mo$Termi.assum)

cast_36mo_no_termites<-castanopsis_36mo[castanopsis_36mo$Termi.assum==0,]

hist(cast_36mo_no_termites$k_ind)
hist(log(cast_36mo_no_termites$k_ind))
#hist((casta_36mo_no_termites$k_ind)^(1/3))

C_k_no_termi_lme1<-lme(data=cast_36mo_no_termites, fixed=k_ind~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_k_no_termi_lme1)##adj r^2 ==
anova(C_k_no_termi_lme1, type="marginal")
AIC(C_k_no_termi_lme1)
plot(C_k_no_termi_lme1, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slight increase in the trend

###Dealing with trend

C_k_no_termi_lme1_2<-lme(data=cast_36mo_no_termites, fixed=k_ind~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varPower(form=~fitted(.)))
summary(C_k_no_termi_lme1_2)##adj r^2 ==
anova(C_k_no_termi_lme1_2, type="marginal")
AIC(C_k_no_termi_lme1_2)
plot(C_k_no_termi_lme1_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slight decrease in the trend


###Dealing with trend

C_k_no_termi_lme1_3<-lme(data=cast_36mo_no_termites, fixed=k_ind~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varExp(form=~fitted(.)))
summary(C_k_no_termi_lme1_3)##adj r^2 ==
anova(C_k_no_termi_lme1_3, type="marginal")
AIC(C_k_no_termi_lme1_3)
plot(C_k_no_termi_lme1_3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slight decrease in the trend


####BEST MODEL no transformation C_k_no_termi_lme1_3
r.squaredGLMM(C_k_no_termi_lme1_3)#### r^2 marg=18.81% r^2cond=79.56%


###With interaction with position core
C_k_no_termi_lme2<-lme(data=cast_36mo_no_termites,fixed=k_ind~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_k_no_termi_lme2)##adj r^2 

anova(C_k_no_termi_lme2, type="marginal")
AIC(C_k_no_termi_lme2)
plot(C_k_no_termi_lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## not bad
r.squaredGLMM(C_k_no_termi_lme2)###

### No position
C_k_no_termi_lme3<-lme(data=cast_36mo_no_termites,fixed=k_ind~Per_WSG_loss,random=~1|PLOT,na.action=na.omit)
summary(C_k_no_termi_lme3)##adj r^2 ==%

anova(C_k_no_termi_lme3, type="marginal")
AIC(C_k_no_termi_lme3)
plot(C_k_no_termi_lme3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## increasing trend

###Dealing with trend
C_k_no_termi_lme3_2<-lme(data=cast_36mo_no_termites,fixed=k_ind~Per_WSG_loss,random=~1|PLOT,na.action=na.omit,weights=varPower(form=~fitted(.)))
summary(C_k_no_termi_lme3_2)##adj r^2 ==%

anova(C_k_no_termi_lme3_2, type="marginal")
AIC(C_k_no_termi_lme3_2)
plot(C_k_no_termi_lme3_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##sligtly decreasing trend 

r.squaredGLMM(C_k_no_termi_lme3_2)###r^2 marg=10.63% r^2cond=45.50%

####BEST MODEL C_k_no_termi_lme1_3
summary(C_k_no_termi_lme1_3)
anova(C_k_no_termi_lme1_3, type="marginal")

r.squaredGLMM(C_k_no_termi_lme1_3)##r^2 marg=18.81% r^2cond=79.61%


###Model for mass loss
hist(cast_36mo_no_termites$ML_percent)
hist(log(cast_36mo_no_termites$ML_percent))
hist((cast_36mo_no_termites$ML_percent)^(1/3))
hist(sqrt(cast_36mo_no_termites$ML_percent))

C_mass_no_termi_lme1<-lme(data=cast_36mo_no_termites, fixed=ML_percent~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_mass_no_termi_lme1)##adj r^2 == %
anova(C_mass_no_termi_lme1, type="marginal")
AIC(C_mass_no_termi_lme1)
plot(C_mass_no_termi_lme1, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Decreasing trend 

r.squaredGLMM(C_mass_no_termi_lme1)###r^2 marg= % r^2cond= %

#Dealing with trend
C_mass_no_termi_lme1_2<-lme(data=cast_36mo_no_termites, fixed=ML_percent~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit, weights=varPower(form=~fitted(.)))
summary(C_mass_no_termi_lme1_2)##adj r^2 == %
anova(C_mass_no_termi_lme1_2, type="marginal")
AIC(C_mass_no_termi_lme1_2)
plot(C_mass_no_termi_lme1_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Decreasing trend 

r.squaredGLMM(C_mass_no_termi_lme1_2)###r^2 marg= % r^2cond= %



#Dealing with trend
C_mass_no_termi_lme1_3<-lme(data=cast_36mo_no_termites, fixed=ML_percent~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit, weights=varExp(form=~fitted(.)))
summary(C_mass_no_termi_lme1_3)##adj r^2 == %
anova(C_mass_no_termi_lme1_3, type="marginal")
AIC(C_mass_no_termi_lme1_3)
plot(C_mass_no_termi_lme1_3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## slight decreasing trend 

r.squaredGLMM(C_mass_no_termi_lme1_3)###r^2 marg=7.99% r^2cond=32.63%

r.squaredLR(C_mass_no_termi_lme1_3)

#Dealing with trend
C_mass_no_termi_lme1_4<-lme(data=cast_36mo_no_termites, fixed=ML_percent~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit, weights=varIdent(form=~1|For_type))
summary(C_mass_no_termi_lme1_4)##adj r^2 == %
anova(C_mass_no_termi_lme1_4, type="marginal")
AIC(C_mass_no_termi_lme1_4)
plot(C_mass_no_termi_lme1_4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Decreasing trend 

r.squaredGLMM(C_mass_no_termi_lme1_4)###r^2 marg=% r^2cond=%

###Remove Position core
C_mass_no_termi_lme2<-lme(data=cast_36mo_no_termites, fixed=ML_percent~Per_WSG_loss,random=~1|PLOT,na.action=na.omit)
summary(C_mass_no_termi_lme2)##adj r^2 == %
anova(C_mass_no_termi_lme2, type="marginal")
AIC(C_mass_no_termi_lme2)
plot(C_mass_no_termi_lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Decreasing trend 

r.squaredGLMM(C_mass_no_termi_lme2)###r^2 marg=% r^2cond=%

###, weights=varIdent(form=~1|For_type)

###BEST MODEL C_mass_no_termi_lme1_3
summary(C_mass_no_termi_lme1_3)##adj r^2 == %
anova(C_mass_no_termi_lme1_3, type="marginal")
r.squaredGLMM(C_mass_no_termi_lme1_3)###r^2 marg=7.99% r^2cond=32.63%


###With interaction with position core
C_mass_no_termi_lme3<-lme(data=cast_36mo_no_termites, fixed=ML_percent~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_mass_no_termi_lme3)##adj r^2 == %
anova(C_mass_no_termi_lme3, type="marginal")
AIC(C_mass_no_termi_lme3)
plot(C_mass_no_termi_lme3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Decreasing trend 


### No position core
C_mass_no_termi_lme5<-lme(data=cast_36mo_no_termites, fixed=ML_percent~Per_WSG_loss,random=~1|PLOT,na.action=na.omit)
summary(C_mass_no_termi_lme5)##adj r^2 == %

anova(C_mass_no_termi_lme5, type="marginal")
AIC(C_mass_no_termi_lme5)
plot(C_mass_no_termi_lme5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## Decreasing trend 





#####Subset data with termites
#####Model for mass loss
cast_36mo_termites<-castanopsis_36mo[castanopsis_36mo$Termi.assum==1,]
hist(cast_36mo_termites$ML_percent)
hist(log(cast_36mo_termites$ML_percent))
hist(sqrt(cast_36mo_termites$ML_percent))
hist((cast_36mo_termites$ML_percent)^(1/3))





###squarred root
C_mass_termi_lme<-lme(data=cast_36mo_termites, fixed=sqrt(ML_percent)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_mass_termi_lme)##adj r^2 ==00 %


anova(C_mass_termi_lme, type="marginal")
AIC(C_mass_termi_lme)
plot(C_mass_termi_lme, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## no trend 

r.squaredGLMM(C_mass_termi_lme)###r^2 marg=0% r^2cond=00%


###No transformation
C_mass_termi_lme2<-lme(data=cast_36mo_termites, fixed=ML_percent~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_mass_termi_lme2)##adj r^2 ==00 %
anova(C_mass_termi_lme2, type="marginal")
AIC(C_mass_termi_lme2)
plot(C_mass_termi_lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## increasing trend 

###No transformation but interaction

C_mass_termi_lme3<-lme(data=cast_36mo_termites, fixed=ML_percent~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_mass_termi_lme3)
anova(C_mass_termi_lme3, type="marginal")
AIC(C_mass_termi_lme3)
plot(C_mass_termi_lme3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##  increasing trend


###Log trannsformation
C_mass_termi_lme4<-lme(data=cast_36mo_termites, fixed=log(ML_percent)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_mass_termi_lme4)
anova(C_mass_termi_lme4, type="marginal")
AIC(C_mass_termi_lme4)
plot(C_mass_termi_lme4, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##  decreasing trend




###Interaction with WSG loss and POsition
C_mass_termi_lme5<-lme(data=cast_36mo_termites,fixed=log(ML_percent)~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_mass_termi_lme5)##adj r^2 ==00 %

anova(C_mass_termi_lme5, type="marginal")
AIC(C_mass_termi_lme5)
plot(C_mass_termi_lme5, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##  




####
###squarred root and interaction
C_mass_termi_lme6<-lme(data=cast_36mo_termites, fixed=sqrt(ML_percent)~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_mass_termi_lme6)


anova(C_mass_termi_lme6, type="marginal")
AIC(C_mass_termi_lme6)
plot(C_mass_termi_lme6, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## no trend 


####BEST MODEL squarred root transformation C_mass_termi_lme
summary(C_mass_termi_lme)
r.squaredGLMM(C_mass_termi_lme)###r^2 marg=0.70% r^2cond=67.57%



#####Model for k

hist(cast_36mo_termites$k_ind)
hist(log(cast_36mo_termites$k_ind))
hist((cast_36mo_termites$k_ind)^(1/3))

###log transformed
C_k_termi_lme<-lme(data=cast_36mo_termites, log(k_ind)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_k_termi_lme)##adj r^2 =00 %
anova(C_k_termi_lme, type="marginal")
AIC(C_k_termi_lme)
plot(C_k_termi_lme, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## increasing trend 


###Dealing with trend 
C_k_termi_lme1_2<-lme(data=cast_36mo_termites, log(k_ind)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit,weights=varExp(form=~fitted(.)))
summary(C_k_termi_lme1_2)##adj r^2 =00 %
anova(C_k_termi_lme1_2, type="marginal")
AIC(C_k_termi_lme1_2)
plot(C_k_termi_lme1_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## increasing trend 
r.squaredGLMM(C_k_termi_lme1_2)


####With interaction
C_k_termi_lme2<-lme(data=cast_36mo_termites,fixed=log(k_ind)~Per_WSG_loss*Pos_soil, random=~1|PLOT,na.action=na.omit)
summary(C_k_termi_lme2)##adj r^2 ==00 %
anova(C_k_termi_lme2, type="marginal")
AIC(C_k_termi_lme2)
plot(C_k_termi_lme2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)## increasing trend

###Dealing with trend
C_k_termi_lme2_2<-lme(data=cast_36mo_termites,fixed=log(k_ind)~Per_WSG_loss*Pos_soil, random=~1|PLOT,na.action=na.omit,weights=varExp(form=~fitted(.)))
summary(C_k_termi_lme2_2)##adj r^2 ==00 %
anova(C_k_termi_lme2_2, type="marginal")
AIC(C_k_termi_lme2_2)
plot(C_k_termi_lme2_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##



###^(1/3)
C_k_termi_lme3<-lme(data=cast_36mo_termites, fixed=(k_ind)^(1/3)~Per_WSG_loss+Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_k_termi_lme3)
anova(C_k_termi_lme3, type="marginal")
AIC(C_k_termi_lme3)
plot(C_k_termi_lme3, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Increasing trend

####With interaction
C_k_termi_lme3_2<-lme(data=cast_36mo_termites, fixed=(k_ind)^(1/3)~Per_WSG_loss*Pos_soil,random=~1|PLOT,na.action=na.omit)
summary(C_k_termi_lme3_2)
anova(C_k_termi_lme3_2, type="marginal")
AIC(C_k_termi_lme3_2)
plot(C_k_termi_lme3_2, abs(resid(., type='pearson'))~fitted(.), type=c('p', 'r'), 
     abline=0)##Increasing trend




#######Best model log transformed C_k_termi_lme1_2
r.squaredGLMM(C_k_termi_lme1_2)### r^2 marg=0.5% r^2cond== 37.86



###Graph
str(castanopsis_36mo)

castanopsis_36mo$Termi.assum<-as.factor(castanopsis_36mo$Termi.assum)
levels(castanopsis_36mo$Termi.assum)
levels(castanopsis_36mo$Termi.assum)<-c("Absence","Presence")
levels(castanopsis_36mo$Termi.assum)

castanopsis_mass_WSG_loss_graph2<-ggplot(data=castanopsis_36mo,aes(x=Per_WSG_loss, y=ML_percent, color=Termi.assum, shape= Pos_soil))+
  geom_point(data=castanopsis_36mo, stat="identity",size=3)+
  geom_smooth(data=castanopsis_36mo[castanopsis_36mo$Termi.assum=="Absence",], aes(lty=Pos_soil), method="lm", formula=y~x, se = T)+
  coord_fixed(ratio = 1)+
  geom_abline(aes(intercept=0, slope=1), lty=2)+
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100))+
  xlab(" Percent wood specific gravity loss") +
  ylab("")+ #Percent mass loss
  labs(color="Termites status")+
  labs(shape="Position of wood core")+
  labs(title="B-Castanopsis mekongensis")+
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.x=element_line(colour = "black"),
        axis.ticks.y=element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size = 25),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.title = element_text(colour = "black", size = 25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        strip.text = element_text(colour = "black", size = 15),
        legend.position = c("top"),
        legend.text=element_text(size=12),
        legend.title = element_text(labs('Position'),size=12.5),
        legend.background = element_blank(),
        plot.title = element_text(size=25, face="italic"),
        legend.key.width=unit(0.4,"cm"),
        legend.key.height=unit(0.4,"cm")) 

castanopsis_mass_WSG_loss_graph2



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(litsea_mass_WSG_loss_graph)

library(ggplot2)
library(gridExtra)
library(reshape)
library(grid)


### Figure 4 Mass loss as a function of wood secific gravity loss
tiff(filename="Figure 4 mass loss against wsg loss up and down with significant relationships2.tiff", res = 400, width=6000, height=5000, compression = "lzw")

grid.arrange(arrangeGrob(litsea_mass_WSG_loss_graph + theme(legend.position="none", plot.margin = unit( c(1,0,0.5,0.5) , units = "lines")),
                         castanopsis_mass_WSG_loss_graph2 + theme(legend.position="none", plot.margin = unit(c(1,0,0.5,0.5) , units = "lines")),
                         ncol=2),mylegend)

dev.off()



###Compute mean percent mass loss
library(dplyr)
library(tidyr)
str(mass.loss.36mo.down)

mass.loss.36mo.down[mass.loss.36mo.down$Species.x=="Litsea_cubeba",]
summary(mass.loss.36mo.down[mass.loss.36mo.down$Species.x=="Litsea_cubeba",])
summary(mass.loss.36mo.down[mass.loss.36mo.down$Species.x=="Castanopsis_mekongensis",])

##Using piping method in R piping is done with %>%
mass.loss.36mo.down_Means<-mass.loss.36mo.down%>%
  group_by(Species.x)%>%
  summarise(mean_ML_percent=mean(ML_percent),
            sd_ML_percent=sd(ML_percent))#,
mass.loss.36mo.down_Means

# A tibble: 2 x 3
#  Species.x                    mean_ML_percent sd_ML_percent
#   <fct>                                <dbl>         <dbl>
#  1 Castanopsis_mekongensis            52.4          17.2
#  2 Litsea_cubeba                      66.5          16.7


## Means per species per forest type
mass.loss.36mo.down_Means2<-mass.loss.36mo.down%>%
  group_by(Species.x,For_type)%>%
  summarise(mean_ML_percent=mean(ML_percent),
            sd_ML_percent=sd(ML_percent))#,
mass.loss.36mo.down_Means2




###Resutls
#> mass.loss.36mo.down_Means2
# A tibble: 6 x 4
# Groups:   Species.x [2]
#  Species.x               For_type            mean_ML_percent sd_ML_percent
#    <fct>                   <fct>                         <dbl>         <dbl>
#  1 Castanopsis_mekongensis Mature forest                  51.1          17.4
#  2 Castanopsis_mekongensis Regenerating Forest            49.6          16.0
#  3 Castanopsis_mekongensis Open land                      62.5          16.6
#  4 Litsea_cubeba           Mature forest                  62.0          17.6
#  5 Litsea_cubeba           Regenerating Forest            67.4          15.4
#  6 Litsea_cubeba           Open land                      78.9          12.0



###SUPPLEMENTARY MATERIALS FIGURES on MICROCLIMATE
#########Climate and microclimate_Plot########
p_temp <-read.csv("Temperature_landscape_government station.csv", header=T, sep=",")
str(p_temp)
head(p_temp)
levels(p_temp$Station)
p_temp$Station <- factor(p_temp$Station, levels = c("Mengsong","Mature forest","Regenerating forest","Open land"))
levels(p_temp$Station)<- c("Climatic station","Mature forest","Regenerating forest","Open land")
aggregate(cbind(Temp_Max) ~ SEASON + Station, data = p_temp, function(x) c(M = mean(x), SD = sd(x)))

Rdate1 <- strptime(as.character(p_temp$Date),"%d-%B-%y")

p_temp1 <- data.frame(p_temp,Rdate1)
head(p_temp1)

library(ggplot2)
library(gridExtra)
library(reshape)
wet = subset (p_temp1, SEASON == "WET" )
ww <- ggplot(data = wet, aes(y=Temp_Max, x=Rdate1, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('blue','red','purple','green'))+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position = "none")

dry = subset (p_temp1, SEASON == "DRY" )
dd <- ggplot(data = dry, aes(y=Temp_Max, x=Rdate1, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('blue','red','purple','green'))+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= c('bottom'),
        legend.text=element_text(size=10),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(2,"cm"),
        legend.key.height=unit(0.4,"cm"))

library(grid)
library(ggplot2)
library(gridExtra)
library(reshape)
library(scales)


###Figure S1 A panel
tiff(filename="Figure S1 A Temperature_landscape government station.tiff", res = 600, width=6000, height=5000, compression = "lzw")
grid.arrange(arrangeGrob(ww + theme(legend.position="none",
                                    plot.margin = unit( c(1,1,0,0) , units = "lines")), 
                         dd + theme(legend.position="bottom",
                                    plot.margin = unit( c(0,1,1,0) , units = "lines")),
                         nrow=2,
                         left = textGrob("Maximum temperature (°C)", 
                                         rot = 90, vjust = 1)))#,
dev.off()


########Temperature landscape######
p_temp <-read.csv("Temp_Landscape_only.csv", header=T, sep=",")
str(p_temp)
head(p_temp)
levels(p_temp$Station)
p_temp$Station <- factor(p_temp$Station, levels = c("Mature forest","Regenerating forest","Open land"))
levels(p_temp$Station)<- c("Mature forest","Regenerating forest","Open land")
aggregate(cbind(Temp_Max) ~ SEASON + Station, data = p_temp, function(x) c(M = mean(x), SD = sd(x)))

Rdate1 <- strptime(as.character(p_temp$Date),"%d-%B-%y")

p_temp1 <- data.frame(p_temp,Rdate1)
head(p_temp1)
library(grid)
library(ggplot2)
library(gridExtra)
library(reshape)
library(scales)
wet1 = subset (p_temp1, Year=="One" & SEASON == "WET")
ww1 <- ggplot(data = wet1, aes(y=Temp_Max, x=Rdate1, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        axis.title.x = element_text(color="black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position = "none")
wet2 = subset (p_temp1, Year=="Two" & SEASON == "WET")
ww2 <- ggplot(data = wet2, aes(y=Temp_Max, x=Rdate1, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        axis.title.x = element_text(color="black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position = "none")

dry1 = subset (p_temp1, Year== "One" & SEASON == "DRY" )
dd1 <- ggplot(data = dry1, aes(y=Temp_Max, x=Rdate1, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        axis.title.x = element_text(color="black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position = "none")
dry2 = subset (p_temp1, Year== "Two" & SEASON == "DRY" )
dd2 <- ggplot(data = dry2, aes(y=Temp_Max, x=Rdate1, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        axis.title.x = element_text(color="black", size=8),
        legend.position= c('bottom'),
        legend.text=element_text(size=10),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(2,"cm"),
        legend.key.height=unit(0.4,"cm"))

###Figure S1 B panel

tiff(filename="Figure S1 B panel Temperature_landscape_only.tiff", res = 600, width=6000, height=5000, compression = "lzw")
grid.arrange(arrangeGrob(ww1 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         ww2 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         dd1 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         dd2 + theme(legend.position="bottom",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         nrow=4,
                         left = textGrob("Maximum temperature (°C)", 
                                         rot = 90, vjust = 1)))#,
dev.off()

##Soil moisture####
soil_moist <-read.csv("Soil moisture and RH_Landscape.csv", header=T, sep=",")
head(soil_moist)
str(soil_moist)

levels(soil_moist$Station)
soil_moist$Station <- factor(soil_moist$Station, levels = c("Mature forest","Regenerating forest","Open land"))

DATE <- strptime(as.character(soil_moist$Date),"%d-%B-%y")
soil_moist1 <- data.frame(soil_moist, DATE)
head(soil_moist1)
str(soil_moist1)
wet1 = subset (soil_moist1, Year=="One" & SEASON == "WET" )
ww1 <- ggplot(data = wet1, aes(y=Avg_Soil_Water_Content, x=DATE, colour=Station)) + ylim(0,0.4) +
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")
wet2 = subset (soil_moist1, Year=="Two" & SEASON == "WET" )
ww2 <- ggplot(data = wet2, aes(y=Avg_Soil_Water_Content, x=DATE, colour=Station)) + ylim(0,0.4) +
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")

dry1 = subset (soil_moist1, Year=="One" & SEASON == "DRY" )
dd1 <- ggplot(data = dry1, aes(y=Avg_Soil_Water_Content, x=DATE, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")
dry2 = subset (soil_moist1, Year=="Two" & SEASON == "DRY" )
dd2 <- ggplot(data = dry2, aes(y=Avg_Soil_Water_Content, x=DATE, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= c('bottom'),
        legend.text=element_text(size=10),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(2,"cm"),
        legend.key.height=unit(0.4,"cm"))        


####Figure S2 Moisture 
tiff(filename="Figure S2 Soil moisture_landscape.tiff", res = 600, width=6000, height=5000, compression = "lzw")
grid.arrange(arrangeGrob(ww1 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         ww2 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         dd1 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         dd2 + theme(legend.position="bottom",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         nrow=4,
                         left = textGrob(expression(Soil~moisture~"("*m^"3"*" "*m^"-3"*")"), 
                                         rot = 90, vjust = 1)))
dev.off()

#######relative humidity####

ww1 <- ggplot(data = wet1, aes(y=Hum_Avg, x=DATE, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")
ww2 <- ggplot(data = wet2, aes(y=Hum_Avg, x=DATE, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")

dd1 <- ggplot(data = dry1, aes(y=Hum_Avg, x=DATE, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")

dd2 <- ggplot(data = dry2, aes(y=Hum_Avg, x=DATE, colour=Station)) + 
  geom_line(aes(group=Station), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= c('bottom'),
        legend.text=element_text(size=10),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(2,"cm"),
        legend.key.height=unit(0.4,"cm"))

####Figure S3 Relative humidity 
tiff(filename="Figure S3 Relative humidity_landscape.tiff", res = 600, width=6000, height=5000, compression = "lzw")
grid.arrange(arrangeGrob(ww1 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         ww2 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         dd1 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         dd2 + theme(legend.position="bottom",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         nrow=4,
                         left = textGrob("Relative humidity (%)", 
                                         rot = 90, vjust = 1)))
dev.off()




#####PAR Photosynthetical active radiation ########
par <-read.csv("PAR_median_plot.csv", header=T, sep=",")
head(par)
str(par)

levels(par$FOR_TYPE)
par$FOR_TYPE <- factor(par$FOR_TYPE, levels = c("Mature forest","Regenerating forest","Open land"))

DATE <- strptime(as.character(par$Date),"%d-%B-%y")
par1 <- data.frame(par, DATE)
head(par1)
str(par1)
wet1 = subset (par1, Year=="One" & SEASON == "WET" )
ww1 <- ggplot(data = wet1, aes(y=median_PAR, x=DATE, colour=FOR_TYPE)) + ylim(0,2200) +
  geom_line(aes(group=FOR_TYPE), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")
wet2 = subset (par1, Year=="Two" & SEASON == "WET" )
ww2 <- ggplot(data = wet2, aes(y=median_PAR, x=DATE, colour=FOR_TYPE)) + ylim(0,2200) +
  geom_line(aes(group=FOR_TYPE), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")

dry1 = subset (par1, Year=="One" & SEASON == "DRY" )
dd1 <- ggplot(data = dry1, aes(y=median_PAR, x=DATE, colour=FOR_TYPE)) + 
  geom_line(aes(group=FOR_TYPE), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= "none")
dry2 = subset (par1, Year=="Two" & SEASON == "DRY" )
dd2 <- ggplot(data = dry2, aes(y=median_PAR, x=DATE, colour=FOR_TYPE)) + 
  geom_line(aes(group=FOR_TYPE), size = 0.8) +
  scale_colour_manual(values = c('red','purple','green'))+
  scale_x_datetime(date_labels = "%b-%Y")+
  xlab(" ") +
  ylab(" ") + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.text = element_text(colour = "black", size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.position= c('bottom'),
        legend.text=element_text(size=10),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.width=unit(2,"cm"),
        legend.key.height=unit(0.4,"cm"))        

####Figure S4 Photosynthetically active radiation
tiff(filename="Figure S4 PAR_PLOT.tiff", res = 600, width=6000, height=5000, compression = "lzw")
grid.arrange(arrangeGrob(ww1 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         ww2 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         dd1 + theme(legend.position="none",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         dd2 + theme(legend.position="bottom",
                                     plot.margin = unit( c(1,1,0,0) , units = "lines")),
                         nrow=4,
                         left = textGrob("PAR (uE)", 
                                         rot = 90, vjust = 1)))
dev.off()


##computing the mean and sd
#####Max temperature by season and by forest type
head(p_temp1)
str(p_temp1)
levels(p_temp1$Station)
levels(p_temp1$SEASON)
####Subset the stations needed [all except Climatic station]
p_temp2<-subset(p_temp1, Station!="Climatic station")
levels(p_temp2$Station)
p_temp2[] <- lapply(p_temp2, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are consistent with the subset dataframe
levels(p_temp2$Station)
levels(p_temp2$SEASON)

##Piping method
library(plyr)
library(dplyr)
library(tidyr)
summary(p_temp2)
##Using piping method in R piping is done with %>%
p_temp1_Means<-p_temp2%>%
  group_by(Station,SEASON)%>%
  summarise(mean_temp_max=mean(Temp_Max, na.rm = T),
            sd_temp_max=sd(Temp_Max, na.rm = T))#,

p_temp1_Means


#####Soil moisture season and by forest type
head(soil_moist)
str(soil_moist)
levels(soil_moist$Station)
soil_moist$Station <- factor(soil_moist$Station, levels = c("Mature forest","Regenerating forest","Open land"))
levels(soil_moist$SEASON)

##Using piping method in R piping is done with %>%
soil_moist_Means<-soil_moist%>%
  group_by(Station,SEASON)%>%
  summarise(mean_soil_moist=mean(Avg_Soil_Water_Content, na.rm = T),
            sd_soil_moist=sd(Avg_Soil_Water_Content, na.rm = T))#,

soil_moist_Means

###Relative humidity 
head(soil_moist)
str(soil_moist)
levels(soil_moist$Station)
soil_moist$Station <- factor(soil_moist$Station, levels = c("Mature forest","Regenerating forest","Open land"))
levels(soil_moist$SEASON)

##Using piping method in R piping is done with %>%
RH_avg_Means<-soil_moist%>%
  group_by(Station,SEASON)%>%
  summarise(mean_RH_avg=mean(Hum_Avg, na.rm = T),
            sd_RH_avg=sd(Hum_Avg, na.rm = T))#,

RH_avg_Means

#####PAR for landscape by season and by forest type
par <-read.csv("PAR_median_plot.csv", header=T, sep=",")
head(par)
str(par)

levels(par$FOR_TYPE)
par$FOR_TYPE <- factor(par$FOR_TYPE, levels = c("Mature forest","Regenerating forest","Open land"))

DATE <- strptime(as.character(par$Date),"%d-%B-%y")
par1 <- data.frame(par, DATE)
head(par1)
str(par1)
levels(par1$FOR_TYPE)
levels(p_temp1$SEASON)
####Subset the stations needed [all except Climatic station]
par2<-subset(par1, SEASON!="")
levels(par2$SEASON)
par2[] <- lapply(par2, function(x) if(is.factor(x)) factor(x) else x)##make sure the factor levels are consistent with the subset dataframe
levels(par2$SEASON)

##Piping method
library(plyr)
library(dplyr)
library(tidyr)
summary(par2)
##Using piping method in R piping is done with %>%
par2_Means<-par2%>%
  group_by(FOR_TYPE,SEASON)%>%
  summarise(mean_median_PAR=mean(median_PAR, na.rm = T),
            sd_median_PAR=sd(median_PAR, na.rm = T))#,

par2_Means


###combine all averages togeter in a dataframe
avg.table<-data.frame(p_temp1_Means,soil_moist_Means[,c(3:4)],RH_avg_Means[,c(3:4)],par2_Means[,c(3:4)])

avg.table
##Write the averages tablein csv
Average_envi_microclimate<-write.csv(avg.table, file="Average_envi_microclimate.csv")
####Table S1
