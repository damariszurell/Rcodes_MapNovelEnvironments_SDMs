# 'Predicting to new environments: tools for visualising model behaviour and impacts on mapped distributions'
# by Damaris Zurell, Jane Elith and Boris Schröder.
# Appendix S1 - artifical species
# last update 2011-09-29

# example with simulated data:
# species' occurrence described by species' tolerance to temperature, woodland;
# three cases with different data coverage of two-dimensional environmental niche:
# (1) species niche entirely encompassed by data;
# (2) species niche truncated, i.e. portions of the niche are not represented in data;
# (3) species niche abuts environmental data, i.e. niche edge coincides with data limits;
# For all cases, different SDMs are estimated on training data, and fitted values are compared.
# Then, predictions are made for changing climate (=warmer temperature while land cover remains constant)
# and again fitted values for future climate are compared.

#**********************************************
# set working directory
setwd("...")

# load libraries
library(Design)
library(boot)
library(gam)
library(gbm)
source("brt.functions.R") #this is extra code provided in Elith et al. (2008) JAnimEcol 77:802-813
# note that the inflated.response() function will additionally require the package 'lhs' to be installed

#********************************************************
# some functions for creating species data, for evaluation and plotting
# (this keeps subsequent code for actual analysis much tidier)
# are contained in 'appendixS1_functions.r'

source('appendixS1_functions.r')
# contains following functions:
# 'species' for creating simulated species;
# 'myplot' for plotting fitted values along with true response curve
# 'inflated.response' for plotting inflated response curves
# 'eo.mask' for calculating environmental overlap

#****************************************************************************
# global variables:
minTemp=3
maxTemp=13
minWood=0
maxWood=70

#****************************************************************************
#assumed response surface (lrm)
temperature<-seq(minTemp,maxTemp,length=25)
woodland<-seq(minWood,maxWood,length=25)

dat=data.frame(expand.grid(woodlandCover=woodland,temp=temperature))
response=inv.logit(-170+40*dat$temp-2.5*dat$temp^2+.35*dat$woodlandCover)
windows()
wireframe(response~dat$woodlandCover*dat$temp,,
          scales=list(arrows=F,tck=.6,distance=.7,z=list(at=c(0,.5,1),labels=c("0.0","","1.0"),cex=2),
          x=list(at=c(10,30,50,70),cex=2),y=list(cex=2),col="black"),
          zlim=c(0,1),zlab=list('True response',rot=94,cex=2),
          xlab=list("Woodland cover [%]",rot=33,cex=2),ylab=list("Temperature [°C]",rot=-26,cex=2),
          par.settings = list(axis.line = list(col = "transparent")), screen=list(z=50,x=-70,y=0),
          alpha.regions=.7)

rm(list=c("woodland","temperature","dat","response"))

#*************************************************************************
#*************************************************************************
#*************************************************************************
# case 1: species niche entirely encompassed by data;

# training data + future
CurrentTemperature=runif(1000,min=minTemp,max=maxTemp)
woodlandCover=runif(1000,min=minWood,max=maxWood)
occurrence=species(CurrentTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
train1=data.frame(occurrence,temp=CurrentTemperature,woodlandCover)
plot(train1)
cor(train1)

FutureTemperature=CurrentTemperature+3
occurrence=species(FutureTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
future1=data.frame(occurrence,temp=FutureTemperature,woodlandCover)

# independent test data for current conditions
CurrentTemperature=runif(1000,min=minTemp,max=maxTemp)
woodlandCover=runif(1000,min=minWood,max=maxWood)
occurrence=species(CurrentTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
test1=data.frame(occurrence,temp=CurrentTemperature,woodlandCover)

#*************************************
# fit models

# generalised additive model
gam1=gam(occurrence~s(temp)+s(woodlandCover),binomial,data=train1)
# boosted regression tree
brt1 <- gbm.step(data=train1,gbm.x = c(2:3),gbm.y = 1,family = "bernoulli",
    tree.complexity = 1,learning.rate = 0.02,bag.fraction = 0.75)

#****************************************************
# predictions under current conditions
pred.gam=predict(gam1,newdata=test1,type="response")
pred.brt=predict.gbm(brt1,newdata=test1,n.trees=brt1$gbm.call$best.trees, type="response")

# plot fitted values
windows()
par(mfrow=c(2,2))
myplot(test1,pred.gam,main="GAM - current")
myplot(test1,pred.brt,main="BRT - current")

# plot response curves
windows()
par(mfrow=c(2,2))
inflated.response(gam1,train1[,2:3],main="GAM",method="stat6")
inflated.response(brt1,train1[,2:3],main="BRT",method="stat6")

#**********************************
# projections into future
pred.gam=predict(gam1,newdata=future1,type="response")
pred.brt=predict.gbm(brt1,newdata=future1,n.trees=brt1$gbm.call$best.trees, type="response")

# plot fitted values
windows()
par(mfrow=c(2,2))
myplot(future1,pred.gam,main="GAM - future")
myplot(future1,pred.brt,main="BRT - future")





#*****************************************************************************
#*****************************************************************************
#*****************************************************************************
# case 2: species niche truncated

# training data + future
CurrentTemperature=runif(1000,min=minTemp,max=maxTemp)
woodlandCover=numeric(1000)
# woodland occurs above 3°C and below 13°C with maximum woodland cover between 7 and 9°C
woodlandCover[CurrentTemperature>3&CurrentTemperature<=7]=
  sapply((CurrentTemperature[CurrentTemperature>3&CurrentTemperature<=7]-3)*maxWood/4,
  function(x){runif(1,min=minWood,max=x)})
woodlandCover[CurrentTemperature>7&CurrentTemperature<9]=
  runif(length(woodlandCover[CurrentTemperature>7&CurrentTemperature<9]),min=minWood,max=maxWood)
woodlandCover[CurrentTemperature>9&CurrentTemperature<13]=
  sapply((13-CurrentTemperature[CurrentTemperature>9&CurrentTemperature<13])*maxWood/4,
  function(x){runif(1,min=minWood,max=x)})
occurrence=species(CurrentTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
train2=data.frame(occurrence,temp=CurrentTemperature,woodlandCover)
plot(train2)
cor(train2)

FutureTemperature=CurrentTemperature+3
occurrence=species(FutureTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
future2=data.frame(occurrence,temp=FutureTemperature,woodlandCover)
plot(future2)
cor(future2)

# test data current conditions
CurrentTemperature=runif(1000,min=minTemp,max=maxTemp)
woodlandCover=numeric(1000)
woodlandCover[CurrentTemperature>3&CurrentTemperature<=7]=
  sapply((CurrentTemperature[CurrentTemperature>3&CurrentTemperature<=7]-3)*maxWood/4,
  function(x){runif(1,min=minWood,max=x)})
woodlandCover[CurrentTemperature>7&CurrentTemperature<9]=
  runif(length(woodlandCover[CurrentTemperature>7&CurrentTemperature<9]),min=minWood,max=maxWood)
woodlandCover[CurrentTemperature>9&CurrentTemperature<13]=
  sapply((13-CurrentTemperature[CurrentTemperature>9&CurrentTemperature<13])*maxWood/4,
  function(x){runif(1,min=minWood,max=x)})
occurrence=species(CurrentTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
test2=data.frame(occurrence,temp=CurrentTemperature,woodlandCover)


#***************************************
# fit models

# generalised additive model
gam2=gam(occurrence~s(temp)+s(woodlandCover),binomial,data=train2)
# boosted regression tree
brt2 <- gbm.step(data=train2,gbm.x = c(2:3),gbm.y = 1,family = "bernoulli",
    tree.complexity = 1,learning.rate = 0.01,bag.fraction = 0.75) 


#*********************************************
# predictions under current conditions
pred.gam=predict(gam2,newdata=test2,type="response")
pred.brt=predict.gbm(brt2,newdata=test2,n.trees=brt2$gbm.call$best.trees, type="response")

# plot fitted values
windows()
par(mfrow=c(2,2))
myplot(test2,pred.gam,main="GAM - current")
myplot(test2,pred.brt,main="BRT - current")

# plot response curves
windows()
par(mfrow=c(2,2))
inflated.response(gam2,train2[,2:3],main="GAM",method="stat6",disp='eo.mask')
inflated.response(brt2,train2[,2:3],main="BRT",method="stat6",disp='eo.mask')

#*************************************************
# projections into future
pred.gam=predict(gam2,newdata=future2,type="response")
pred.brt=predict.gbm(brt2,newdata=future2,n.trees=brt2$gbm.call$best.trees, type="response")

# plot fitted values
windows()
par(mfrow=c(2,2))
myplot(future2,pred.gam,main="GAM - future")
myplot(future2,pred.brt,main="BRT - future")




#****************************************************************************
#****************************************************************************
#****************************************************************************
# case 3: edge niche

# training data + future
CurrentTemperature=runif(1000,min=minTemp,max=maxTemp)
woodlandCover=numeric(1000)
woodlandCover[CurrentTemperature<10.5]=
  runif(length(woodlandCover[CurrentTemperature<10.5]),min=minWood,max=maxWood)
occurrence=species(CurrentTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
train3=data.frame(occurrence,temp=CurrentTemperature,woodlandCover)
plot(train3)
cor(train3)

FutureTemperature=CurrentTemperature+3
occurrence=species(FutureTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
future3=data.frame(occurrence,temp=FutureTemperature,woodlandCover)
plot(future3)
cor(future3)

# test data current conditions
CurrentTemperature=runif(1000,min=minTemp,max=maxTemp)
woodlandCover=numeric(1000)
woodlandCover[CurrentTemperature<10.5]=
  runif(length(woodlandCover[CurrentTemperature<10.5]),min=minWood,max=maxWood)
occurrence=species(CurrentTemperature,woodlandCover)
occurrence<-sapply(occurrence,function(x){rbinom(1,1,x)})
test3=data.frame(occurrence,temp=CurrentTemperature,woodlandCover)


#***************************************
# fit models

# generalised additive model
gam3=gam(occurrence~s(temp)+s(woodlandCover),binomial,data=train3)
# boosted regression tree
brt3 <- gbm.step(data=train3, gbm.x = c(2:3), gbm.y = 1, family = "bernoulli",
    tree.complexity = 1, learning.rate = 0.02, bag.fraction = 0.75)

#*********************************************
# predictions under current conditions
pred.gam=predict(gam3,newdata=test3,type="response")
pred.brt=predict.gbm(brt3,newdata=test3,n.trees=brt3$gbm.call$best.trees, type="response")

# plot fitted values
windows()
par(mfrow=c(2,2))
myplot(test3,pred.gam,main="GAM - current")
myplot(test3,pred.brt,main="BRT - current")

# plot response curves
windows()
par(mfrow=c(2,2))
inflated.response(gam3,train3[,2:3],main="GAM",method="stat6",disp='eo.mask')
inflated.response(brt3,train3[,2:3],main="BRT",method="stat6",disp='eo.mask')

#*************************************************
# projections into future
pred.gam=predict(gam3,newdata=future3,type="response")
pred.brt=predict.gbm(brt3,newdata=future3,n.trees=brt3$gbm.call$best.trees, type="response")

# plot fitted values
windows()
par(mfrow=c(2,2))
myplot(future3,pred.gam,main="GLM - future")
myplot(future3,pred.brt,main="BRT - future")
