# 'Predicting to new environments: tools for visualising model behaviour and impacts on mapped distributions'
# by Damaris Zurell, Jane Elith and Boris Schroeder.
# Appendix S1 - accompanying functions
# last update 2011-09-29

#-------------------------------------------------------------------------------
# define species (lrm)
species=function(temp,wood,sdev=3){

  # creates artifical species with two-dimensional niche;
  # species shows positive linear response to woodland cover and unimodal response to temperature
  #
  # written by Damaris Zurell, last update 2011-02-17

  return(inv.logit(-170+40*temp-2.5*temp^2+
  .35*wood+
  rnorm(max(length(temp),length(wood)),0,sdev)))}

#-------------------------------------------------------------------------------
# plot fitted values for all three predictors + true response curve
myplot<- function(x,y,main,thresh=F,ylab="Occurrence probability"){

  # plot fitted values along with true species response
  #
  # written by Damaris Zurell, last update 2011-02-17
  
  plot(x$temp,y,xlab="Temperature",ylab=ylab,main=main,ylim=c(0,1),xlim=c(3,16))
  i=seq(3,16,length=100)
  lines(i,species(i,70,sdev=0),lty="solid",col="grey80",lwd=2)
  plot(x$woodlandCover,y,xlab="Woodland cover",ylab="",main=main,ylim=c(0,1))
  i=seq(min(x$woodlandCover),max(x$woodlandCover),length=100)
  lines(i,species(8,i,sdev=0),lty="solid",col="grey80",lwd=2)
  }
  

#-------------------------------------------------------------------------------
# plot inflated response curves - inflated partial dependence plots
inflated.response=function(object,predictors,select.columns=NULL,label=NULL,len=50,lhsample=100,lwd=1,
    ylab="Occurrence probabilities",method="stat3",disp="all",overlay.mean=T,
    col.curves='grey',col.novel='grey',col.mean='black',lwd.known=2,lwd.mean=2,...){
    
  # plot inflated response curves;
  # plot effect of one variable over range of other predictors;
  # method determines at which values the other predictors are held constant:
  # method='mean' corresponds to conventional partial dependence plots,
  # method='stat3' (Default) considers minimum, mean and maximum values of predictors,
  # method='stat6' considers min,mean,median,max and quartiles.
  # for 'stat3' and 'stat6' effects of one variables is plotted for all possible
  # combinations of remaining predictors - as the number of combinations increases exponentially,
  # the maximum number of combinations can be set with lhsample. Whenever lhsample is exceeded,
  # candidate combinations are drawn by latin hypercube sampling.
  # len determines the number of intervals along the environmental gradient plotted,
  # i.e. smoothness of response curves.
  # disp can take options 'all' or 'eo.mask' - in the latter case, eo.mask() is used
  # to distinguish between areas of the estimated environmental niche / plotting areas
  # that are supported by data and those that require extrapolation.
  # if overlay.mean is true, then the mean response curve is overlaid on the inflated plot.
  #
  # written by Damaris Zurell, last update 2011-09-17
  
  if (is.null(select.columns)) select.columns=1:ncol(predictors)
  
  require(lhs,quietly=T)
  for (i in select.columns)
  {
  summaries=data.frame(matrix(0,6,ncol(predictors)))
  for (iz in 1:ncol(predictors)) summaries[,iz]=summary(predictors[,iz])
  if (method=="stat3") {summaries.j=as.matrix(summaries[c(1,4,6),-i],ncol=(ncol(predictors)-1));comb=min(lhsample,3^(ncol(predictors)-1));nc=3} else
  if (method=="stat6") {summaries.j=as.matrix(summaries[,-i],ncol=(ncol(predictors)-1));comb=min(lhsample,6^(ncol(predictors)-1));nc=6} else
  if (method=="mean") {summaries.j=as.matrix(summaries[4,-i],ncol=(ncol(predictors)-1));comb=1;nc=1;overlay.mean=F}
  dummy.j=as.matrix(predictors[1:len,-i],ncol=(ncol(predictors)-1))
  if (comb<lhsample) {
    mat=vector("list",ncol(dummy.j))
    for (m in 1:ncol(dummy.j)) mat[[m]]=1:nc
    mat=expand.grid(mat)
    } else
  mat=round(qunif(randomLHS(lhsample,ncol(dummy.j)),1,nrow(summaries.j)),0)
  if (is.null(label)) label=names(predictors)

  for (r in 1:nrow(mat))
    {
      for (j in 1:ncol(dummy.j))
      {
      dummy.j[,j]=as.vector(rep(summaries.j[mat[r,j],j],len))
      }

    dummy=data.frame(seq(min(predictors[,i]),max(predictors[,i]),length=len),dummy.j)
    names(dummy)[-1]=names(predictors)[-i]
    names(dummy)[1]=names(predictors)[i]

    if (is(object,"gbm")) curves<-predict.gbm(object, dummy,n.trees=object$gbm.call$best.trees, type="response") # when using brt code from Elith et al. (2008) JAnimEcol
    else if (is(object,"glm")) curves<-predict(object, dummy, type="response")
    else if (is(object,"randomForest")) curves<-predict(object,dummy)
    else if (is(object,"tree")) curves<-predict(object,dummy)
    else if (is(object,"list")) curves<-mars.predict(object, dummy)$prediction[[1]]   #when using mars code from Elith & Leathwick (2007) Div Distr
    else if (is(object,"fda")) curves<-predict(object,dummy,type="post")[,2]
    else if (is(object,"nnet")) curves<-predict(object,dummy,type="raw")
    else {print("SDM class unknown");break}
    
    # display all lines in same type
    if (disp=='all')
    {
    if (r==1)
    {
    if (i==1) plot(dummy[,names(predictors)[i]],
      curves,type="l",ylim=c(0,1),xlab=label[i],ylab=ylab,
      lwd=lwd,col=col.curves,...)
    else plot(dummy[,names(predictors)[i]],
      curves,type="l",ylim=c(0,1),xlab=label[i],ylab="",lwd=lwd,col=col.curves,...)
    }
    else lines(dummy[,names(predictors)[i]],
      curves,lwd=lwd,col=col.curves,...)
    }
    
    # highlight extrapolation to novel environmental conditions
    if (disp=='eo.mask')
    {
    novel=eo.mask(predictors,dummy)
    curves.known=curves
    curves.known[novel==1]=NA
    curves.novel=curves
    curves.novel[novel==0]=NA
    
    if (r==1)
    {
    if (i==1) {plot(dummy[,names(predictors)[i]],
      curves.known,type="l",ylim=c(0,1),xlab=label[i],ylab=ylab,
      lwd=lwd.known,col=col.curves,...)
      lines(dummy[,names(predictors)[i]],
      curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
    else {plot(dummy[,names(predictors)[i]],
      curves.known,type="l",ylim=c(0,1),xlab=label[i],ylab="",lwd=lwd.known,col=col.curves,...)
      lines(dummy[,names(predictors)[i]],
      curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
    }
    else {lines(dummy[,names(predictors)[i]],
      curves.known,lwd=lwd.known,col=col.curves,...)
      lines(dummy[,names(predictors)[i]],
      curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
    }
    }
    
    #-------------------------------------------------
    # now, this is for overlaying mean response curve
    if (overlay.mean==T)
    {
    dummy=predictors[1:len,]
    dummy[,i]=seq(min(predictors[,i]),max(predictors[,i]),length=len)
    for (j in 1:ncol(predictors))
      {
      if (j!=i) 
        {
        dummy[,j]=rep(mean(predictors[,j]),len)
        }
      }
    
    if (is(object,"gbm")) curves<-predict.gbm(object, dummy,n.trees=object$gbm.call$best.trees, type="response")
    else if (is(object,"glm")) curves<-predict(object, dummy, type="response")
    else if (is(object,"randomForest")) curves<-predict(object,dummy)
    else if (is(object,"tree")) curves<-predict(object,dummy)
    else if (is(object,"list")) curves<-mars.predict(object, dummy)$prediction[[1]]
    else if (is(object,"fda")) curves<-predict(object,dummy,type="post")[,2]
    else if (is(object,"nnet")) curves<-predict(object,dummy,type="raw")
    else {print("SDM class unknown");break}

    lines(dummy[,names(predictors)[i]],
      curves,lwd=lwd.mean,col=col.mean,...)
    }    
  }}


#-------------------------------------------------------------------------------
# calculate environmental overlap mask
# extension of MESS that was proposed by Elith et al. 2010 MethodsEcolEvol 1:330-342.

eo.mask=function(traindata,newdata,nbin=5,type="EO")
  {
  # a bin size of one corresponds to MESS
  # type 'EO' returns a vector of zeros and ones for analog(0) and novel(1) environments
  # type 'ID' returns a character vector defining the combination of bins each data entry 
  # belongs to - this may help finding the problem maker parts of the prediction space
  
  train.minima=apply(traindata,2,min)
  train.maxima=apply(traindata,2,max)
  
  train.ids=apply(apply(ceiling(apply(round(
    sweep(sweep(traindata, 2, train.minima, "-"), 2, train.maxima - train.minima, "/")*nbin,4),
    c(1,2),FUN=function(x){if(x==0)x=1 else x=x})),
    c(1,2),FUN=function(x){if(x<1)x=0 else if(x>nbin)x=nbin+1 else x=x}),1,paste,collapse=".")
  
  new.ids=apply(apply(ceiling(apply(round(
    sweep(sweep(newdata[,names(train.minima)], 2, train.minima, "-"), 2, train.maxima - train.minima, "/")*nbin,4),
    c(1,2),FUN=function(x){if(x==0)x=1 else x=x})),
    c(1,2),FUN=function(x){if(x<1)x=0 else if(x>nbin)x=nbin+1 else x=x}),1,paste,collapse=".")
    
  if (type=="ID") return(new.ids)
  else if (type=="EO") return(sapply(new.ids%in%train.ids,FUN=function(x){if(x==T) x=0 else if(x==F)x=1}))    
  }  
  
