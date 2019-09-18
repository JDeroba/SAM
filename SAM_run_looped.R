rm(list=ls(all=T))
library(stockassessment)
source(paste("C:\\Users\\jonathan.deroba\\Documents\\GitHub\\SAM-ICES-WG-SR-one-off","functions.R",sep="\\")) #shouldn't have to touch
#library(TMB)

direct<-"B:\\jderoba\\ICES_AMWG_SR\\sa\\SAM_BH"

files<-list.files(direct) 

for(r in 4:4){ #loop over OM level folders
  r.folder<-paste(direct,files[r],sep="\\") #OM level (r) folder
  iters<-list.files(r.folder,pattern="iter") #datasets within OM folder
  converge.count<-0 #count the number of converged runs
  SRparms<-c() #matrix(NA,2,8)  #container for all Stock recruit param estimates
  for(i in 1:1){  #loop over datasets within OM    length(iters)
    write(i,paste(direct,paste(files[r],"loopcount.txt",sep="\\"),sep="\\"))
    iter<-paste(r.folder,iters[i],sep="\\") #single dataset folder
    #fit<-readRDS(file=paste(r.folder,iters[2],"SAMfit.RData",sep="\\")) #read old result back-in; I noticed some plots aren't made correctly when you read in old results.
    #read in single dataset
    cn <- read.ices(paste(iter,"sim-LANUM.txt",sep="\\"))
    cw <- read.ices(paste(iter,"sim-WELAND.txt",sep="\\"))
    mo <- read.ices(paste(iter,"sim-MATPROP.txt",sep="\\"))
    nm <- read.ices(paste(iter,"sim-NATMOR.txt",sep="\\"))
    pf <- read.ices(paste(iter,"sim-FPROP.txt",sep="\\"))
    pm <- read.ices(paste(iter,"sim-MPROP.txt",sep="\\"))
    sw <- read.ices(paste(iter,"sim-WEST.txt",sep="\\"))
    surveys <- read.ices(paste(iter,"sim-TUNE.txt",sep="\\")) 
    
    cn.p<-rowSums(cn[,8:ncol(cn)])
    cn<-cbind(cn[,1:7],cn.p)
    colnames(cn)[8]<-8
    cw.p<-rowMeans(cw[,8:ncol(cw)])
    cw<-cbind(cw[,1:7],cw.p)
    colnames(cw)[8]<-8
    mo.p<-rowMeans(mo[,8:ncol(mo)])
    mo<-cbind(mo[,1:7],mo.p)
    colnames(mo)[8]<-8
    nm.p<-rowMeans(nm[,8:ncol(nm)])
    nm<-cbind(nm[,1:7],nm.p)
    colnames(nm)[8]<-8
    pf.p<-rowMeans(pf[,8:ncol(pf)])
    pf<-cbind(pf[,1:7],pf.p)
    colnames(pf)[8]<-8
    pm.p<-rowMeans(pm[,8:ncol(pm)])
    pm<-cbind(pm[,1:7],pm.p)
    colnames(pm)[8]<-8
    sw.p<-rowMeans(sw[,8:ncol(sw)])
    sw<-cbind(sw[,1:7],sw.p)
    colnames(sw)[8]<-8
    temp.p<-rowSums(surveys[[1]][,8:ncol(surveys[[1]])])
    temp<-cbind(surveys[[1]][,1:7],temp.p)
    colnames(temp)[8]<-8
    attr(temp,"time")<-c(0,0)
    #attributes(temp)
    surveys<-temp
    
    ages<-seq(1,ncol(cn),1) #just ages
    
    #create SAM dat
    dat <- setup.sam.data(surveys=surveys,
                          residual.fleet=cn, 
                          prop.mature=mo, 
                          stock.mean.weight=sw, 
                          prop.f=pf, 
                          prop.m=pm, 
                          natural.mortality=nm,
                          #land.frac=lf,
                          #dis.mean.weight = dw,
                          #land.mean.weight = lw,
                          catch.mean.weight=cw)
    
    
    conf<-defcon(dat) #a default configuration for SAM; 
    conf$maxAgePlusGroup<-1
    conf$stockRecruitmentModelCode<-2 #estimate BH SR
    conf$fbarRange<-c(4,7)
    conf$keyLogFsta[1,]<-c(0,1,2,3,3,3,3,3)
    saveConf(conf,file=paste(iter,"ModelConf.txt",sep="\\"),overwrite=T)
    
    par<-defpar(dat,conf) #some default starting values
    #turn off survival process variance and set to 0
    par$logSdLogN[2]<-1 #set sd of logN age2-20 to 0 (1 in log space)
    #par$logSdLogObs[1]<--7 #set catch sd to some value (e.g., zero, truth)
    #par$logSdLogObs[2]<-1
    par$itrans_rho<-exp(1)
    fit<-try(sam.fit(dat,conf,par,run=T,map=list("itrans_rho"=factor(c(NA)),"logSdLogN"=factor(c(1,NA))))) #fit the model
    
    if(is.null(attr(fit,"condition"))){
    saveRDS(fit,file=paste(iter,"SAMfit.RData",sep="\\" )) #save the results
    modelTable<-modeltable(fit) #AIC and number of params
    write.csv(modelTable,file=paste(iter,"ModelTable.csv",sep="\\"))
    
    if(fit[[6]][3]==0){  #check for model convergence and record in SRparms below
    converge.count<-converge.count+1
    converge<-"yes"
    } else {
    converge<-"no"
    }
    
    if(fit$conf$stockRecruitmentModelCode==2){
    alpha<-partable(fit)["rec_loga_0","exp(par)"]
    beta<-partable(fit)["rec_logb_0","exp(par)"]
    alpha.flr<-alpha/beta
    beta.flr<-1/beta
    
    sigmaR<-partable(fit)["logSdLogN_0","exp(par)"]
    sigmacatch<-partable(fit)["logSdLogObs_0","exp(par)"]
    natMor<- fit$data$natMor[nrow(fit$data$natMor),]
    StWt<-fit$data$stockMeanWeight[nrow(fit$data$stockMeanWeight),]
    Mature<-fit$data$propMat[nrow(fit$data$propMat),]
    SRparmsa<-SR.parms(nage=max(fit$data$maxAgePerFleet),M=natMor,Wt=StWt,Mat=Mature,alpha=alpha,beta=beta,type=fit$conf$stockRecruitmentModelCode)
    SRparms<-rbind(SRparms,cbind(SRparmsa,alpha.flr,beta.flr,files[r],iters[i],"converge"=converge,sigmaR,sigmacatch) ) #
    } else {
      if(fit$conf$stockRecruitmentModelCode==0){
        sigmaR<-partable(fit)["logSdLogN_0","exp(par)"]
        sigmacatch<-partable(fit)["logSdLogObs_0","exp(par)"]
        SRparms<-rbind(SRparms,cbind(files[r],iters[i],"converge"=converge,sigmaR,sigmacatch) ) #
      }
    } #if Bev Holt
    
    ##Will save fit details and diagnostic plots pdf if TRUE
    if(FALSE){
      qplotnames<-getnames(confa=conf$keyLogFpar[2:nrow(conf$keyLogFpar),],agesa=ages,looptoa=(max(conf$keyLogFpar)+1),looptob=(length(names(surveys))),namesa=names(surveys))
      Fpronames<-getnames(confa=t(as.matrix(conf$keyVarF[1,])),agesa=ages,looptoa=(max(conf$keyVarF)+1),looptob=1,namesa="SD LogF process")
      Npronames<-getnames(confa=t(as.matrix(conf$keyVarLogN)),agesa=ages,looptoa=(max(conf$keyVarLogN)+1),looptob=1,namesa="SD LogN process")
      Obsnames<-getnames(confa=conf$keyVarObs,agesa=ages,looptoa=(max(conf$keyVarObs)+1),looptob=nrow(conf$keyVarObs),namesa=c("Catch",names(surveys)))
      Obsnames<-paste("SD",Obsnames)
      varplotnames<-c(Fpronames,Npronames[1],Obsnames) #only do first N pro because other set to zero
      Fstanames<-getnames(confa=t(as.matrix(conf$keyLogFsta[1,])),agesa=ages,looptoa=(max(conf$keyLogFsta)+1),looptob=1,namesa="F")
      #make the plots and save to pdf in Run sub-directory
      plotfxn(afit=fit,datdirect=r.folder,run=iters[i],confa=conf,qplotnamesa=qplotnames,varplotnamesa=varplotnames,fstanamesa=Fstanames,catchmultnamesa=catchmultnames,retroyrs=2) #as.vector(seq(2011,2014)))
    } #end if for diagnostic plots
    #code to extract and record var covar courtesy of anders 9/18/19
    sdr <- TMB:::sdreport(fit$obj)
    nam <- names(sdr$value)
    idx <- nam%in%c("logR","logssb")
    cov <- sdr$cov[idx,idx]
    rownames(cov) <- nam[idx]
    colnames(cov) <- nam[idx]
    write.csv(cov,paste(direct,paste(files[r],iters[i],"SRCovar.csv",sep="\\"),sep="\\"))
    
    } #close if/try
  } #i loop
  write.csv(SRparms,paste(direct,paste(files[r],"SRparms.csv",sep="\\"),sep="\\"))
  write.csv((converge.count/length(iters)),paste(direct,paste(files[r],"converge_count.csv",sep="\\"),sep="\\"))
  
  quickplot<-function(true=NA,xwant=NA,xlab=NA,dat=NA){
    plot(y=seq(1:nrow(dat)),x=dat[,xwant],ylab="Iter",xlab=xlab)
    abline(v=true,lty=2,lwd=2)
    abline(v=median(dat[,xwant]),lty=3,col="red",lwd=2)
    legend("topright",legend=c("True","Median"),text.col=c("black","red"))
  }
  
  forpdf<-paste(files[r],".pdf",sep="")
  pdf(paste(direct,paste(files[r],forpdf,sep="\\"),sep="\\")) #create pdf for graph storage
  
  quickplot(true=0.65,xwant="steep",xlab="Steepness - all runs",dat=SRparms)
  quickplot(true=0.65,xwant="steep",xlab="Steepness - converged",dat=SRparms[SRparms$converge %in% c("yes"),])
  if(nrow(SRparms[SRparms$converge %in% c("no"),])>0) {
  quickplot(true=0.65,xwant="steep",xlab="Steepness - not converged",dat=SRparms[SRparms$converge %in% c("no"),])
  }
  
  quickplot(true=1000,xwant="B0",xlab="Unfished Biomass - all runs",dat=SRparms)
  quickplot(true=1000,xwant="B0",xlab="Unfished Biomass - converged",dat=SRparms[SRparms$converge %in% c("yes"),])
  if(nrow(SRparms[SRparms$converge %in% c("no"),])>0){
  quickplot(true=1000,xwant="B0",xlab="Unfished Biomass -  not converged",dat=SRparms[SRparms$converge %in% c("no"),])
  }
  
  quickplot(true=0.3,xwant="sigmaR",xlab="SigmaR - all runs",dat=SRparms)
  quickplot(true=0.3,xwant="sigmaR",xlab="SigmaR - converged",dat=SRparms[SRparms$converge %in% c("yes"),])
  if(nrow(SRparms[SRparms$converge %in% c("no"),])>0) {
  quickplot(true=0.3,xwant="sigmaR",xlab="SigmaR - not converged",dat=SRparms[SRparms$converge %in% c("no"),])
  }
  
  quickplot(true=0.0,xwant="sigmacatch",xlab="Catch Obs Std. - all runs",dat=SRparms)
  quickplot(true=0.0,xwant="sigmacatch",xlab="Catch Obs Std. - converged",dat=SRparms[SRparms$converge %in% c("yes"),])
  if(nrow(SRparms[SRparms$converge %in% c("no"),])>0) {
  quickplot(true=0.0,xwant="sigmacatch",xlab="Catch Obs Std. - not converged",dat=SRparms[SRparms$converge %in% c("no"),])
  }
  dev.off() #close pdf
} # r loop








