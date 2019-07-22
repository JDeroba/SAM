rm(list=ls(all=T))
library(stockassessment)
source(paste("C:\\Users\\jonathan.deroba\\Documents\\GitHub\\SAM-ICES-WG-SR-one-off","functions.R",sep="\\")) #shouldn't have to touch
#library(TMB)

direct<-"B:\\jderoba\\ICES_AMWG_SR\\sa\\SAM_BH"

files<-list.files(direct) 



for(r in 1:1){ #loop over OM level folders
  r.folder<-paste(direct,files[r],sep="\\") #OM level (r) folder
  iters<-list.files(r.folder,pattern="iter") #datasets within OM folder
  converge.count<-0 #count the number of converged runs
  SRparms<-c() #matrix(NA,2,8)  #container for all Stock recruit param estimates
  for(i in 1:2){  #length(iters)){  #loop over datasets within OM
    
    iter<-paste(r.folder,iters[i],sep="\\") #single dataset folder
    #read in single dataset
    cn <- read.ices(paste(iter,"sim-LANUM.txt",sep="\\"))
    cw <- read.ices(paste(iter,"sim-WELAND.txt",sep="\\"))
    mo <- read.ices(paste(iter,"sim-MATPROP.txt",sep="\\"))
    nm <- read.ices(paste(iter,"sim-NATMOR.txt",sep="\\"))
    pf <- read.ices(paste(iter,"sim-FPROP.txt",sep="\\"))
    pm <- read.ices(paste(iter,"sim-MPROP.txt",sep="\\"))
    sw <- read.ices(paste(iter,"sim-WEST.txt",sep="\\"))
    surveys <- read.ices(paste(iter,"sim-TUNE.txt",sep="\\")) 
    
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
    conf$stockRecruitmentModelCode<-2 #estimate BH SR
    conf$fbarRange<-c(4,7)
    saveConf(conf,file=paste(iter,"ModelConf.txt",sep="\\"),overwrite=T)
    
    par<-defpar(dat,conf) #some default starting values
    #turn off survival process variance and set to 0
    par$logSdLogN[2]<-1 #set sd of logN age2-20 to 0 (1 in log space)
    fit<-try(sam.fit(dat,conf,par,run=T,map=list("logSdLogN"=factor(c(1,NA))))) #fit the model
    
    if(is.null(attr(fit,"condition"))){
    saveRDS(fit,file=paste(iter,"SAMfit.RData",sep="\\" )) #save the results
    converge.count<-converge.count+1
    modelTable<-modeltable(fit) #AIC and number of params
    write.csv(modelTable,file=paste(iter,"ModelTable.csv",sep="\\"))
    if(fit$conf$stockRecruitmentModelCode==2){
    alpha<-partable(fit)["rec_loga_0","exp(par)"]
    beta<-partable(fit)["rec_logb_0","exp(par)"]
    alpha.flr<-alpha/beta
    beta.flr<-1/beta
    
    natMor<- fit$data$natMor[nrow(fit$data$natMor),]
    StWt<-fit$data$stockMeanWeight[nrow(fit$data$stockMeanWeight),]
    Mature<-fit$data$propMat[nrow(fit$data$propMat),]
    SRparmsa<-SR.parms(nage=max(fit$data$maxAgePerFleet),M=natMor,Wt=StWt,Mat=Mature,alpha=alpha,beta=beta,type=fit$conf$stockRecruitmentModelCode)
    SRparms<-rbind(SRparms,cbind(SRparmsa,alpha.flr,beta.flr,files[r],iters[i]) )
    } 
    
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
    } #close if/try
  } #i loop
  write.csv(SRparms,paste(direct,paste(files[r],"SRparms.csv",sep="\\"),sep="\\"))
  write.csv(converge.count,paste(direct,paste(files[r],"converge_count.csv",sep="\\"),sep="\\"))
} # r loop
