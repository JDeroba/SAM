rm(list=ls(all=T))
library(stockassessment)
direct.liz<-"H:\\ICES_AMWG_SR\\sa\\ASAP_Liz"
directs<-c("H:\\ICES_AMWG_SR\\sa\\SAM_BH","H:\\ICES_AMWG_SR\\sa\\SAM_RICKER","H:\\ICES_AMWG_SR\\sa\\SAM_RW")
source(paste("C:\\Users\\jonathan.deroba\\Documents\\GitHub\\SAM-ICES-WG-SR-one-off","functions.R",sep="\\")) #shouldn't have to touch

#load(file=paste(direct.liz,"asap_output_meanSR_fit.RData",sep="\\"))

for(d in 1:2) { #length(directs)){ #loop over BH and Ricker and RW fits
  direct<-directs[d]
  files<-list.files(direct) 
  
  for(r in 1:2) { #length(files)){ #loop over OM level folders
    r.folder<-paste(direct,files[r],sep="\\") #OM level (r) folder
    iters<-list.files(r.folder,pattern="iter") #datasets within OM folder
    run<-files[r]
    for(i in 1:5) { #length(iters)){  #loop over datasets within OM    
      iter<-paste(r.folder,iters[i],sep="\\") #single dataset folder
      fit<-readRDS(file=paste(r.folder,iters[i],"SAMfit.RData",sep="\\")) #read old result back-in;
      
      iter.b<-iters[i]
      AIC<-modeltable(fit)[3]
      nll<--1*modeltable(fit)[1]
      npar<-modeltable(fit)[2]
      max.grad<-max(fit$sdrep$gradient.fixed)
      
      sigmaR<-partable(fit)["logSdLogN_0","exp(par)"]
      
      natMor<- fit$data$natMor[nrow(fit$data$natMor),]
      StWt<-fit$data$stockMeanWeight[nrow(fit$data$stockMeanWeight),]
      Mature<-fit$data$propMat[nrow(fit$data$propMat),]
      select<-matrix(nrow=nrow(faytable(fit)),ncol=ncol(faytable(fit)))
      for(s in 1:nrow(faytable(fit))){
        select[s,]<-faytable(fit)[s,]/max(faytable(fit)[s,])
      }
      
      if(fit[[6]][3]==0){  #check for model convergence and record in SRparms below
        conv<-1
      } else {
        conv<-NA
      }
      
      if(fit$conf$stockRecruitmentModelCode==2){  #BH
        sr.code<-1
        alpha<-partable(fit)["rec_loga_0","exp(par)"]
        beta<-partable(fit)["rec_logb_0","exp(par)"]
        alpha.flr<-alpha/beta
        beta.flr<-1/beta
        SRparmsa<-SR.parms(nage=max(fit$data$maxAgePerFleet),M=natMor,Wt=StWt,Mat=Mature,alpha=alpha,beta=beta,type=fit$conf$stockRecruitmentModelCode)
        MSYparms<-MSY.find(M=natMor,Wt=StWt,Mat=Mature,steep=SRparmsa$steep,selectivity.F = select[nrow(select),],B0=SRparmsa$B0,type=fit$conf$stockRecruitmentModelCode,med.recr=0,SSBR0 = SRparmsa$SSBR0,B.MSYflag = F)
        Bmsy<-max.ypr(F=MSYparms$F.MSY ,M=natMor,Wt=StWt,Mat=Mature,steep=SRparmsa$steep,selectivity.F = select[nrow(select),],B0=SRparmsa$B0,type=fit$conf$stockRecruitmentModelCode,med.recr=0,SSBR0 = SRparmsa$SSBR0,B.MSYflag = T)
      } else if(fit$conf$stockRecruitmentModelCode==1) {  #Ricker 
        sr.code<-2
        alpha<-partable(fit)["rec_loga_0","exp(par)"]
        beta<-partable(fit)["rec_logb_0","exp(par)"] 
        SRparmsa<-SR.parms(nage=max(fit$data$maxAgePerFleet),M=natMor,Wt=StWt,Mat=Mature,alpha=alpha,beta=beta,type=fit$conf$stockRecruitmentModelCode)
        MSYparms<-MSY.find(M=natMor,Wt=StWt,Mat=Mature,steep=SRparmsa$steep,selectivity.F = select[nrow(select),],B0=SRparmsa$B0,type=fit$conf$stockRecruitmentModelCode,med.recr=0,SSBR0 = SRparmsa$SSBR0,B.MSYflag = F)
        Bmsy<-max.ypr(F=MSYparms$F.MSY ,M=natMor,Wt=StWt,Mat=Mature,steep=SRparmsa$steep,selectivity.F = select[nrow(select),],B0=SRparmsa$B0,type=fit$conf$stockRecruitmentModelCode,med.recr=0,SSBR0 = SRparmsa$SSBR0,B.MSYflag = T)
      } else if(fit$conf$stockRecruitmentModelCode==0) {  #RW recruitment
        sr.code<-3
        idx <- names(fit$sdrep$value) == "logR"
        recruitment <- exp(fit$sdrep$value[idx])
        alpha<-mean(recruitment) #not really alpha, just mean recruitment; makes SR.parms function usable
        beta<-NA
        SRparmsa<-SR.parms(nage=max(fit$data$maxAgePerFleet),M=natMor,Wt=StWt,Mat=Mature,alpha=alpha,beta=NA,type=fit$conf$stockRecruitmentModelCode)
        SRparmsa$steep<-NA
        MSYparms<-MSY.find(M=natMor,Wt=StWt,Mat=Mature,steep=NULL,selectivity.F = select[nrow(select),],B0=SRparmsa$B0,type=fit$conf$stockRecruitmentModelCode,med.recr=alpha,SSBR0 = SRparmsa$SSBR0,B.MSYflag = F)
        Bmsy<-max.ypr(F=MSYparms$F.MSY ,M=natMor,Wt=StWt,Mat=Mature,steep=NULL,selectivity.F = select[nrow(select),],B0=SRparmsa$B0,type=fit$conf$stockRecruitmentModelCode,med.recr=alpha,SSBR0 = SRparmsa$SSBR0,B.MSYflag = T)
      }
      
      idxr <- names(fit$sdrep$value) == "logR"
      recruitment <- exp(fit$sdrep$value[idxr])
      idxb <- names(fit$sdrep$value) == "logssb"
      ssb <- exp(fit$sdrep$value[idxb])
      idxf <- names(fit$sdrep$value) == "logfbar"
      F.ave <- exp(fit$sdrep$value[idxf])
      years<-fit$data$years
      srest.a<-cbind(rep(sr.code,length(years)),rep(run,length(years)),rep(iter.b,length(years)),years,ssb,recruitment,F.ave)
      colnames(srest.a)<-c("sr.code","run","iter","year","ssb","recr","F.ave")
      rownames(srest.a)<-NULL
      
      vcv.a<-read.csv(file=paste(r.folder,iters[i],"SRCovar.csv",sep="\\"))
      vcv.a<-vcv.a[,!names(vcv.a) %in% c("X")]
      vcv.a$run<-rep(run,length(years)*2)
      vcv.a$iter<-rep(iter.b,length(years)*2)
      vcv.a$sr.code<-rep(sr.code,length(years)*2)
      
      sr.a<-cbind(run,iter.b,alpha,beta,SRparmsa$B0,SRparmsa$SSBR0,SRparmsa$steep,MSYparms$F.MSY,Bmsy,MSYparms$MSY,sr.code,nll,npar,max.grad,conv,AIC,sigmaR)
      colnames(sr.a)<-c("run","iter","alpha","beta","SSB0","ssb.per.rec0","steepness","Fmsy","SSBmsy","MSY","sr.code","nll","npar","max.grad","conv","AIC","sigmaR")
      if(i==1 & r==1) {
        sr<-sr.a
        srest<-srest.a
        vcv<-vcv.a
      } else {
        sr<-rbind(sr,sr.a)
        srest<-rbind(srest,srest.a)
        vcv<-rbind(vcv,vcv.a)
      }
    } #close i
  } #close r
  ReadMe<-"vcv is the variance covariance matrix for logR and logSSB.  For each run and iter, there is
   a matrix 2*nyears X 2*nyears with logssb[year1:nyears]...logR[year1:nyears].  sr contains the sr parameter
   estimates, MSY values, and other self explanatory output from each dataset.  srest contains the
   SSB, recruitment, and F.ave (ages 4-7) estimates for each dataset.  sr.code columns are similar to Liz
   1=BH, 2=Ricker, but 3=random walk for these SAM fits."
  
  filename<-paste0("SAM_fits_srcode",d,".RData")
  save(ReadMe,sr,srest,vcv,file=paste(directs[d],filename,sep="\\"))
  
} #close d
#load(file=paste(directs[2],"SAM_fits_srcode2.RData",sep="\\"))
#head(vcv)
