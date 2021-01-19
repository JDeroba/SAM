rm(list=ls(all=T))
#devtools::install_github("fishfollower/SAM/stockassessment")
library(stockassessment)
directs<-c("H:\\ICES_AMWG_SR\\sa\\SAM_BH","H:\\ICES_AMWG_SR\\sa\\SAM_RICKER","H:\\ICES_AMWG_SR\\sa\\SAM_RW")

#load(file=paste(directs[3],"SAM_fits_srcode3.RData",sep="\\"))

#AIC.BH<-data.frame(run=sr[,"run"],iter=sr[,"iter"],BH.AIC=sr[,"AIC"],BH.nll=sr[,"nll"],BH.max.grad=sr[,"max.grad"],BH.conv=sr[,"conv"],BH.npar=sr[,"npar"])
#AIC.Ricker<-data.frame(run=sr[,"run"],iter=sr[,"iter"],Ricker.AIC=sr[,"AIC"],Ricker.nll=sr[,"nll"],Ricker.max.grad=sr[,"max.grad"],Ricker.conv=sr[,"conv"],Ricker.npar=sr[,"npar"])
#AIC.RW<-data.frame(run=sr[,"run"],iter=sr[,"iter"],RW.AIC=sr[,"AIC"],RW.nll=sr[,"nll"],RW.max.grad=sr[,"max.grad"],RW.conv=sr[,"conv"],RW.npar=sr[,"npar"])

#rm(sr)
#rm(srest)
#rm(vcv)

AIC.ALL<-merge(AIC.BH,AIC.Ricker,by=c("run","iter"),all=TRUE)
AIC.ALL<-merge(AIC.ALL,AIC.RW,by=c("run","iter"),all=TRUE)
#head(AIC.ALL)

AIC.Winner<-lapply(X=1:nrow(AIC.ALL),FUN=function(x){
  temp<-data.frame(AIC.ALL[x,"BH.AIC"],AIC.ALL[x,"Ricker.AIC"],AIC.ALL[x,"RW.AIC"])
  if(which.min(temp)==1){
    winner="BH"}
  if(which.min(temp)==2){
    winner="Ricker"}
  if(which.min(temp)==3){
    winner="RW"}
  return(winner)
  }
  )
AIC.Winner= do.call(rbind.data.frame, AIC.Winner)
names(AIC.Winner)=c("Winner")
AIC.ALL$Winner=AIC.Winner
#save(AIC.ALL,file=paste("H:\\ICES_AMWG_SR\\sa","SAM_AIC_TABLE.RData",sep="\\"))
table(AIC.ALL$Winner)

########################################################
load(file=paste("H:\\ICES_AMWG_SR\\sa","SAM_AIC_TABLE.RData",sep="\\"))
AIC.ALL$Winner=as.vector(AIC.ALL$Winner$Winner) #quick data munge
for.truth=read.csv(file=paste("H:\\ICES_AMWG_SR\\sa","runs_liz.csv",sep="\\"))
for.truth$run=paste0("r",for.truth$case)
for.truth$model=ifelse(for.truth$model=="bevholt","BH",for.truth$model)
for.truth$model=ifelse(for.truth$model=="ricker","Ricker",for.truth$model)

for.truth$model_jjd=ifelse(for.truth$model=="BH","BH",for.truth$model)
for.truth$model_jjd=ifelse(for.truth$model=="Ricker","Ricker",for.truth$model_jjd)
for.truth$model_jjd=ifelse(for.truth$model=="segreg","BH",for.truth$model_jjd)
for.truth$model_jjd=ifelse(for.truth$model=="mean","RW",for.truth$model_jjd)

AIC.ALL=merge(AIC.ALL,for.truth,by="run")
AIC.ALL$correct_jjd=ifelse(AIC.ALL$model_jjd==AIC.ALL$Winner,1,0)
AIC.ALL$correct_om=ifelse(AIC.ALL$model==AIC.ALL$Winner,1,0)

#table(AIC.ALL$model_jjd,AIC.ALL$correct_jjd)
#table(AIC.ALL$model,AIC.ALL$correct_om)

library(janitor)
model.winner=tabyl(AIC.ALL,model,Winner)
model.winner$model=paste0(model.winner$model,".om")
model.winner.perc=round(model.winner[,2:4]/rowSums(model.winner[,2:4]),3)
model.winner.perc$model=model.winner$model

model.winner.jjd=tabyl(AIC.ALL,model_jjd,Winner)
model.winner.jjd$model=paste0(model.winner.jjd$model,".om")
model.winner.perc.jjd=round(model.winner.jjd[,2:4]/rowSums(model.winner.jjd[,2:4]),3)
model.winner.perc.jjd$model=model.winner.jjd$model

#function to calcuate proportion winner for three way contingency
threeway=function(dat=NULL){
dat.prop=lapply(1:length(dat),function(x) {
  tempsums=rowSums(dat[[x]][,2:4])
  props=dat[[x]][,2:4]/tempsums
  rownames(props)=paste0(dat[[x]][,1],".om")
  return(props)
})
names(dat.prop)=names(dat)
return(dat.prop)
}
#by sigmaR
model.winner.sigR=tabyl(AIC.ALL,model,Winner,sigmaR)
model.winner.sigR.prop=threeway(model.winner.sigR)
#by devs
model.winner.devs=tabyl(AIC.ALL,model,Winner,devs)
model.winner.devs.prop=threeway(model.winner.devs)

playchi=chisq.test(model.winner.devs[[3]])
