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
for.truth=read.csv(file=paste("H:\\ICES_AMWG_SR\\sa","runs_liz.csv",sep="\\"))
#truth=data.frame("case"=for.truth$case,"model"=for.truth$model)
for.truth$run=paste0("r",for.truth$case)
for.truth$model=ifelse(for.truth$model=="bevholt","BH",for.truth$model)
for.truth$model=ifelse(for.truth$model=="ricker","Ricker",for.truth$model)

for.truth$model_jjd=ifelse(for.truth$model=="BH","BH",for.truth$model)
for.truth$model_jjd=ifelse(for.truth$model=="Ricker","Ricker",for.truth$model_jjd)
for.truth$model_jjd=ifelse(for.truth$model=="segreg","BH",for.truth$model_jjd)
for.truth$model_jjd=ifelse(for.truth$model=="mean","RW",for.truth$model_jjd)

AIC.ALL=merge(AIC.ALL,for.truth,by="run")
AIC.ALL$correct_jjd=ifelse(AIC.ALL$model_jjd==AIC.ALL$Winner$Winner,1,0)
AIC.ALL$correct_om=ifelse(AIC.ALL$model==AIC.ALL$Winner$Winner,1,0)

#table(AIC.ALL$model_jjd,AIC.ALL$correct_jjd)
#table(AIC.ALL$model,AIC.ALL$correct_om)

model.winner=table(AIC.ALL$model,AIC.ALL$Winner$Winner)
rownames(model.winner)=paste0(rownames(model.winner),".om")
model.winner.perc=round(model.winner/rowSums(model.winner),3)

model.winner.jjd=table(AIC.ALL$model_jjd,AIC.ALL$Winner$Winner)
rownames(model.winner.jjd)=paste0(rownames(model.winner.jjd),".om")
model.winner.perc.jjd=round(model.winner.jjd/rowSums(model.winner.jjd),3)

table(AIC.ALL$model,AIC.ALL$Winner$Winner,AIC.ALL$sigmaR)

