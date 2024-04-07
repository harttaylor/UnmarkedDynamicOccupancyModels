# ----------------------------------------------------------------------------
# Use AUC (psi estimates vs. actual site-year detections) to validate MSOMs
# ----------------------------------------------------------------------------
library(pROC)
names = c("BTNW")

for (i in names){
  #psi.estim<-read.csv(paste0("2_outputs/",i,"/",i,".preds"), header=T)
  #str(psi.estim)
  #need to change to matrix, transpose matrix, then stack the transposed columns?
 # psi.estimB<-as.matrix(psi.estim[,3:74])
  #str(psi.estimB)
  #psi.estim.stacked<-stack(psi.estim[,c(3:74)])
  #str(psi.estim.stacked)
  #occ.pred<-psi.estim.stacked$values 
  #write.csv(psi.estim.stacked, file=paste0("2_outputs/",i,"/",i,".psi.estim.stacked.csv"))
  occ.pred <- ifelse(P.occ>0.5, 1, 0)

  site_det<-read.csv(paste0("0_data/processed/species weights/",i,"/",i,"site_det.csv"), header=T)
  str(site_det)
  site_det.m<-as.matrix(site_det[,2:26])
  site_det.t<-t(site_det.m)
  site_det.tdf<-as.data.frame(site_det.t)
  site_det.stacked<-stack(site_det.tdf)
  str(site_det.stacked)
  occ.actual<-site_det.stacked$values 
  write.csv(site_det.stacked, file=paste0("0_data/processed/species weights/",i,"/",i,"site_det.stacked.csv"))
  
  validate.df<-data.frame(occ.pred, occ.actual)
  z<-roc(occ.actual, occ.pred)
  auc<-as.numeric(z$auc)
  
  tiff(paste0("2_outputs/",i,"/ROC_",i,".tiff"), units="in", width=12, height=8, res=300)
  plot.roc(z, xlim=c(1,0), main=paste0("AUC Pred. ",i," Occup. at Calling Lake, 1994-2016 (AUC=",round(auc,2),")"), asp=NA)
  dev.off()
}

