print.logistiX<-function(x,...)
{
  res<-x
  f<-cbind(res$estout,res$ciout)[,-4]
  colnames(f)[2]<-"method.est"
  colnames(f)[4]<-"method.ci"
  print(f)
  invisible(res)
}

summary.logistiX<-function(object, esttype="LX", citype="exact", cimethod=NULL, pmid=FALSE, testtype="TST", cilevel=0.95, ...){
  obj<-object
  estimates<-coef.logistiX(obj, type=esttype)
  if(is.null(cimethod)) cimethod<-""
#  if(citype != "scores") 
  ci<-confint.logistiX(obj, type=citype, method=cimethod, level=cilevel)
#  else if(!pmid) ci<-obj$ciout[obj$ciout[,2]=="SC",3:4]
#      else ci<-obj$ciout[obj$ciout[,2]=="SC-Pmid",3:4]
  pvalues<-test(obj, varnum=min(as.numeric(names(estimates))):max(as.numeric(names(estimates))),type=testtype)[,c(2:3,5)]
  cat("Exact logistic regression\n\nCall:\n")
  print(obj$call)
  cat("\nEstimation method:          ",esttype,
      "\nCI method:                  ",citype,cimethod,
      "\nTest method:                ",testtype)
      if(pmid)cat("-PMid")
      cat("\n")
  cat("\nSummary of estimates, confidence intervals and parameter hypotheses tests:\n\n")
  print(cbind(estimates,ci,pvalues) )
}

