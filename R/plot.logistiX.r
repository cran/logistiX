## plot null and alternative distribution of sufficient statistic, compute (penalized) likelihood
plot.logistiX<-function(x,var,beta=0, plot=TRUE, firth=FALSE, ...){
    res<-x
    if(attr(res,"class") != "logistiX") stop("Object is not of class logistiX.\n")
    ci<-ncol(res$distout)
    v1<-res$distout[(res$distout[,1]==var) & (res$distout[,ci]!=-1),]
    if(nrow(v1)==0) stop("Variable ",var," is not in distribution.\n")
    tobs<-res$tobs
  
    score<-exp(beta * v1[,2]) * v1[,ci]
    score<-score/sum(score)
    score[order(v1[,2])]<-score
    t.stat<-v1[,2]
    t.stat[order(t.stat)]<-t.stat
    if(firth){
      e_t2<-t((t.stat)^2) %*% score
      et_2 <- (t(t.stat)%*% score)^2
      penalty <- sqrt(e_t2 - et_2)
    } else penalty<-1
    likeall<-score*penalty
    if(plot) barplot(score, names.arg=t.stat)
  
    return(list(dist=data.frame(t.stat=t.stat, score=score, penalized.score=likeall),likelihood=likeall[t.stat==tobs[1+var]]))
  }



like<-function(beta, res, var, firth){
 plot.logistiX(res, var, beta, plot=FALSE, firth=firth)$likelihood
 }


maximize<-function(res, var){
  optim(0, function(X,...) -like(X,...), res=res, var=var)$par
  }


## exact and mid-p null hypothesis tests for parameters
test<-function(obj, varnum, pmid=FALSE,  type="scores"){

 #method=c("scores","TST","probability") ### "PCLR" yet to be implemented!
 method<-type
 mind<-min(obj$distout$varnum)
 maxd<-max(obj$distout$varnum)
 if(min(varnum)<mind | max(varnum)>maxd) stop("Variables in output object (", mind, ":", maxd, ")do not match varnum argument (",varnum,").\n")
 if(pmid) pmidfactor<-0.5
 else pmidfactor<-1
 cardinality<-stat<-pval<-vn<-1:length(varnum)
 for(i in 1:length(varnum)){
   counts<-obj$distout$counts[obj$distout$varnum==varnum[i]]
   t.stat<-obj$distout$t.stat[obj$distout$varnum==varnum[i]]
   obs.stat<-obj$tobs[1+varnum[i]]
   if(method=="scores"){
     m.stat<-sum(t.stat*counts)/sum(counts)
     ss.stat<-sum(t.stat*t.stat*counts)/sum(counts)
     v.stat<-ss.stat-m.stat^2
     scores<-(t.stat-m.stat)^2/v.stat
     statistic<-(obs.stat-m.stat)^2/v.stat
     pvalue<-(sum(counts[scores>statistic])+sum(counts[scores==statistic]*pmidfactor)) /sum(counts)
   } else
   if(method=="TST"){
    tge<-(sum(counts[t.stat>obs.stat])+sum(counts[t.stat==obs.stat]*pmidfactor))/sum(counts)
    tle<-(sum(counts[t.stat<obs.stat])+sum(counts[t.stat==obs.stat]*pmidfactor))/sum(counts)
    smaller.tail<-min(tge,tle)
    pvalue=min(1,2*smaller.tail)
    statistic<-obs.stat
   } else
   if(method=="probability"){
     probs<-counts/sum(counts)
     statistic<-probs[obs.stat==t.stat]
     pvalue<-(sum(probs[probs<statistic])+sum(probs[probs==statistic]*pmidfactor))/sum(probs)
   }   else 
   if(method=="PCLR"){
     fulllike<-like(obj$estout[(varnum[i]-mind)*4+4,3], obj, varnum[i], firth=TRUE)
     nulllike<-like(0, obj, varnum[i], firth=TRUE)   
     statistic<-2*log(fulllike/nulllike)
     pvalue<-1-pchisq(statistic,1)
   } else stop("Method must be either scores, TST, PCLR or probability.\n")
   stat[i]<-statistic
   pval[i]<-pvalue
   vn[i]<-i
   cardinality[i]<-length(counts)
 }
 return(data.frame(varnum=vn,statistic=stat, pvalue=pval, method=method, cardinality=cardinality))
 }
 
 
coef.logistiX<-function(object, type="LX", ...){
  obj<-object
  matchi<-substr(as.character(obj$estout[,2]),1,2)==substr(type,1,2)   
  vstart<-min(obj$estout[,1])
  vstop<-max(obj$estout[,1])
  res<-(1:(max(obj$estout[,1])-min(obj$estout[,1])+1))*0
  names(res)<-min(obj$estout[,1]) : max(obj$estout[,1])
  res<-(obj$estout[matchi,3])
  res[res==999]<-sign(res[res==999])*Inf
  attr(res,"method")<-type
  names(res)<-vstart:vstop
  return(res)
}

coefficients.logistiX<-function(object, type="LX", ...){
  obj<-object
  matchi<-substr(as.character(obj$estout[,2]),1,2)==substr(type,1,2)   
  vstart<-min(obj$estout[,1])
  vstop<-max(obj$estout[,1])
  res<-(1:(max(obj$estout[,1])-min(obj$estout[,1])+1))*0
  names(res)<-min(obj$estout[,1]) : max(obj$estout[,1])
  res<-(obj$estout[matchi,3])
  res[res==999]<-sign(res[res==999])*Inf
  attr(res,"method")<-type
  names(res)<-vstart:vstop
  return(res)
}


cilimit<-function(obj, varnum, pmidfactor, low=TRUE, level, interval=c(-19,19), tol=0.00001){
      if(low) side<-"Lower"
      else side<-"Upper"
      tobs<-obj$tobs[1+varnum]
      objlike<-plot.logistiX(obj,varnum,0,plot=FALSE)
      # check if tobs is at border of distribution
      if(tobs==max(objlike$dist$t.stat) & !low) return(Inf)
      else  if(tobs==min(objlike$dist$t.stat) & low) return(-Inf)
      
      target<-function(beta, obj, varnum, pmidfactor, low, level,tobs){
          objlike<-plot.logistiX(obj,varnum,beta,plot=FALSE)
          if(!low)  target<-abs(max(sum(objlike$dist$penalized.score[objlike$dist$t.stat<tobs]),0)+objlike$dist$penalized.score[objlike$dist$t.stat==tobs]*pmidfactor - (1-level)/2)
          else target<-abs(max(sum(objlike$dist$penalized.score[objlike$dist$t.stat>tobs]),0)+objlike$dist$penalized.score[objlike$dist$t.stat==tobs]*pmidfactor - (1-level)/2)
          #cat(beta,target,"\n")
          return(target)
          }
      ciobj<-optimize(f=target, interval=interval, obj=obj, varnum=varnum, pmidfactor=pmidfactor, low=low, level=level,tobs=tobs)
      if(ciobj$objective>tol) {
        newint<-interval/2
        ciobj<-optimize(f=target, interval=newint, obj=obj, varnum=varnum, pmidfactor=pmidfactor, low=low, level=level,tobs=tobs)
        if(ciobj$objective>tol)       warning(side, " confidence limit of variable ", varnum, "not converged. Abs diff in confidence level =",ciobj$objective,"\n") 
        }
      cilim<-ciobj$minimum
      return(cilim)
    }

PCCL<-function(obj=obj, varnum=varnum, low=low, level=level, interval=interval, tol=0.0001){
  # implements the profile penalized completely conditional likelihood method proposed by Heinze and Puhr, SiM, 2010
  tobs<-obj$tobs[1+varnum]
  betamax<-coef.logistiX(obj,"CCFL")
  betamax<-betamax[names(betamax)==varnum]
  objlike<-plot.logistiX(obj,varnum,betamax,plot=FALSE, firth=TRUE)
  dist<-objlike$dist
  maxloglike<-log(objlike$likelihood)
  targetloglike<-maxloglike-qchisq(level,1)/2
  target<-function(beta,obj,varnum,targetloglike){
      return(abs(log(plot.logistiX(obj,varnum,beta,plot=FALSE, firth=TRUE)$likelihood)-targetloglike))
      }
  if(low) int<-c(interval[1],betamax)
  else int<-c(betamax,interval[2])
  cilimit<-optimize(target,   interval=int, obj=obj, varnum=varnum, targetloglike=targetloglike)
  if(cilimit$objective>tol) {
    if(low) side<-"lower"
    else side<-"upper"
    warning("PPCCL ", side, " confidence limit for variable ",varnum, " not converged, abs diff in confidence level=",cilimit$objective,".\n")
    }
  return(cilimit$minimum)
  }



confint.logistiX<-function(object, parm, level=0.95, type="exact", method="TST", ran.steps=seq(5/100,0.95,0.01), interval=c(-19,19), vstart=NULL, vstop=NULL, ...){
  # type=c("exact","pmid","randomized","PCCL")
  # method=c("TST","scores")
  # "exact" sets a factor of 1 for the observed value
  # "pmid" sets a factor of 0.5 for the observed value
  # "randomized" numerically integrates the confidence limit over pmid-factors from 0 to 1
  # ran.steps: steps for "randomized" method (the default is 0.05, 0.06, ..., 0.95 to compute a 10% trimmed mean to approximate the integrated randomized CL)
  # interval: interval to search for confidence limits
  
  ### some combinations not yet implemented, such as:
  ### randomized - scores
  
  obj<-object
  if(type=="randomized" & method=="scores") stop("Randomized method with scores not yet implemented.\n")
  
  
  if(type=="PCCL") method=""
  if(is.null(vstart)) vstart<-min(obj$ciout[,1])
  if(is.null(vstop))  vstop<-max(obj$ciout[,1])
  ci<-matrix(0,vstop-vstart+1,2)
  labs<-c(paste(100*(1-level)/2,"%"),paste(100*(1-(1-level)/2),"%"))
  colnames(ci)<-labs
  rownames(ci)<-vstart:vstop
  if(type=="pmid") pmidfactor<-0.5
  else if(type=="exact") pmidfactor<-1
  
  # call plot.logistiX and set beta such that the sum of penalized.score at t LE observed T is 0.025 or 0.975
  
  indv<-1
  for(varnum in vstart:vstop){
    for(low in c(TRUE, FALSE)){
      index<-1+!low
      if(method == "scores") {
       if(type != "pmid") ci[indv,index]<-obj$ciout[obj$ciout[,2]=="SC",3:4][indv,low==c(TRUE,FALSE)]
       else ci[indv,index]<-obj$ciout[obj$ciout[,2]=="SC-Pmid",3:4][indv,low==c(TRUE,FALSE)]
       }
      else {
          if(type=="exact" | type=="pmid") {
              ci[indv,index]<-cilimit(obj=obj, varnum=varnum, pmidfactor=pmidfactor, low=low, level=level,interval=interval)
              }
          if(type=="randomized")  {
            ci[indv,index]<-mean(sapply(seq(5/100,0.95,0.01),function(X) cilimit(obj=obj, varnum=varnum, pmidfactor=X, low=low, level=level, interval=interval)))
            }
          if(type=="PCCL") {
            ci[indv,index]<-PCCL(obj=obj, varnum=varnum, low=low, level=level, interval=interval)
            }
          }
      }
    indv<-indv+1
    }
   if(method!="") method<-paste("-",method,sep="")
   attr(ci,"CL method")<-paste(type, method, sep="")
   ci[ci==999]<-Inf
   ci[ci==-999]<--Inf
   return(ci)
  
  }

  
profile.logistiX<-function(fitted, parm=1, firth=FALSE, type="loglike", steps=101, betaseq=NULL, normalize=FALSE, ...){
    # type should be one of c("loglike", "like")
    obj<-fitted
    varnum<-parm
    if(length(varnum)>1) stop("Argument varnum should be scalar.\n")
    est<-coef.logistiX(obj)
     # determine suitable range by looking at exact CI minus 10% of the difference to the estimate
     ci<-confint.logistiX(obj, "PCCL",vstart=varnum, vstop=varnum)
     if(is.null(betaseq)){
      lo<-ci[1]-0.1*(est[names(est)==varnum]-ci[1])
      up<-ci[2]+0.1*(ci[2]-est[names(est)==varnum])
      stepsize<-(up-lo)/(steps-1)
      cat("Lo,up,stepsize", lo, up, stepsize,"\n")
      betaseq<-seq(lo,up,stepsize)
      }
     prof.tmp<-sapply(betaseq,function(X) plot.logistiX(obj,varnum,X, firth=firth, plot=FALSE)$likelihood)
     if(normalize) prof.tmp<-prof.tmp/max(prof.tmp)
     if(type=="loglike") prof.tmp<-log(prof.tmp)
     outmat<-cbind(betaseq,prof.tmp)
     colnames(outmat)<-c("beta","profile")
     attr(outmat,"class")<-"profile.logistiX"
     return(outmat)
     }   

  
  

