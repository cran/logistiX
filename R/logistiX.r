logistiX<-function(x=NULL,y=NULL,strat=NULL,tst=1, sc=1,  pmid=2,option=2,noint=0, vstart=1, vstop=NA, alpha=0.05, details = 0)
{
 # x should be matrix of 0 and 1, y is vector of 0 and 1

 lr=0     ### currently not supported
 if(is.null(x) | is.null(y)) stop("ERROR: No x or y.\n")
 if(nrow(x) != length(y)) stop("ERROR: length of y does not match nrow(x).\n")
 for(i in 1:ncol(x)) {
  xtab<-sort(unique(x[,i]))
  if((xtab[1]!= 0 & xtab[2]!=1) | length(xtab) > 2) stop("ERROR: currently, program can only process binary x (0/1 coded).\n")
  }
 n<-length(y)

 if(is.null(strat)) strat<-rep(1,n)

 n_j<-strattab<-table(strat)
 n_by<-sum(n_j)
 maxstrat=length(strattab)
 IP<-ncol(x)+1
 if(is.na(vstop)) vstop<-IP-1
 IPv<-(vstop-vstart+1)

estout<-(1:(2*3*IPv))*0
nmethods<-(tst+sc+lr)*(1+(pmid==2))

ciout<-matrix(0,nmethods*IPv,8)
conddist<-numeric(0)

msa.ind<-t(cbind(x,y))
mode(msa.ind) <- "double"
mode(maxstrat) <- "integer"
mode(IP) <- "integer"
mode(noint) <- "integer"
mode(option) <- "integer"
mode(tst) <- "integer"
mode(sc) <- "integer"
mode(lr) <- "integer"
mode(pmid) <- "integer"
mode(vstart) <- "integer"
mode(vstop) <- "integer"
mode(n_j) <- "integer"
mode(n_by) <- "integer"
mode(estout)<- "double"
mode(ciout)<-"double"


para<-c(maxstrat,
IP,
n_j,
noint,
option,
tst,
sc,
lr,
pmid,
vstart,
vstop,
alpha)
#cat("\n", para,"\n\n")

mint<-vstart:vstop
maxt<-vstart:vstop
zeilenDist<-vstart:vstop
for(j in vstart:vstop){
  index<-j-vstart+1
  if(j !=0) {
    mint[index] <- max(sum(y) - sum(1-x[,j]), 0)
    maxt[index] <- min(sum(x[,j]), sum(y))
    } else {
    mint[index] <-  0
    maxt[index] <- n
    }

  zeilenDist[index] <- maxt[index]-mint[index]+1
  #cat("Variable ", j, "\n")
  #cat("  mint: ", mint[index], " maxt: ", maxt[index], "zeilenDist: ", zeilenDist[index], "\n")
  }
sumzeilenDist<-sum(zeilenDist)+vstop-vstart+1
spaltenDist <- IP-noint+1+1

distout<-matrix(0,sumzeilenDist,spaltenDist)
mode(distout) <- "double"

#double *by_data_array, int maxstrat, int IP, int noint, int options, int tst, int sc, int lr, int pmid,
#	int vstart, int vstop, double alpha, int * n_j, int n_by, double *ciout, double *estout, double *distout, int details


interim <- .C("xlr",
    as.double(msa.ind),
    as.integer(maxstrat),
as.integer(IP),
as.integer(noint),                                          
as.integer(option),
as.integer(tst),
as.integer(sc),
as.integer(lr),
as.integer(pmid),
as.integer(vstart),
as.integer(vstop),
as.double(alpha),
as.integer(n_j),
as.integer(n_by),
ciout=double(length(ciout)),
estout=double(length(estout)),
distout=double(length(distout)),
as.integer(details),
package="logistiX"
)

#estout


varMax <- length(interim$estout)/(3*2)  #3 Spalten, 2 x pro Variable
estout <- matrix(NA, 4*varMax, 3)


dreiVar <- 0
for (i in 1:(3*varMax)) {
 estout[i+dreiVar, 3] <- interim$estout[i*2]
  if ((i%%3)==0)
   dreiVar <- dreiVar+1
 }

 estmeth<-character(4*varMax)

for(i in 0:varMax-1)
{
 for(j in 1:4)
 {
  estout[4*i+j, 1] <- i+1
 }
 estmeth[4*i+1] <- "MUE"
 estmeth[4*i+2] <- "MLE"
 estmeth[4*i+3] <- "LX "
 estmeth[4*i+4] <- "CCFL"
}

estout <- data.frame(estout[,1], estmeth, estout[,3])
colnames(estout) <- c("varnum", "method", "estimate")

#ciout

ciout <- t(matrix(interim$ciout,8,nmethods*((vstop-vstart+1))))



method <- matrix(0, varMax*4, 1)
for(i in 0:varMax-1)
{
 method[4*i+1, 1] <- "TST"
 method[4*i+2, 1] <- "TST-Pmid"
 method[4*i+3, 1] <- "SC"
 method[4*i+4, 1] <- "SC-Pmid"
}

ciout <- data.frame(ciout[,1], method[,1], ciout[,2], ciout[,3], ciout[,4], ciout[,5], ciout[,6], ciout[,7], ciout[,8])
colnames(ciout) <- c("var", "method", "lower", "upper", "p-value (2-sided)", "p-value (LE)", "p-value (GE)", "chi2", "z")

ciout[ciout[,2] ==  "TST" | ciout[,2] ==  "TST-Pmid",8:9] <- NA
#distout

interim$distout = t(matrix(interim$distout, spaltenDist, sumzeilenDist))

#cat("counts:", counts, "\n")

tobs <- interim$distout[1,2:(ncol(interim$distout)-1)]
distout <- interim$distout[interim$distout[,ncol(interim$distout)] != -1,]
counts <- distout[,ncol(distout)]
varnum <- distout[,1]
indexsp<-2+distout[,1]
t.stat <- numeric()

for(i in 1:nrow(distout))
 t.stat[i] <- distout[i,indexsp[i]]

distout <- data.frame(varnum=varnum, t.stat=t.stat, counts=counts)
colnames(distout) <- c("varnum", "t.stat", "counts")

distout<-distout[distout$counts!=0,]
#distout[order(distout[,1], partial=order(distout[,2])),] <- distout
distout <- distout[order(distout$t.stat),]
distout <- distout[order(distout$varnum),]

####put results into list
if(option>=2) {
  res<-list(estout=estout, ciout=ciout, distout = distout, tobs = tobs, 
  # inter = interim$ciout, 
  call=match.call())
  }
else res<-list(distout=distout, tobs=tobs)

attr(res,"class")<-"logistiX"

if(option>=2) {
  for(j in 0:(varMax-1))  {
   # cat("j: ",j," j+1: ", j+1, "\n")
    beta.fc<-optimize(function(X,...) -like(X,...), interval=c(-20,20),res=res, var=j+1, firth=T)$minimum
    res$estout[4*j+4,3]<-beta.fc  
    }
    }

attr(res,"class")<-"logistiX"
res
}


