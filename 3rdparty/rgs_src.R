###########################
#to generate R functions for RGS method
#copy all of the followings and paste on R console 
###########################

##
#Beginning of functions
###########################
########################### 


f.runEcen<-function(x,k=41,p=.2,pc=.05) 
{
#x: ISIs of spike train 
#Moving extreme center of pc-quantile with window size max(k,length(x)*p) 
#output: e-center of log transformed isi

#k: minimum window size
#p: window size (100*p=percentage of data point in a window)
#pc: high/low extreme percentage 

#Example: f.runEcen(x=spiketrain)

	f.Extcenter<- function(x,pc=pc){ #begin Extcenter

	#center of top pc and bottom pc of log transformed data
	
	mean(quantile(log(x),  probs = c(pc, 1-pc),type=8))

	#median unbiased estimate
	} #end Extcenter

 n <- length(x)
 k <- max(k,length(x)*p)
 k <- as.integer(k)
 res=rep(NA,n)

    if (k < 0) 
        stop("'k' must be positive")
    if (k%%2 == 0) 
        warning("'k' must be odd!  Changing 'k' to ", k <- as.integer(1 + 
            2 * (k%/%2)))
    if (k > n) 
        warning("'k' is bigger than 'n'!  Changing 'k' to ", 
            k <- as.integer(1 + 2 * ((n - 1)%/%2)))

kd=(k-1)/2 #half window size
	if (n==k) res[1:n]=f.Extcenter(x,pc=pc)  
   	
	if (n>k) {  #L1
#end rule constant
	tx1=f.Extcenter(x[1:k],pc=pc)
	tx2=f.Extcenter(x[(n-k+1):n],pc=pc)
	res[1:kd]=tx1
   	res[(n-kd+1):n]=tx2 

#moving center
   	for(i in (kd+1):(n-kd)) { #L2
	res[i]=f.Extcenter(x[(i-kd):(i+kd)],pc=pc)
	} #L2

		}   #L1

 return(res) #moving e-center of log(isi) 
	}  #end f.Extcenter
#############


#############
f.runmed<-function(x,k=41,p=.2){
#x: ISIs 
#Moving median with window size max(k,length(x)*p) 
#output: running median of the log transformed isi
#k: minimum window size
#p: window size (percentage of data point)

tx=runmed(x,k=max(k,floor(length(x)*p/2)*2+1),endrule="constant")
as.vector(tx)
} #end runmed
#############

######################
#modified normalization function
######################

f.norm<-
function(x,pc=.05,delta=0.1,...){

#this function generates the normalized log(isi) using
#three step central curve for log(isi)

#step 1: starting line e-center 
#step 2: central moving median of the isis who falls in centeral area
#from the e-center
#step 3: if difference between e-center and moving median is large
#another central moving median based on isis redefine central area 
#from the previous central moving median 
   
#input: x #isi data
#output: normalized ISI  (log(isi)-central curve)
#
#step 1 E-center
	tyc<-f.runEcen(x,pc=pc,...) #e-center of log(isi).

#step 2 moving median on central portiion
#standard deviation of normalized isi 
#and running median for +-sd*qnorm(.95)

	smad=mad(log(x)-tyc)
	lx=log(x)[abs(log(x)-tyc)<=qnorm(.95)*smad]
	tyc2=f.runmed(lx,k=41)

	#interpolate the tyc2 to all the area
	txs=seq(log(x))[abs(log(x)-tyc)<=qnorm(.95)*smad]
	txs=c(1,txs,length(x));tyc2=c(tyc2[1],tyc2,tyc2[length(tyc2)])
	tyc3=approx(txs,tyc2, xout=seq(x))$y # linear interpolation 

#Step 3 requires another run of step 2
#

		smad=mad(log(x)-tyc3)
		lx=log(x)[abs(log(x)-tyc3)<=qnorm(.95)*smad]
		tyc2=f.runmed(lx,k=41)

		#interpolate the tyc2 to all the area
		txs=seq(log(x))[abs(log(x)-tyc3)<=qnorm(.95)*smad]
		txs=c(1,txs,length(x));tyc2=c(tyc2[1],tyc2,tyc2[length(tyc2)])
		tyc3=approx(txs,tyc2, xout=seq(x))$y #  linear interpolation

return(log(x)-tyc3)
} #end f.norm


#######################
#initial seeds of bursts
########################
f.Nb.seed<-
function(x,pc=0.05,central=b.Cont,thresh=qnorm(0.005),...){

#generates the burst seed intervals whose length <=thresh 
#x: ISIs 
#p: window size of running center or median
#pc: top or bottom proportions used in E-center
#central: data.frame to determind the central portion
#combined normalized ISIs from a group of spike trains
#must have the variable "norm.length"
 
#Output: burst seeds ISIs with normalized ISIs
#data.frame(name=name,id=te1,cluster=clus2,
#interval=interval,norm.length=tz[te1])

tz<-f.norm(x,pc=.05,...) #normalized interval length (log unit) 

central=central[,"norm.length"] #Combined normalized ISIs
thresh=thresh*mad(central)

#if the normalized interval length is less than threshold
#burst-seeds: the normalized ISIs less than threshold + center 

center=median(tz[abs(tz)<=abs(thresh)])
all=data.frame(id=seq(x),interval=x,norm.length=tz) #original data
burst=all[tz<=thresh+center,]

#if burst is null, smallest interval to the burst
if (sum(tz<=thresh+center)==0) burst=all[all$norm.length==min(all$norm.length),]

tdata=data.frame(burst,clusid=seq(length(burst$id)))
list(tdata,all)
} #End f.Nb.seed
###############################



#####################
#pause seed intervals
#####################
f.Np.seed<-
function(x,pc=0.05,central=a.control,thresh=qnorm(0.995),...){

#detect the pause seeds whose normalized ISI>=thresh*mad(norm.length) 
#x: ISIs 
#p: window size of running center or median
#pc: top or bottom proportions used in E-center
#central: data.frame to determind the central portion
#combined normalized ISIs from a group of spike trains
#must have the variable "norm.length"

#Output: tdata=data.frame(name=name,id=te1,cluster=clus2,
#interval=interval,norm.length=tz[te1])

tz<-f.norm(x,pc=.05,...) #normalized interval length (log unit) 
central=central[,"norm.length"]
thresh=thresh*mad(central)

#seeds: the normalized interval length > threshold + center

center=median(tz[abs(tz)<=thresh])
all=data.frame(id=seq(x),interval=x,norm.length=tz)#original data
pause=all[tz>=thresh+center,]

#if pause is null, largest interval to the pause
if (sum(tz>=thresh+center)==0) pause=all[all$norm.length==max(all$norm.length),]

tdata=data.frame(pause,clusid=seq(length(pause$id)))
list(tdata,all)
} #End f.Np.seed
#########################

##############
#P-values for the bursts
##############

f.Np.burst<-function(bcluster,mu=0,sigma=1){ # begin f.NP.burst

#calculates the p-value of the burst cluster bcluster
#input: bcluster: id, length, norm.length #burst candidate
#output: log(P.value)

#parameters of the central distribution
#mu=median(central) 
#sigma=mad(central)

#Probability of sum of q Normal ISIs being less than time.2
#sum of q Normal ISIs ~ N(mean=q*mu,sd=sigma*sqrt(q)) 
	bid1<-bcluster[,"id"]
	time.2=sum(bcluster[,"norm.length"]) #sum of normalized ISI
	q<-length(bid1)
	burst.p1=pnorm(time.2,mean=q*mu,sd=sigma*sqrt(q),lower.tail = T,log.p = TRUE)
	burst.p1
  }  #End f.Np.burst

##############
#P-values for the pauses
##############

f.Np.pause<-function(pcluster,mu=0,sigma=1){ #Begin f.Np.pause

#calculates the p-value of the pause cluster pcluster
#input: pcluster: id, length, norm.length
#cenral: combined ISIs 
#output: log(P.value)

#parameters of the central distribution
#mu=median(central) 
#sigma=mad(central)


#Probability of sum of q Normal ISIs being greater than time.2
#sum of q Normal ISIs ~ N(mean=q*mu,sd=sigma*sqrt(q)) 
	bid1<-pcluster[,"id"]
	time.2=sum(pcluster[,"norm.length"]) #sum of normalized ISI
	q<-length(bid1)

	pause.p1=pnorm(time.2,mean=q*mu,sd=sigma*sqrt(q),lower.tail = FALSE,log.p = TRUE)

	pause.p1
  } #End f.Np.pause

#############
#Burst and Pause detection by Robust Gassian Surprise
#############

f.Nbp.RGS<-
function(ab=ab,ap=ap,adata,central=a.control,cthresh=qnorm(.95),...){

#from the initial burst or pause list (ab or ap) this function extends bursts
#or pauses by adding ISIs that decreases P-value of burst or pause cluster
#ab: list of burst seed
#ap: list of pause seed
#adata: ISI data from the spike train 
#central: combined isi data.frame from which central parameters to be estimated
#the column norm.length is normalized isi

#initial cut-off points for the central distribution
	central1=central[,"norm.length"]
	central1=central1[abs(central1-median(central1))<cthresh*mad(central1)]

#parameters of the central distribution
	mu=median(central1) 
	sigma=mad(central1)

	burst=NULL
	N=dim(adata)[1]
	clusterno<-unique(ab$clusid)

if(length(clusterno)>0) { #initial cluster candidates not empty

	for (jj in clusterno) { #jj
	p1=1
	p2=0
	bid1<-ab[ab$clusid==jj,"id"]

#extend to left
	while (p2<p1&&bid1[1]>1) { #
	p1<-f.Np.burst(bcluster=adata[bid1,],mu=mu,sigma=sigma)
	bid2<-c(bid1[1]-1,bid1)
	p2<-f.Np.burst(bcluster=adata[bid2,],mu=mu,sigma=sigma)
	if (p2<p1) bid1<-bid2
	} #

	p1=1
	p2=0
#extend to right
	while (p2<p1&&bid1[length(bid1)]<N) { #
	p1<-f.Np.burst(bcluster=adata[bid1,],mu=mu,sigma=sigma)
	bid2<-c(bid1,bid1[length(bid1)]+1)
	p2<-f.Np.burst(bcluster=adata[bid2,],mu=mu,sigma=sigma)

	if (p2<p1) bid1<-bid2
	} #

#P-value of each cluster
	pv1<-f.Np.burst(bcluster=adata[bid1,],mu=mu,sigma=sigma)
	burst0<-data.frame(adata[bid1,],clusid=clusterno[jj],P=pv1)
	burst<-rbind(burst,burst0)
	} #jj 

}  #cluster candidates not empty  

	pause=NULL
	N=dim(adata)[1]
	clusterno<-unique(ap$clusid)

if(length(clusterno)>0) { #for non empty pause set

	for (jj in clusterno) { #jj
	p1=1
	p2=0
	bid1<-ap[ap$clusid==jj,"id"]

#extend to left
	while (p2<p1&&bid1[1]>1) { #
	p1<-f.Np.pause(pcluster=adata[bid1,],mu=mu,sigma=sigma)

	bid2<-c(bid1[1]-1,bid1)
	p2<-f.Np.pause(pcluster=adata[bid2,],mu=mu,sigma=sigma)

	if (p2<p1) bid1<-bid2
	} #

	p1=1
	p2=0
#extend to right
	while (p2<p1&&bid1[length(bid1)]<N) { #
	p1<-f.Np.pause(pcluster=adata[bid1,],mu=mu,sigma=sigma)

	bid2<-c(bid1,bid1[length(bid1)]+1)
	p2<-f.Np.pause(pcluster=adata[bid2,],mu=mu,sigma=sigma)

	if (p2<p1) bid1<-bid2
	} #
	
#P-value
	pv2<-f.Np.pause(pcluster=adata[bid1,],mu=mu,sigma=sigma)
	pause0<-data.frame(adata[bid1,],clusid=clusterno[jj],P=pv2)
	pause<-rbind(pause,pause0)
	} #jj

}  #for non empty pause set 



#eliminate overlapped burst-clusters with larger P values 
	clusterno<-unique(ab$clusid)
	if(length(clusterno)>0) { #non empty initial cluster candidates

	bp2=burst
	for (w in 1:20) { #w: 20 could be changed to any number >10 
	clid=unique(bp2$clusid)
	for(ii in 1:(length(clid)-1)) { #ii
	d1<-bp2[bp2$clusid==clid[ii],"id"]
	d2<-bp2[bp2$clusid==clid[ii+1],"id"]
	overlap<-intersect(d1,d2)
	if (length(overlap)>=1) {  
	if(mean(bp2[bp2$clusid==clid[ii],"P"])<= mean(bp2[bp2$clusid==clid[ii+1],"P"]))
	bp2<-bp2[bp2$clusid!=clid[ii+1],] else bp2<-bp2[bp2$clusid!=clid[ii],]
			}
		} #ii
	} #w

	burst<-bp2
	burst$S=-burst$P
	burst$P=exp(burst$P)

  } # non empty initial cluster candidates

#eliminate overlapped pause-clusters with larger P values
	clusterno<-unique(ap$clusid)
	if(length(clusterno)>0) { #non empty pause set  

	bp2=pause
	for (w in 1:20) {  #w
	clid=unique(bp2$clusid)
	for(ii in 1:(length(clid)-1)) { #ii
	d1<-bp2[bp2$clusid==clid[ii],"id"]
	d2<-bp2[bp2$clusid==clid[ii+1],"id"]
	overlap<-intersect(d1,d2)
	if (length(overlap)>=1) {  
	if(mean(bp2[bp2$clusid==clid[ii],"P"])<= mean(bp2[bp2$clusid==clid[ii+1],"P"]))
	bp2<-bp2[bp2$clusid!=clid[ii+1],] else bp2<-bp2[bp2$clusid!=clid[ii],]
			}
		} #ii
	} #w

	pause<-bp2

	pause$S=-pause$P
	pause$P=exp(pause$P)
  } #non empty pause set 


#select bursts with P<0.05 and rename cluster id
	burst=burst[burst$P<=0.05,]
	pause=pause[pause$P<=0.05,]
	if (dim(burst)[1]>0) burst$clusid<-rep(seq(length(unique(burst$clusid))),table(burst$clusid))
	if (dim(pause)[1]>0) pause$clusid<-rep(seq(length(unique(pause$clusid))),table(pause$clusid))

	list(burst=burst,pause=pause)
 } #f.Nbp.RGS

#################################

##################################
#This function summarizes the bursts and pauses for a group of spike trains
#Group of Spike trains (ISIs) are read and stored as list object
#Output is list of bursts and pauses of each spike trains
###################################
   
f.BPsummary<-function(data=SIM1,mm=length(data),thresh1=qnorm(0.95),thresh0=qnorm(0.05)
,k0=41,p0=0.2,pc0=0.05,Pthresh=0.01) { #f.BPsummary Begin

#This function summarizes the bursts and pauses for a group of spike trains
#generated by RGS method

#INPUT: data: list of spike trains (isi intervals) 
#data[[kk]] is a data frame that has the column "isi"
#Output: list of burst and pauses
#bpout=f.BPsummary(data=SIM1)
#bpout$burst[[kk]]: bursts of spike train kk, SIM1[[kk]]
#bpout$pause[[kk]]: pauses of spike train kk, SIM1[[kk]]	

	Sim=data
#generate combine normalized isi 
	b.Cont<-NULL
	for (kk in 1:mm){ #kk=1
	zd=Sim[[kk]]$isi

	#normalized isi
	tnorm=f.norm(zd,k=k0,p=p0,pc=pc0)
	b.df<-data.frame(isi=zd,norm.length=tnorm)
	b.Cont<-rbind(b.Cont,b.df)
	}

#output format
	Sim.burst=list()
	Sim.pause=list()

	for(kk in 1:mm){	#kk
	zd=Sim[[kk]]$isi

#initial bursts and pauses calculation

	tbmed1.3<-f.Nb.seed(zd,thresh=thresh0,central=b.Cont,pc=pc0,k=k0,p=p0)
	ab2=tbmed1.3[[1]] #burst seeds for zd 
	adata=tbmed1.3[[2]] #all data 
	tpcen1.3<-f.Np.seed(zd,thresh=thresh1,central=b.Cont,pc=pc0,k=k0,p=p0)
	ap2=tpcen1.3[[1]] #pause seeds

#RGS bursts and pauses detection

	nbp2=f.Nbp.RGS(ab=ab2,ap=ap2,adata,central=b.Cont,cthresh=qnorm(0.95))

#initial bursts and pauses with P < 0.05
	bp.b=nbp2[[1]] 
	bp.p=nbp2[[2]]

	n=length(zd)+1

#number of bursts and pauses for Bonferroni adjustment
	Bfactor=length(unique(bp.b$clusid))
	Pfactor=length(unique(bp.p$clusid))
	#Bfactor=dim(ab2)[1]
	#Pfactor=dim(ap2)[1]

#Bonferroni adjustment
	bp.b$adjP=ifelse(bp.b$P*Bfactor>=1,1,bp.b$P*Bfactor)
	bp.p$adjP=ifelse(bp.p$P*Pfactor>=1,1,bp.p$P*Pfactor)

#Surprise value for adjusted P
	bp.b$adjS=-log(bp.b$adjP)
	bp.p$adjS=-log(bp.p$adjP)

	alpha=Pthresh  #definition of significance

	tn=names(table(bp.b$clusid)) 
	Csize=as.vector(table(bp.b$clusid))
	TN=as.numeric(tn[Csize>=1]) #correct to get cluster with one interval
#One may change the rule such as more than 1 for example 
#TN=as.numeric(tn[Csize>=2])
#TN=as.numeric(tn[Csize>=3])

#Select bursts and pauses with adjusted P-value < alpha
	bp.p$S=bp.p$adjS;bp.p$P=bp.p$adjP
	bp.b$S=bp.b$adjS;bp.b$P=bp.b$adjP 

	if(length(is.element(bp.b$clusid,TN)&bp.b$adjP<=alpha)>0) {
	bp.b=bp.b[is.element(bp.b$clusid,TN)&bp.b$adjP<=alpha,]
	}

	tn=names(table(bp.p$clusid)) 
	Csize=as.vector(table(bp.p$clusid))
	TN=as.numeric(tn[Csize>=1]) 
	if(length(is.element(bp.p$clusid,TN)&bp.p$adjP<=alpha)>0) {
	bp.p=bp.p[is.element(bp.p$clusid,TN)&bp.p$adjP<=alpha,]
	}

	tzd=cumsum(c(0,zd[-length(zd)])) ##spike time

#other summaries
	bp.b$start=tzd[bp.b$id]
	bp.b$end=tzd[bp.b$id]+bp.b$interval
	bp.p$start=tzd[bp.p$id]
	bp.p$end=tzd[bp.p$id]+bp.p$interval

	if (dim(as.data.frame(bp.b))[1]>0) {
	bp.b=bp.b[,-c(5,6)]
	}
	if (dim(as.data.frame(bp.p))[1]>0) {
	bp.p=bp.p[,-c(5,6)]
	}

	if (dim(as.data.frame(bp.b))[1]>0){
	bp.b$clusid=rep(1:length(unique(bp.b$clusid)),table(bp.b$clusid))
	}
	if (dim(as.data.frame(bp.p))[1]>0) {
	bp.p$clusid=rep(1:length(unique(bp.p$clusid)),table(bp.p$clusid))
	}

	Sim.burst[[kk]]=bp.b #bursts for spike train data[[kk]]
	Sim.pause[[kk]]=bp.p #pauses for spike train data[[kk]]

	} #kk

#output for data, a set of spike trains 
list(burst=Sim.burst,pause=Sim.pause)
} #End f.BPsummary

######################################
######################################
f.Combine<-function(Data=Sim,mm=length(Data)){

#this function generates normalized isi and combines isi from a set of spike trains
#input spike trains Data: List object
#Data[[kk]] isi vector for spike train kk
#or data frame that has isi component
#mm: the number of spike trains
#output: isi, normalized isi

	b.Cont<-NULL
	for (kk in 1:mm){ #kk=1
	zd=Data[[kk]]
	if(!is.vector(zd)) zd=Data[[kk]]$isi  
	k0=41
	p0=0.2
	pc0=0.05
	tnorm=f.norm(zd,k=41,p=.2,pc=.05)
	b.df<-data.frame(isi=zd,norm.length=tnorm)
	b.Cont<-rbind(b.Cont,b.df)
	} #kk
  b.Cont
 } #f.Combine

###########################################
###########################################

f.BPsummary2<-function(data=SIM1,mm=lenth(data),thresh1=qnorm(0.95),
thresh0=qnorm(0.05),k0=41,p0=0.2,pc0=0.05,Pthresh=0.01,
b.Cont=f.Combine(data0)) { #function Begin

#This function summarizes the bursts and pauses for a group of spike trains
#generated by RGS method while using the central distribution obtained from
#other data, say data0. 
#When data0=data, this function is the same as f.BPsummary

#INPUT: data: list of spike trains (isi intervals)
#Output: list of burst and pauses
#bpout=f.BPsummary(data=SIM1) 
#Combined isi and normalized isi for all
#May use other list
#bpout$burst[[kk]]: bursts of spike train kk, SIM1[[kk]]
#bpout$pause[[kk]]: pauses of spike train kk, SIM1[[kk]]	

	Sim=data

#output format

	Sim.burst=list()
	Sim.pause=list()

	for(kk in 1:mm){	#kk
	zd=Sim[[kk]]$isi


#initial bursts and pauses calculation
	tbmed1.3<-f.Nb.seed(zd,thresh=thresh0,central=b.Cont,pc=pc0,k=k0,p=p0)
	ab2=tbmed1.3[[1]] #bursts
	adata=tbmed1.3[[2]] #all data
	tpcen1.3<-f.Np.seed(zd,thresh=thresh1,central=b.Cont,pc=pc0,k=k0,p=p0)
	ap2=tpcen1.3[[1]] #pauses

	nbp2=f.Nbp.RGS(ab=ab2,ap=ap2,adata,central=b.Cont,cthresh=qnorm(0.95))

#initial bursts and pauses
	bp.b=nbp2[[1]] 
	bp.p=nbp2[[2]]

#number of bursts and pauses for Bonferroni adjustment
	Bfactor=length(unique(bp.b$clusid))
	Pfactor=length(unique(bp.p$clusid))
	#Bfactor=dim(ab2)[1]
	#Pfactor=dim(ap2)[1]
#Bonferroni adjustment
	bp.b$adjP=ifelse(bp.b$P*Bfactor>=1,1,bp.b$P*Bfactor)
	bp.p$adjP=ifelse(bp.p$P*Pfactor>=1,1,bp.p$P*Pfactor)

#Surprise value for adjusted P
	bp.b$adjS=-log(bp.b$adjP)
	bp.p$adjS=-log(bp.p$adjP)

	alpha=Pthresh  #definition of significance

	tn=names(table(bp.b$clusid)) 
	Csize=as.vector(table(bp.b$clusid))
	TN=as.numeric(tn[Csize>=1]) #correct to get cluster with one interval
#One may change the rule such as more than 1 for example 
#TN=as.numeric(tn[Csize>=2])
#TN=as.numeric(tn[Csize>=3])

#Select bursts and pauses with adjusted P-value < Threshold alpha
	bp.p$S=bp.p$adjS;bp.p$P=bp.p$adjP
	bp.b$S=bp.b$adjS;bp.b$P=bp.b$adjP 

	if(length(is.element(bp.b$clusid,TN)&bp.b$adjP<=alpha)>0) {
	bp.b=bp.b[is.element(bp.b$clusid,TN)&bp.b$adjP<=alpha,]
	}

	tn=names(table(bp.p$clusid)) 
	Csize=as.vector(table(bp.p$clusid))
	TN=as.numeric(tn[Csize>=1]) 
	if(length(is.element(bp.p$clusid,TN)&bp.p$adjP<=alpha)>0) {
	bp.p=bp.p[is.element(bp.p$clusid,TN)&bp.p$adjP<=alpha,]
	}

#other summaries
	tzd=cumsum(c(0,zd[-length(zd)])) #spike time
	bp.b$start=tzd[bp.b$id]
	bp.b$end=tzd[bp.b$id]+bp.b$interval
	bp.p$start=tzd[bp.p$id]
	bp.p$end=tzd[bp.p$id]+bp.p$interval

	if (dim(as.data.frame(bp.b))[1]>0) {
	bp.b=bp.b[,-c(5,6)]
	}
	if (dim(as.data.frame(bp.p))[1]>0) {
	bp.p=bp.p[,-c(5,6)]
	}

	if (dim(as.data.frame(bp.b))[1]>0){
	bp.b$clusid=rep(1:length(unique(bp.b$clusid)),table(bp.b$clusid))
	}
	if (dim(as.data.frame(bp.p))[1]>0) {
	bp.p$clusid=rep(1:length(unique(bp.p$clusid)),table(bp.p$clusid))
	}

	Sim.burst[[kk]]=bp.b
	Sim.pause[[kk]]=bp.p
	} #kk

   list(burst=Sim.burst,pause=Sim.pause)
 } #End f.BPsummary2



###########################
###########################
########################### 
#End of functions



