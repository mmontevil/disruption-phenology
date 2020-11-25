#CRAN R code to analyze disruption of plant-pollinator networks for the article: Disruption of biological processes in the Anthropocene: the case of phenological mismatch
#Author Maël Montévil
###download data 
#Phenology  https://doi.org/10.5061/dryad.rp321
#interaction webs www.web-of-life.es
#historical trends require a file references.csv of the network used, with the added column of time of data collection date2

#cNODF (not used below) https://doi.org/10.5061/dryad.dv1gq

##########################################
#initialize global variables and load functions
##########################################

options(browser = "firefox")
library(bipartite)
library(cocor)

#number of  CPU cores
coreNb=14

basedir="/home/kamome/Work/Articles/inprocess/disruption-phenology/analysis/disruption-phenology/"
sourcedir=basedir
datadir=paste( basedir,"data/",sep="") 
outdir=paste( basedir,"output/",sep="") 
dataRefs=paste( basedir,"data/references6.csv",sep="") 



refs<- read.csv(file =dataRefs,header=T)
beephen <- read.csv(file = paste(datadir,"dryad - bee phen.csv",sep=""))
plantphen<-read.csv(file =paste(datadir,"dryad - plant phen.csv",sep="" ))
source(paste(sourcedir,'toolbox.r',sep="" ))

################### Basic use#######################


###load and initiate network
#load interaction network
interactionweb<- read.csv(file =paste(datadir,"M_PL_017.csv",sep=""),header=T)

#initialize phenology
abb=initializeb(interactionweb>0)




###randomize phenologies

#number of samples
 n=10000
#randomize uniformly
resultat=  entropizeb(abb, iterr2=1,iterr=n,thres=0, alter=F,effect=c(1,1),act=0)
#randomize following a gaussian distribution (based on phenological data)
resultat=  entropizeb(abb, iterr2=1,iterr=n,thres=0, alter=T,effect=c(1,1),act=0)

#build a treemap of the results
resultat=resultat[[3]]
resultat$loss= resultat$disp+resultat$disb
resultat$count=1

library(treemap)
dataset000=aggregate(resultat,by=list(resultat$loss),sum)
dataset001=aggregate(resultat,by=list(resultat$loss),mean)
dataset001$count=dataset000$count/n
treemap(dataset001  ,index=c("loss"),vSize="count",sortID="loss", palette = "Reds",type="manual",vColor=c("count"), title.legend = "Proportion of the microscopic space", title = "Number of disappearing species" )
   
   
   

###compute properties of macrostates and disruptions
###makes computations only in initial space. thres: parameter R in the article, act: parameter a
resultat=  entropize(abb, iterr2=10,iterr=20,thres=0, testlist=F,destroy=T,act=0)

#plot entropy
plot(logwiC~loss,data=resultat[["resbp"]])
#plot specificity
plot(logWi-logwiC~loss,data=resultat[["resbp"]])

###compute properties of macrostates and disruptions
###makes computations in both initial and final space
resultat= entropizefull (abb, thres=0, iter2wf=10, iiter1wf=10, iter2main=10 , iterrm=20,act=0,full=F)
#final specificity
plot(logWf-logwf~loss,data=resultat[["resbp"]])
#specificity change
plot((logWf-logwf - (logWi-logwiC))~loss,data=resultat[["resbp"]])
#specificity per capita change
n=abb$np+abb$nb
plot(( (logWf-logwf)/(n -loss+0.0001) - (logWi-logwiC)/n)~ loss,data=resultat[["resbp"]])
#total specificity change
logwiC0=resultat[["resbp"]]$logwiC[(resultat[["resbp"]]$loss==0)]
plot((logWf-logwf - (logWi-logwiC0))~loss,data=resultat[["resbp"]])
#total specificity change per capita
plot(( (logWf-logwf)/(n-loss+0.0001) - (logWi-logwiC0)/n)~loss,data=resultat[["resbp"]])



     
###time series
#compute for original network (saves in a file)
mm= entropizeseries2 ( iter2=100,randomize=FALSE, thres=0, iterr=100,act=0,effectm=1,effectv=1,alter=T,alterand=T)
#compute for null model (saves in a file)
mmr= entropizeseries2 ( iter2=100,randomize=TRUE, thres=0, iterr=100,act=0,effectm=1,effectv=1,alter=T,alterand=T)

#load saved results 
rescum0= entropizeseries2load ( randomize=FALSE, thres=0, act=0,effectm=1,effectv=1,alter=T,alterand=T)
rescum0r= entropizeseries2load (randomize=TRUE, thres=0, act=0,effectm=1,effectv=1,alter=T,alterand=T)

#load saved results and combines them (original network - null model)
dataset0= entropizeseries2loadcombined ( randomize=T, thres=0, act=0,effectm=1,effectv=1,alter=F,alterand=T,refs)
       
plot( disappb ~time, dataset0)
x=lmrob( disappb ~lat* time, data=dataset0); summary(x)


