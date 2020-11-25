



library(Rcpp)
library(compiler)
library("Brobdingnag")
library(foreach)
library(doMC) 
library(robustbase)
registerDoMC(coreNb)
enableJIT(1)

#compile and load cpp
sourceCpp(paste(sourcedir,'source.cpp',sep="" ))

#global variables
meanplant=mean(plantphen$endRob-plantphen$startRob)
sdplant=sd(plantphen$endRob-plantphen$startRob)
meanstartb=  mean(  beephen$start.now -beephen$Rob.start ,na.rm=T)
meanstartp=  mean( plantphen$startnow -plantphen$startRob   ,na.rm=T)
sdstartb=  sd( beephen$start.now -beephen$Rob.start ,na.rm=T)
sdstartp=  sd( plantphen$startnow -plantphen$startRob   ,na.rm=T)
 
 
meanspec<- function (x){
    return (mean(x, na.rm=T))
}
 
# compute initial phenologies
initialize0<- function (web0){
    web=web0>0
    np=dim(web)[1]
    nb=dim(web)[2]
    web=web[sample(1:np),sample(1:nb)]
  
    startp= floor(runif(np, min = min(plantphen$startRob), max = max(plantphen$startRob)))
    endp=startp + ceiling(abs(rnorm(np, mean =meanplant , sd = sdplant)))
    diffp=endp-startp
    startb=1:nb *0
    endb=1:nb *0
    zeeross=1:nb *0

    for (j in 1:nb){
        pm1=startp[web[,j]]
        pm2=endp[web[,j]]

        minmin=min(pm1)
        minmax=min(max(pm1),min(pm2)-1)

        maxmax=max(pm2)
        maxmin=max(min(pm2),max(pm1)+1)

        res=day_presence(pm1,pm2)
        v1=(1):(minmax)
        v2=(maxmin):(maxmax+1)

        prob1=log(res[v1]+1)
        prob2=log(res[v2]+1)
        startb[j]=sample(v1, 1, prob=prob1)
        endb[j]=sample(v2, 1, prob=prob2)

        zeeross[j]=min(res[startb[j]:endb[j]])
    }

    for (j in 1:nb){
        pm1=startp[web[,j]]
        pm2=endp[web[,j]]
        if(length(pm1)==0){
            zeeross[j]=0
        }else{
            res=day_presence(pm1,pm2)
            zeeross[j]=min(res[startb[j]:(endb[j])])
        }
    }

    web=web[,zeeross>0]

    startb=startb[zeeross>0]
    endb=endb[zeeross>0]
    np=dim(web)[1]
    nb=dim(web)[2]
    if(length(np)>0){
        zeeross=1:np *0

        for (j in 1:np){
            pm1=startb[web[j,]]
            pm2=endb[web[j,]]
            res=day_presence(pm1,pm2)

            zeeross[j]=max(res[startp[j]:endp[j]])
        }

        web<-web[zeeross>0,]
        startp<-startp[zeeross>0]
        endp<-endp[zeeross>0]
    }else{
        web=c()
    }
    np=dim(web)[1]
    nb=dim(web)[2]
  
    diffb=endb-startb
    diffp=endp-startp

    return(list(web=as.matrix(web>0),startb=startb,endb=endb,startp=startp,endp=endp,np=np,nb=nb,diffb=diffb,diffp=diffp))  
}

initialize<-cmpfun(initialize0)


# compute initial phenologies with thres % of surviving species

initializeb0<- function (web0, thres=0.99){
            testt=TRUE
            while(testt){
            xxx1<-mclapply(1:coreNb,function(i0) { 
                        testt2=TRUE
                        itera=20
                        while(testt2 &(itera>0)){
                                abb=initialize(web0>0)
                                if(length(abb[["nb"]])>0){
                                    testt2=(abb[["nb"]]/dim(web0)[2] <thres)||((abb[["np"]]/dim(web0)[1] <thres))
                                }
                                itera=itera-1
                        }
                return(list(testt2,abb))
            },mc.preschedule=T,mc.cores = coreNb)
            
            for (i in xxx1){
                if(!i[[1]]){
                        abb=i[[2]]
                        testt=FALSE
                    }
                }
            }
            return(abb)
}
initializeb<-cmpfun(initializeb0)

##########################################

# macrostate description in inital space
   
entropize0<- function (abb, iterr2=10,iterrm=100, alpha=1,thres=0, testlist=FALSE,distest=FALSE, disiterr=10,destroy=T,act=0,listtarget=-1){
    web=abb[["web"]]
    startb=abb[["startb"]]
    endb=abb[["endb"]]
    startp=abb[["startp"]]
    endp=abb[["endp"]]
    np=abb[["np"]]
    nb=abb[["nb"]]
    diffb=abb[["diffb"]]
    diffp=abb[["diffp"]]
 
    variables=c("loss", "logwi", "logwiC","logWi","logwf0","logwf","logWf","distt","other")
 
    tempp=0
    resi=c()
    jjj=c()
    kkk=c()

    x1=rep(NaN,(nb+np+1)*iterr2)
    rescumbp=as.data.frame(cbind(x1,x1,x1,x1,x1,x1,x1,x1,x1))
    x1=rep(NaN,(nb+1)*iterr2)
    rescumb=as.data.frame(cbind(x1,x1,x1,x1,x1,x1,x1,x1,x1))
    x1=rep(NaN,(np+1)*iterr2)
    rescump=as.data.frame(cbind(x1,x1,x1,x1,x1,x1,x1,x1,x1))
    names(rescumbp) <- variables
    names(rescumb) <- variables
    names(rescump) <- variables
    rescumbp$loss=rep(0:(nb+np),times=iterr2)
    rescumb$loss=rep((0:nb),times=iterr2)
    rescump$loss=rep((0:np),times=iterr2)

    np=dim(web)[1]
    nb=dim(web)[2]
  
    for (j in 0:181){
        for (k in 1:max(5,(100-j)*10)){
            kk=k/max(5,(100-j)*10)
            dp=floor(np*kk)
            db=floor(nb*kk)
            diam=(dp+db)*log((j+1)*2+1)+ (np-dp+nb-db) * log(j*2+1) 
            if((tempp - diam <log(alpha) -log(iterrm))){
                resi=c(resi,rep(diam, iterrm))
                tempp<-diam
                jjj=c(jjj,rep(j, iterrm))
                kkk=c(kkk,rep(kk, iterrm))
            }
        }
    }

    resi=c(resi[1:(length(resi)-iterrm)],rep(diam, iterrm))
    jjj=c(jjj[1:(length(jjj)-iterrm)],rep(j, iterrm)) 
    kkk=c(kkk[1:(length(kkk)-iterrm)],rep(kk, iterrm))

    for (ll in 0:(iterr2-1)){
        print(ll)
        res=c()
        res0=c()
        gc()
        xxx1<-mclapply(1:length(jjj),function(i0) { 
            j=jjj[i0]
            kk=kkk[i0]
            dp=floor(np*kk)
            db=floor(nb*kk)

            shiftp=c(floor(runif(dp, min = -(j+1), max = j+1+1)),floor(runif(np-dp, min = -j, max = j+1))) 
            shiftb=c(floor(runif(db, min = -(j+1), max = j+1+1)),floor(runif(nb-db, min = -j, max = j+1)))

            expob=c(rep(1,db),rep(0,nb-db))
            expop=c(rep(1,dp),rep(0,np-dp))
            startpt=(startp-shiftp) %% 365 +1
            endpt=startpt+diffp
            hop=startp-startpt %% 365 +1
            largs=(hop)>182.5
            smalls=(hop)< -182.5
            shiftp=(hop) - 365 * largs +365*smalls

            startbt=((startb-shiftb)+destroy*182*(i0 %% 2 ==1)) %% 365 +1
            endbt=startbt+diffb

            hop=startb-startbt %% 365 +1
            largs=(hop)>182.5
            smalls=(hop)< -182.5
            shiftb=(hop) - 365 * largs +365*smalls
            starrrb=startb 
            endddb=endb
            starrrp=startp 
            endddp=endp
            if(destroy &0  & (i0 %% 2 ==1)){
                starrrb=startbt
                endddb=endbt
                starrrp=startpt 
                endddp=endpt
            }
            jjj0=max(0,j-1)

            if(distest){
                hopp=iterate_dist(web,startpt,endpt,startbt,endbt,thres,expob,expop,starrrp,endddp,starrrb,endddb,j0=jjj0,iterrm=disiterr,distest=distest,act=act)
            }else{
                if(testlist){
                    hopp0=iterated(web,startpt,endpt,startbt,endbt,thres,expob,expop,starrrp,endddp,starrrb,endddb ,act=act);
                    hopp=hopp0[[1]]
                }else{
                    hopp=iterateb(web,startpt,endpt,startbt,endbt,thres,expob,expop,act=act);
                }
            }

            if(testlist){
                tempp=hopp0[[2]]
                if(listtarget>=0){
                    if((((hopp[1]+hopp[2])!=0)&&(hopp[1]+hopp[2])!=listtarget)){
                        tempp=c()
                    }
                }
                        
                return(list(hopp[1]+hopp[2],jjj0 ,hopp[1],hopp[2],log((j+1)*2+1)*(db+dp) + log(j*2+1)*(np+nb-db-dp) ,tempp))
               }else{ 
                return(c(hopp,mean(abs(shiftp))*0.5+0.5*mean(abs(shiftb)),log((j+1)*2+1)*(hopp[3]) + log(j*2+1)*(hopp[4]) ))
            }
        },mc.preschedule=T,mc.cores = coreNb)

        if(testlist){
            return(list(xxx1)) 
        }
        xxx1<-simplify2array(xxx1)
        res=rbind(res,t(xxx1))
        res0=rbind(res0,res)

        res0<-as.data.frame(res0)
        names(res0)<-c("disp","disb","expoplus","expomoins","logwf","dist","diam")

        resii=(resi)

        reducediam=log((jjj+1)*2+1)*(res0$expoplus) + log(jjj*2+1)*(res0$expomoins) 
        reducediam2=log((181+1)*2+1)*(res0$expoplus+res0$expomoins) 

        for ( i in 0:(np+nb)){
            index0=ll*(np+nb+1)+ i+1
            trindex= res0$disp+res0$disb==i
            rescumbp$logwi[index0]=log(sum(exp(as.brob( tail(sort( resii[trindex]),iterrm*2))))/iterrm )
            wmean<-function(x){return( weighted.mean(x[trindex], exp(resii[trindex] - rescumbp$logwi[index0]), na.rm = TRUE))}
            rescumbp$logwf0[index0]=wmean(reducediam) 
            temp=wmean(reducediam2 )
            if(is.infinite(temp))
                {temp=NaN}
            rescumbp$logWf[index0]=temp
            rescumbp$logwf[index0]=wmean(res0$logwf)
            rescumbp$distt[index0]=wmean(res0$dist)
            rescumbp$other[index0]=wmean(np-res0$disp)
        }
        #corrections
         index0=ll*(np+nb+1)+ 0:(np+nb)+1
        xx=c(1/sort(colSums(web>0)),rowSums(web>0)*0)
        xxx=c(0)
        for(i in 1:(length(index0)-1)){
            xxx=c(xxx,sum(xx[1:i]) )
        }
        rescumbp$logwiC[index0]=(rescumbp$logwi[index0])*(np+nb-xxx)/(np+nb) + log((365))*xxx

        for ( i in 0:(np)){
            index0=ll*(np+1)+ i+1
            trindex= res0$disp==i
            rescump$logwi[index0]=log(sum(exp(as.brob( tail(sort( resii[trindex]),iterrm*2))))/iterrm )
            wmean<-function(x){return( weighted.mean(x[trindex], exp(resii[trindex] - rescump$logwi[index0]), na.rm = TRUE))}
            rescump$logwf0[index0]=wmean(reducediam) 
            temp=wmean(reducediam2 )
            if(is.infinite(temp))
                    {temp=NaN}
            rescump$logWf[index0]=temp
            rescump$logwf[index0]=wmean(res0$logwf)
            rescump$distt[index0]=wmean(res0$dist)
            rescump$other[index0]=wmean(nb-res0$disb)
        }
        #corrections

        for ( i in 0:(nb)){
            index0=ll*(nb+1)+ i+1
            trindex= res0$disb==i
            rescumb$logwi[index0]=log(sum(exp(as.brob( tail(sort( resii[trindex]),iterrm*2))))/iterrm )
            wmean<-function(x){return( weighted.mean(x[trindex], exp(resii[trindex] - rescumb$logwi[index0]), na.rm = TRUE))}
            rescumb$logwf0[index0]=wmean(reducediam) 
            temp=wmean(reducediam2 )
            if(is.infinite(temp))
                    {temp=NaN}
            rescumb$logWf[index0]=temp
            rescumb$logwf[index0]=wmean(res0$logwf)
            rescumb$distt[index0]=wmean(res0$dist)
            rescumb$other[index0]=wmean(np-res0$disp)
        }
        index0=ll*(nb+1)+ 0:(nb)+1
        xx=c(1/sort(colSums(web>0)))
        xxx=c(0)
        for(i in 1:(length(index0)-1)){
            xxx=c(xxx,sum(xx[1:i]) )
        }
        rescumb$logwiC[index0]=(rescumb$logwi[index0])*(nb-xxx)/(nb) + log((365))*xxx

    }

    rescumbp$logWi = (np+nb)*log(365)
    rescumb$logWi =(np+nb)*log(365)
    rescump$logWi =(np+nb)*log(365)
  
    resultatbp= as.data.frame(aggregate(rescumbp,by=list(rescumbp$loss),meanspec))
    resultatp  = as.data.frame(aggregate(rescump,by=list(rescump$loss),meanspec))
    resultatb  = as.data.frame(aggregate(rescumb,by=list(rescumb$loss),meanspec))
    resultatbp=resultatbp[,variables]
    resultatb=resultatb[,variables]
    resultatp=resultatp[,variables]
  
    return (list("resbp"=resultatbp,"raw"=res0,"diam"=resi,"resitbp"=rescumbp,
            "resp"=resultatp,"resb"=resultatb))  
}
entropize<-cmpfun(entropize0)


 

# computations in final space

entropizewf0<- function (abb, thres=0, iter1=5, iter2=5 , iterrm=100,target="bp",act=0,destroy=T){
    np=abb[["np"]]
    nb=abb[["nb"]]   
    web=abb[["web"]]    
    rescum0=c()

    if(target=="pb" | target=="bp"){
        onlyp=FALSE;
        nn=np+nb;
        indexx=1
    }
    if(target=="p"){
        onlyp=TRUE;
        nn=np;
        indexx=3
    }
    if(target=="b"){
        onlyp=FALSE;
        nn=nb;
        indexx=4
    }

    for(hop in 1:iter2){
        xxx00 =  entropize(abb, iterr2=1,iterr=iterrm, alpha=2,thres=thres, testlist=TRUE,distest=FALSE, disiterr=10,destroy=destroy,act=act)
        for (i in 1:10){gc(full=TRUE);}
        
        xxx00=xxx00[[1]]
        hop00=simplify2array(lapply(xxx00, 
            function(i) {
                return(i[[indexx]]) }))
        hop0=simplify2array(lapply(xxx00, 
            function(i) {
                return(i[[5]]) }))

        idd=c()
        for (j1 in 0: nn){
            if (sum(hop00==j1) > iter1) {
                xx=sort(hop0[hop00==j1],index.return = T )
                idd=c(idd, sample((1:length(hop00))[hop00==j1],iter1))
            }else{
                idd=c(idd, (1:length(hop00))[hop00==j1])
            }
        }

        xxx001<-xxx00[idd];
        rm(xxx00)
        for (aa in 1:10){gc(full=TRUE);}
        ressss<-t(simplify2array(mclapply(1:length(idd),
            function(j2) {
                i0= xxx001[[j2]][[2]]
                webpheno=xxx001[[j2]][[6]]
                hopp=iterate_distc(webpheno$web ,thres,webpheno$startp, webpheno$endp, 
                        webpheno$startb, webpheno$endb, j0=i0 ,iterrm=floor(iterrm*1),onlyp=onlyp,act=act)
                return(c(hopp, xxx001[[j2]][[indexx]], xxx001[[j2]][[5]] ))
            },mc.preschedule=T,mc.cores =coreNb)))

        for (aa in 1:10){gc(full=TRUE);}
        res22=c()
        for (j1 in 0: nn){
            if((act>0)&FALSE){
                res22=c(res22, quantile( ressss[ ressss[,2]==j1  ,1],na.rm=T)[4])
            }else{
                res22=c(res22, mean( ressss[ ressss[,2]==j1  ,1],na.rm=T))
            }
        }

        rescum0=rbind(rescum0,res22)
    }

    logwf=colMedians(rescum0,na.rm=T)
    Lloss=rep(0: nn,each=iter2)
    Llogwf =as.vector(rescum0)
    fulll=cbind(Lloss,Llogwf)
     
    return (list(logwf,fulll,ressss))  
}

entropizewf<-cmpfun(entropizewf0)

  

  
  
# computations in both initial and final space
#full: compute also for d limited to plants or pollinators
  
entropizefull<- function (abb, thres=0, iter2wf=5, iiter1wf=10, iter2main=20 , iterrm=20,act=0,full=TRUE,destroy=T){
    resultat=  entropize(abb, iterr2=iter2main,iterr=iterrm, alpha=2,thres=thres, testlist=FALSE,distest=F,act=act,destroy=destroy)
    for (aa in 1:10){gc(full=TRUE);}
    res1=  entropizewf (abb, thres=thres, iter1=iiter1wf, iter2=iter2wf , iterrm=iterrm,target="bp",act=act,destroy=destroy)
    resultat[["resbp"]]$logwf=res1[[1]]
    for (aa in 1:10){gc(full=TRUE);}
    if(full){
        res1=  entropizewf (abb, thres=thres, iter1=iiter1wf, iter2=iter2wf , iterrm=iterrm,target="b",act=act,destroy=destroy)
        resultat[["resb"]]$logwf=res1[[1]]
        for (aa in 1:10){gc(full=TRUE);}
        res1=  entropizewf (abb, thres=thres, iter1=iiter1wf, iter2=iter2wf , iterrm=iterrm,target="p",act=act,destroy=destroy)
        resultat[["resp"]]$logwf=res1[[1]]
        for (aa in 1:10){gc(full=TRUE);}
    }
    return (resultat)  
}
  

##########################################

  
# randomize phenologies (uniform (alter=F)  or gaussian (alter=T), with parameters based on data, modulated but the effect parameter ) )

entropizeb0<- function (abb, iterr2=1,iterrm=100, thres=0,  alter=FALSE, effect=c(1,1), act=0,distest=F){
    web=abb[["web"]]
    startb=abb[["startb"]]
    endb=abb[["endb"]]
    startp=abb[["startp"]]
    endp=abb[["endp"]]
    np=abb[["np"]]
    nb=abb[["nb"]]
    diffb=abb[["diffb"]]
    diffp=abb[["diffp"]]
 
    resi=cbrob(1,2)
    jjj=c()
    kkk=c()
    rescum0=c()
    rescum0b=c()
    rescum0bb=c()

    np=dim(web)[1]
    nb=dim(web)[2]
    jjj=1:iterrm
  
    for (lll in 1:iterr2){
        print(lll)

        res0=c()
        gc()
        xxx1<-mclapply(1:length(jjj),function(i0) { 
            if(!alter){
                j=182
                #(floor(runif(np, min = -1, max = 1+1))) *sample(c(rep(1,floor(np*kk)),rep(0,np-floor(np*kk))))
                shiftp=floor(runif(np, min = -(j+1), max = j+1+1)) 
                shiftb=floor(runif(nb, min = -(j+1), max = j+1+1))
            }else{
                shiftp= round((rnorm(np, mean =meanstartp*effect[1], sd = sdstartp*effect[2]))) 
                shiftb= round((rnorm(nb, mean = meanstartb*effect[1]  , sd = sdstartb*effect[2]))) 
            }

            startpt=(startp-shiftp) %% 365 +1
            endpt=startpt+diffp
            hop=startp-startpt %% 365 +1
            largs=(hop)>182.5
            smalls=(hop)< -182.5
            shiftp=(hop) - 365 * largs +365*smalls

            #startbt=((startb-shiftb)) %% 365 +1
            startbt=((startb-shiftb)) %% 365 +1
            endbt=startbt+diffb

            hop=startb-startbt %% 365 +1
            largs=(hop)>182.5
            smalls=(hop)< -182.5
            shiftb=(hop) - 365 * largs +365*smalls

            if(distest){
                hopp=iterate_dist(web,startpt,endpt,startbt,endbt,thres,expob,expop,starrrp,endddp,starrrb,endddb,j0=jjj0,iterrm=disiterr,distest=distest,act=act)
            }else{
                hopp=iteratec(web,startpt,endpt,startbt,endbt,thres,act)    
            }
            
            return(c(hopp,mean(abs(shiftp))*0.5+0.5*mean(abs(shiftb))))
        },mc.preschedule=T,mc.cores =coreNb)

        gc()

        res0=rbind(res0,t(simplify2array(xxx1)) )
        res0<-as.data.frame(res0)
        names(res0)<-c("disp","disb","expoplus","expomoins","logwf","dist")
        rescumb=c()

        for ( i in 0:(dim(web)[2]+dim(web)[1])){
            rescumb=c(rescumb,sum(res0$disp+res0$disb==i)/iterrm)
        }
        rescum0b=rbind(rescum0b,rescumb)
        rescum0=rbind(rescum0,c(mean(res0$disp+res0$disb)/(np+nb),mean(res0$disp)/(np),mean(res0$disb)/(nb),
            sd(res0$disp+res0$disb)/(np+nb),sd(res0$disp)/(np),sd(res0$disb)/(nb) ))
    } 

    densitymax =colMeans(rescum0b,na.rm=T)
    disappeard =colMeans(rescum0,na.rm=T)

    loss=0:(dim(web)[2]+dim(web)[1])
    resultat=as.data.frame(cbind(loss, densitymax))
    return (list(disappeard,rescum0,res0))  
}

entropizeb<-cmpfun(entropizeb0)

 
  
  
    
    
###################time series analysis################

# computes for all networks in the data base
# same parameters than for entropizeb
#iter2 : number of iteration per network; randomize: original network or null model
entropizeseries2<- function ( iter2=100,randomize=FALSE, thres=0, iterr=10000,act=0,effectm=1,effectv=1,alter=F,alterand=F){
    
    conf=data.frame(c(1))
    conf$iter2=iter2
    conf$randomize =randomize
    conf$thres=thres
    conf$iterr=iterr
    conf$act  =act
    conf$effectm   =effectm
    conf$effectv   =effectv
    conf$alter   =alter
    if(alterand){randt="r2_"}else{randt=""}
      
    name=paste(outdir,"rand_",randt,randomize,"_thres_",thres,"_alter_",alter,"_act_",act,"_em_",effectm,"_ev_",effectv,".csv",sep="");
    name2=paste(outdir,"rand_",randt,randomize,"_thres_",thres,"_alter_",alter,"_act_",act,"_em_",effectm,"_ev_",effectv,"_settings.csv",sep="");
    write.csv(conf,  name2)

    resentrop=c()
    for (i in 1:dim(refs)[1]){
        robertson<- read.csv(file =paste(datadir, refs$ID[i],".csv",sep="" ),header=T);
        r1=c()
        for (j in 1:iter2){
            testt=TRUE
            while(testt){
                robertson2=robertson
                if (randomize){
                    if(alterand){
                        robertson2=(r2dtable(1,rowSums(robertson>0),colSums(robertson>0)))[[1]]
                    }else{
                        robertson2=(shuffle.web(robertson>0,N=1, legacy=F))[[1]]
                    }
                }else{
                robertson2=robertson
                }

                abb=initializeb(robertson2>0)
 
                if(length(abb[["nb"]])>0){
                    testt=abb[["nb"]]==0}
                }
 
            resultat=  entropizeb(abb, iterr2=1,iterr=iterr,thres=thres, alter=alter,effect=1*c(effectm,effectv),act=act)
            r1=rbind(r1,c(i,resultat[[1]]))
        }
        resentrop=rbind(resentrop,r1)
    }
    resentrop=as.data.frame(resentrop)
    names(resentrop) <-c("x1","disappb","disapp","disapb","varpb","varp","varb")
    write.csv(resentrop,  name)
    return(resentrop);
}

# loads saved computation
entropizeseries2load<- function ( randomize=FALSE, thres=0,act=0,effectm=1,effectv=1,alter=F,alterand=F){
    if(alterand){randt="r2_"}else{randt=""}
    name=paste(outdir,"rand_",randt,randomize,"_thres_",thres,"_alter_",alter,"_act_",act,"_em_",effectm,"_ev_",effectv,".csv",sep="");
    rescum0=read.csv(file = name,header=T)
    return(rescum0)
}
    
# loads saved computation combining original network and null model
entropizeseries2loadcombined<- function (randomize=FALSE, thres=0,act=0,effectm=1,effectv=1,alter=F,alterand=F,refs){
    rescum0= entropizeseries2load ( randomize=F, thres=thres, act=act,effectm=effectm,effectv=effectv,alter=alter& T,alterand=alterand& T)
    rescum0r= entropizeseries2load (randomize=T, thres=thres,act=act,effectm=effectm,effectv=effectv,alter=alter& T,alterand=alterand& T)
    dataset0=aggregate(rescum0,by=list(rescum0$x1),meanspec)
    
    if(randomize){
        dataset0r=aggregate(rescum0r,by=list(rescum0r$x1),meanspec)
        dataset0=dataset0-dataset0r
    }
    conn=refs$Connectance
    lat=abs(refs$Latitude)
    time=refs$date2
    dataset0=as.data.frame(cbind(dataset0,conn,lat,time))

    climate=ifelse(abs(refs$Latitude)<=23,"Tropical","temperate")
    climate=(ifelse(abs(refs$Latitude)<60,climate,"Cold"))
    climate=as.factor(climate)

    dataset0$climate=climate
    timee=time>=1990
    dataset0$timee=timee
  
    return(dataset0)
} 
    
    
    
