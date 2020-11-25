//CRAN R code to analyze disruption of plant-pollinator networks for the article: Disruption of biological processes in the Anthropocene: the case of phenological mismatch
//Author Maël Montévil
//Cite as Montévil, M. 2020, _code for: Disruption of biological processes in the Anthropocene: the case of phenological mismatch_ DOI: 10.5281/zenodo.4290412


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp14)]]
//[[Rcpp::plugins(unwindProtect)]] 


//convenience and optimization functions

// [[Rcpp::export]]
LogicalMatrix submat_rcpp(LogicalMatrix X, LogicalVector condition, bool row) { 
    if(is_true(all(condition))){
        return(X);
    }
    
    int n=X.nrow(), k=X.ncol();
    if(row){
        LogicalMatrix out(sum(condition),k);
        for (int i = 0, j = 0; i < n; i++) {
            if(condition(i)) {
                out(j,_) = X(i,_);
                j = j+1;
            }
        }
        return(out);
    }else{
        LogicalMatrix out(n,sum(condition));
        for (int i = 0, j = 0; i < k; i++) {
            if(condition(i)) {
                out(_,j) = X(_,i);
                j = j+1;
            }
        } 
     return(out);
    }
}

// [[Rcpp::export]]
int fastSum(IntegerVector pm1) {
    int n = pm1.size();
    int res=0;
    for(int i = 0; i < n; ++i) {
        res+=pm1(i);
    }
    return res ;
}

 IntegerVector do_conc_(IntegerVector x, IntegerVector y){
    int nx=x.size(), n=x.size()+y.size(),i,j;
    IntegerVector out=no_init(n);
    for (i=0; i<nx; ++i){ 
        out[ i ] = x[ i ];
    }
   for (j=i, i=0; j<n; ++j, ++i){
        out[ j ] = y[i] ;
    }
   return out;
}


void randomThrs(int nn, int minn, int maxx,NumericVector rand,IntegerVector res) {
     rand=runif(nn,minn,maxx);
     for(int i = 0; i < nn; ++i) {
        res[i]=std::floor(rand[i]);
    }   
}

IntegerVector vfloor(NumericVector rand) {
    int nn=rand.size();
    IntegerVector res(nn);
    for(int i = 0; i < nn; ++i) {
        res[i]=std::floor(rand[i]);
    }   
    return res;
}
int fastbSum(LogicalVector vecc) {
    int res=0;
    for(int i = 0; i < vecc.size(); ++i) {
        if(vecc[i]){
            res++;
        }
    }   
    return res;
}


//data type combining phenology and interaction matrix
class webPheno{
    public:
        LogicalMatrix web;
        IntegerVector startp;
        IntegerVector endp; 
        IntegerVector startb;
        IntegerVector endb;

    webPheno(){}
    
    webPheno( LogicalMatrix web_,    IntegerVector startp_,    IntegerVector endp_,    IntegerVector startb_,    IntegerVector endb_ ):
        web(web_), startp(startp_), endp(endp_) , startb(startb_) , endb(endb_)  {}
        
    
 
    void removeb(LogicalVector remainb){
        web=submat_rcpp(web,remainb,false);
        startb=startb[remainb];
        endb=endb[remainb];  
    }

    void removep(LogicalVector remain){
        web=submat_rcpp(web,remain,true);
        startp=startp[remain];
        endp=endp[remain];
    }
       
    List returnL(){
        return(List::create(Named("web")=web, Named("startp")=startp, Named("endp")=endp, Named("startb")=startb,Named("endb")=endb));
    }
};


webPheno webPheno2(List ll) {
    return  webPheno(ll(0),(ll[1]),(ll[2]),(ll[3]),(ll[4]) );
}



//transform start/end, into the sum of presence every day
// [[Rcpp::export]]
IntegerVector day_presence(IntegerVector pm1,IntegerVector pm2) {
    int n = pm1.size();
    IntegerVector res(365*2);
    res=rep(0,365*2);

    for(int i = 0; i < n; ++i) {
        int  ll= pm2(i)-pm1(i);
        for(int j = 0; j < ll +1; ++j) {
            res[pm1(i) +j-1] += 1;
        }
    }
    for(int i = 0; i < 365; ++i) {
        res[i] +=res(i+365);
        res[365+i] =res(i);
    }
    return res ;
}


// [[Rcpp::export]]
double combineb(IntegerVector availp, int startb, int endb) {
    double res=0;

    for(int i = startb; i < endb+1; ++i) {
        res+= (availp(i-1)==0)*1;
    }
    return res/(endb-startb+1) ;
}




// [[Rcpp::export]]
LogicalVector day_presence_for(IntegerVector pm1,IntegerVector pm2, int startb, int endb) {
    int n = pm1.size();
    int len=endb-startb+1;
    LogicalVector  res(len);
    res.fill(true);

    int st1,st2, st, en;
    for(int i = 0; i < n; ++i) {
        // res(i) +=pm2[i];
        //int  ll= pm2(i)-pm1(i);
        st1=pm1[i];
        st2=pm2[i];
        st=(std::max(st1,startb) -startb);
        en=(std::min(st2,endb)-startb+1);
        for(int j =st; j < en; ++j) {
            res[j] =false;
        }
        st=(std::max(st1-365,startb) -startb);
        en=(std::min(st2-365,endb)-startb+1);
        for(int j =st; j < en; ++j) {
            res[j] =false;
        }
        st=(std::max(st1+365,startb) -startb);
        en=(std::min(st2+365,endb)-startb+1);
        for(int j =st; j < en; ++j) {
            res[j] =false;
        }
    }
    return res;
}

// [[Rcpp::export]]
double  ratio_day_presence_for(IntegerVector pm1,IntegerVector pm2, int startb, int endb) {
    LogicalVector res=day_presence_for( pm1, pm2,  startb,  endb);
    double res0=0;
    int n = res.size();
    for(int i = 0; i < n; ++i) {
        if (res[i]) {
            res0 ++;
        }
   }
    return res0/(endb-startb+1) ;
}



//returns changed web after agency effects
LogicalMatrix activity_(webPheno w1, NumericVector zeeross, double thres, double act) {
//    webPheno res=w1;
    int st1,en1,timee;
    IntegerVector pm1,pm2,res;
    int nb=w1.startb.size();
    int np=w1.startp.size();
  
    for(int j = 0; j < nb ; ++j) {
        if(!(zeeross[j]  <=thres)){
            if((runif(1, 0, 1))[0]<0.4){
                LogicalVector vp=w1.web(_, j);
                IntegerVector pm1=w1.startp[vp];
                IntegerVector pm2=w1.endp[vp];
                int sbj=w1.startb[j];
                LogicalVector res=day_presence_for( pm1, pm2,  sbj,  w1.endb[j]);
                for (int k = 0; k < res.size()  ; ++k) {
                    if(res[k]){
                        for (int i = 0; i <  np  ; ++i) {
                            if(!w1.web(i, j)){
                                timee=sbj +k;
                                st1=w1.startp[i];
                                en1=w1.endp[i];
                                //    if(!((std::min(w1.endp(i),w1.endb[j])-w1.startb[j]+1)<=(std::max(w1.startp(i),w1.startb(j)) -w1.startb(j))) ||
                                if(  ((st1<=timee) && (timee<=en1))||
                                        ((st1<=timee-365) && (timee-365<=en1))||
                                        ((st1<=timee+365) && (timee+365<=en1))){
                                    if((runif(1, 0, 1))[0]<act){
                                        w1.web(i, j)=true;
                                    }
                                }
                            }
                        } 
                    }
                }
            }
        }
    }
    return w1.web ;
}






struct iterate_res {
    NumericVector res;
    webPheno web;
    List varia;
};



//computes the result of a phenology change
iterate_res iterate_(webPheno w0, double thres, IntegerVector expop,IntegerVector expob ,webPheno w1, double  actb  , bool distest=false ) {
    List Lres;

    IntegerVector startp=w0.startp;
    IntegerVector startb=w0.startb;
    IntegerVector endp=w0.endp;
    IntegerVector endb=w0.endb;
    IntegerVector startp1, startb1,endp1,endb1;
    if(distest){
         startp1=w1.startp;
         startb1=w1.startb;
         endp1=w1.endp;
         endb1=w1.endb;
    }
    IntegerVector expopp=expop;
    IntegerVector expobb=expob;
    LogicalMatrix web;
    if(actb>0){    
         web=clone(w0.web);
    }else{
         web=w0.web;
    }
            
    int np = web.rows();
    int nb = web.cols();
    
    int dispb=0;
    int dispp=0;
    int dispb_temp=1;
    int dispp_temp=1;
    int stb, enb,enp,stp;
    double restemp;
    NumericVector res ={0,0,0,0,0};
    //iteration to draw all consequences of the phenological change
    while((dispp_temp>0))
    {
        if(nb>0)
        {
          NumericVector zeeross(nb);
            for(int i = 0; i < nb; ++i) {
                stb=startb[i];
                enb=endb[i];
                LogicalVector vp=web(_, i);
                restemp=sum(day_presence_for( startp[vp], endp[vp],  stb,  enb) );
                zeeross[i]= restemp/(enb-stb+1);
            }
                        
            
            if(actb>0){
                webPheno wt(web,startp,endp,startb,endb);
                web=activity_( wt,  zeeross,  thres,actb); 
                for(int i = 0; i < nb; ++i) {
                    stb=startb[i];
                    enb=endb[i];
                    LogicalVector vp=web(_, i);
                    restemp=sum(day_presence_for( startp[vp], endp[vp],  stb,  enb));
                    zeeross[i]= restemp/(enb-stb+1);
                }               
            }
            LogicalVector remainb= (zeeross  <=thres) ;
            
            expobb=expobb[remainb]; 
            web=submat_rcpp(web,remainb,false);

            startb=startb[remainb];
            endb=endb[remainb];
            if(distest){
                startb1=startb1[remainb];
                endb1=endb1[remainb];           
            }

            dispb_temp=sum(1-remainb);
            dispb+=dispb_temp;

            np = web.rows();
            nb = web.cols();
        }else{
            dispb_temp=0;
        }
        dispb_temp=0;
    
        if(np>0)
        {
            NumericVector zeeross(np);
            for(int i = 0; i < np; ++i) {
                stp=startp[i];
                enp=endp[i];
                LogicalVector vb=web(i,_);
                restemp=sum(day_presence_for( startb[vb], endb[vb],  stp,  enp));
                zeeross[i]= restemp/(enp-stp+1);
            }
            LogicalVector remain= (zeeross  <1) ;

       
            expopp=expopp[remain];
            web=submat_rcpp(web,remain,true);
            
            startp=startp[remain];
            endp=endp[remain];
            if(distest){
                startp1=startp1[remain];
                endp1=endp1[remain];                 
            }

            dispp_temp=sum(1- remain);
            dispp+=dispp_temp;
        
            np = web.rows();
            nb = web.cols();
        }else{
            dispp_temp=0;
        }
        dispp_temp=0;
    }
     
    res[0]=dispp;
    res[1]=dispb;
    res[2]=sum(expopp)+sum(expobb);
    res[3]=sum(1-expopp)+sum(1-expobb);
 
    webPheno wres;

    if(distest){
            Lres=List::create(  thres);
            wres=webPheno (web, startp1, endp1,  startb1, endb1);
    }else{
        res[4]=NA_REAL;
        Lres=List(res);   
    }
    
    iterate_res result;
    result.web=wres;  
    result.res= res;   
    result.varia=Lres;
    return  result;
}




 

NumericVector iterate_b(webPheno w0, double thres, IntegerVector expop,IntegerVector expob, double act){
    webPheno w1;
    struct iterate_res  Lres=iterate_ ( w0,  thres,  expop, expob , w1, act);
    NumericVector res=Lres.res;
    return(  res  );
} 

NumericVector iterate_c(webPheno w0, double thres, double act){
    webPheno w1;
    struct iterate_res Lres=iterate_ ( w0,  thres,  w0.startp*0, w0.startb*0 , w1, act);
    NumericVector res=Lres.res;
    return(  res  );
} 



// assess w_0^f 
double disappear_(webPheno w0, double thres, int j0 = 0 , int iterrm=10 ,bool onlyp=false, double act=0) {

    IntegerVector diffpp=w0.endp - w0.startp;
    IntegerVector diffbb=w0.endb-w0.startb;
    NumericVector resim (iterrm+1,-1.0);
    LogicalVector disapv (iterrm+1,false);
    
    double res=0;
    int itter=std::floor(iterrm /2);
    int  np = w0.web.rows();
    int nb = w0.web.cols();
    if((nb==0) | (np==0)){
        return(0);
    }
    double tempp=0;
    int size=0;
    double km, kk,diam,disap;
    int dpm,dbm;
    for(int j = 0; j < 182; ++j) {
        km=std::max(5,(20-j)*10+1);
        for (int k =0; k< km; ++k)
        {
            kk=k*0.1/(0.1*km);
            dpm=std::floor(np*kk);
            dbm=std::floor(nb*kk);
            diam=(dpm+dbm)*log((j+1)*2+1)+ (np-dpm+nb-dbm) * log(j*2+1) ;
            if((tempp - diam < -log(iterrm))|((kk==1)&(j==181)))
            {
                tempp=diam;
                size+=1;
            }
        }
    }

    IntegerVector jv(size);
    IntegerVector dpmv(size);
    IntegerVector dbmv(size);
    NumericVector diamv(size);
    size=0;
    tempp=0;
    int ij0=-1;
    for(int j = 0; j < 182; ++j) {
        km=std::max(5,(20-j)*10+1);
        for (int k =0; k< km; ++k)
        {
            kk=k*0.1/(0.1*km);
            dpm=std::floor(np*kk);
            dbm=std::floor(nb*kk);
            diam=(dpm+dbm)*log((j+1)*2+1)+ (np+nb-dpm-dbm) * log(j*2+1) ;
            if((tempp - diam < -log(iterrm))|((kk==1)&(j==181)))
                {
                        tempp=diam;
                        jv[size]=j;
                        dpmv[size]=dpm;
                        dbmv[size]=dbm;  
                        diamv[size]=diam;                       
                        size+=1;
                }
        }
        if(j==j0){
            ij0=(size-1);
        }
    }

    int incr=12;
    int bott=0;
    int topt=size-1;
    bool test1=true;
    int itt=std::floor((topt -bott)*0.1);
    int count=0;
    int count2=incr;

    IntegerVector startpt (np);
    IntegerVector startbt (nb);
    IntegerVector endpt (np);
    IntegerVector endbt (nb);
    IntegerVector shiftp (np);
    IntegerVector shiftb (nb);
    IntegerVector expob(nb);
    IntegerVector expop(np);

    IntegerVector rand11, rand12, randb11, randb12;
    NumericVector disap0;
    
    bool testt;
    while(test1)
    {
        if((size-itt)<1)
        {         
            Rcout << "The value of shiftp : " <<  act   << "\n";}
            testt=false;
            for(int ii = 0; ii < itter+1; ++ii) 
            {
                expob=do_conc_(rep(1,dbmv[itt]),rep(0,nb-dbmv[itt]));
                expop=do_conc_(rep(1,dpmv[itt]),rep(0,np-dpmv[itt])); 
                rand11 = floor(runif(dpmv[itt], -(jv[itt]+1),  jv[itt]+1+1));
                rand12 = floor(runif(np-dpmv[itt],  -jv[itt],  jv[itt]+1));   
                shiftp=do_conc_(rand11, rand12) ;

                startpt = w0.startp - shiftp ;
                
                for (int i=0; i<np; ++i){ 
                    startpt[ i ] = startpt[ i ] % 365 ;
                    int temp=startpt[ i ];      
                    if(temp<0){
                        startpt[ i ] =temp +365;
                    }
                }
                startpt=startpt+1;
                endpt = startpt+ diffpp;
                                                          
                randb11 =floor(runif(dbmv[itt], -(jv[itt]+1),  jv[itt]+1+1));
                randb12=floor(runif(nb-dbmv[itt],  -jv[itt],  jv[itt]+1));

                shiftb=do_conc_(randb11, randb12); 
                startbt = w0.startb - shiftb ; 
                
                for (int i=0; i<nb; ++i){ 
                    startbt[ i ] = startbt[ i ] % 365 ;
                    int temp=startbt[ i ];
                    if(temp<0){
                        startbt[ i ] =temp +365;
                    }
                }
                
                startbt=startbt+1;             
                endbt = startbt + diffbb;

                webPheno w1(w0.web,startpt,endpt,startbt,endbt);
                disap0=(iterate_b(w1,thres,expop,expob,act));
                diamv[itt]=(disap0[2])*log((jv[itt]+1)*2+1)+ (disap0[3]) * log(jv[itt]*2+1) ;
                            
                disap=disap0[0]+disap0[1];
                
                if(onlyp){  
                    disap=disap0[0];
                }
           
                if(disap ==0)
                {
                    if(count2<incr)
                    {
                        count2=incr-1;
                    }
                    testt=true ;
                    resim[count]=diamv[itt];
                    disapv[count]=true;
                    count+=1;
                    if(count>=iterrm)
                    {
                        count=0;
                    }
                }
        }
        if(testt){
            bott=itt;
        }else{
            topt=itt;
        }
               
        if(count2==incr)
        {
            int temp=itt;
            itt=bott+std::floor((topt -bott)*0.5);
            if(itt==temp){
                    count2=count2-1;
                    itt=std::min((itt+1),size-1);
            }                                    
        }else{
            test1= !(count2==0) & !(itt==size-1);
            count2=count2-1;
            itt=std::min((itt+1),size-1);
        }
    }

    NumericVector maxx= resim[disapv];

    res=max(maxx);
    res=res+1*log(sum((maxx > res- log(iterrm)))*1) -log(iterrm);
    res=res*res/res;
    
    return (res);
}



// [[Rcpp::export]]
NumericVector iterate_dist(LogicalMatrix web0,  IntegerVector startpt0,IntegerVector endpt0, IntegerVector startbt0,IntegerVector endbt0, double thres, IntegerVector expob,IntegerVector expop ,IntegerVector startp0,  IntegerVector endp0, IntegerVector startb0, IntegerVector endb0, int  j0=0 ,int iterrm=10, bool distest=false, bool onlyp=false, double act=0) {
    webPheno w0(web0,  startpt0, endpt0, startbt0, endbt0);
    webPheno w1(web0,  startp0, endp0, startb0, endb0);
    struct iterate_res Lres=iterate_( w0,  thres, expop, expob , w1,act ,  distest=distest); 
    NumericVector res=Lres.res;
    if(distest){
        res[4]=disappear_(Lres.web, thres, j0,  iterrm, onlyp=onlyp,act=act);
        return(res);
    }
    else{
       return(res);
    }
}

// [[Rcpp::export]]
NumericVector iterate_distb(SEXP web0,  SEXP startpt0,SEXP endpt0, SEXP startbt0,SEXP endbt0, double thres,SEXP startp0,  SEXP endp0, SEXP startb0, SEXP endb0, int  j0=0 ,int iterrm=10, bool distest=false, bool onlyp=false, double act=0) {
    webPheno w0(web0,  startpt0, endpt0, startbt0, endbt0);
    webPheno w1(web0,  startp0, endp0, startb0, endb0);
    struct iterate_res Lres=iterate_( w0,  thres,  w0.startp*0, w0.startb*0 , w1, act ,   distest=distest); 
    NumericVector res=Lres.res;
    if(distest){
        res[4]=disappear_(Lres.web, thres, j0,  iterrm, onlyp=onlyp,act=act);
        return(res);
    }
    else{
       return(res);
    }
}

// [[Rcpp::export]]
double iterate_distc(SEXP web0, double thres,SEXP startp0,  SEXP endp0, SEXP startb0, SEXP endb0, int  j0=0 ,int iterrm=10, bool onlyp=false, double act=0) {
    webPheno w1(web0,  startp0, endp0, startb0, endb0);
    return(disappear_(w1, thres, j0,  iterrm, onlyp=onlyp,act=act));
}




// [[Rcpp::export]]
NumericVector iterateb(SEXP web0,  SEXP startpt0,SEXP endpt0, SEXP startbt0,SEXP endbt0, double thres, SEXP expob,SEXP expop, double act=0){
    webPheno w0(web0,  startpt0, endpt0, startbt0, endbt0);
    return(  iterate_b( w0,  thres,  expop, expob,act));
} 

// [[Rcpp::export]]
NumericVector iteratec(SEXP web0,  SEXP startpt0,SEXP endpt0, SEXP startbt0,SEXP endbt0, double thres, double act=0){
    webPheno w0(web0,  startpt0, endpt0, startbt0, endbt0); 
    return(  iterate_c( w0,  thres,act=act)  );
} 

// [[Rcpp::export]]
List iterated( SEXP web0,  SEXP startpt0,SEXP endpt0, SEXP startbt0,SEXP endbt0, double thres, SEXP expob,SEXP expop ,IntegerVector startp0,  IntegerVector endp0, IntegerVector startb0, IntegerVector endb0, double act=0){
    webPheno w0(web0,  startpt0, endpt0, startbt0, endbt0);
    webPheno w1(web0,  startp0, endp0, startb0, endb0);
    bool distest=true;
    struct iterate_res  Lres=iterate_ ( w0,  thres,  expop, expob , w1, act,distest=distest);
    return( List::create(Named("res")=Lres.res, Named("webpheno")= Lres.web.returnL() ) );
} 


