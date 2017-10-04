#CORES<-4#0 ##update for flux
JOBS<-100#00 ##update for flux
nodefile <- Sys.getenv("PBS_NODEFILE")
## are we on a cluster? ##nchar: count the number of characters
CLUSTER <- nchar(nodefile)>0
## if CLUSTER=-FALSE, assume we are on a multicore machine

require(doParallel)
CORES<-Sys.getenv("PBS_NUM_PPN")
if (CORES=="") CORES=4 else CORES=as.numeric(CORES)
registerDoParallel(CORES)

tic <- Sys.time()
mpar <- foreach(
  i=1:JOBS,
  .packages=c('pomp'),
  .inorder=FALSE) %dopar% 
  {
    Sys.sleep(i*.1)
    setwd("~/R")
    #setwd("~/pulse_1wk") 
    #setwd("~/Daphnia_pulse") ##turn on for flux
    NMIF<-200#0 ##update for flux
    NP<-10000#0 ##update for flux
    METHOD="mif2"
    source("Daphnia_1.R")
    param.tab <- read.table("params.csv",sep=",",row.names=1,header=TRUE)
    LV.pars <- c("taudem","tauenv","log.betaz1","log.betaz2","log.betaz3","log.betaz4","alpha", "mu","cc","eta","meanp","sigmap","sdbeta")
    LV.ivps <- c("prey.0")
    LV.rw.sd<- rw.sd(taudem=0.02, tauenv=0.02,log.betaz1=0.02,log.betaz2=0.02,log.betaz3=0.02,log.betaz4=0.02,alpha=0.02,
                mu=0.02,cc=0.02,eta=0.02,meanp=0.02,sigmap=0.02,sdbeta=0.02,prey.0=ivp(0.1))
    #LV.rw.sd <- rw.sd(rep(0.02,length(LV.pars)),rep(0.1,length(LV.ivps)))
    #names(LV.rw.sd) <- c(LV.pars,LV.ivps)
    
    LV.hyperparams <-
      list(min=unlist(param.tab["lower.bound",]),max=unlist(param.tab["upper.bound",]))
    
    LV.rprior <- function(hyperparams, ...)
    {
      r<-runif(length(hyperparams$min),min=hyperparams$min,max=hyperparams$max)
      names(r) <- names(hyperparams$min)
      return(r)
    }
    set.seed(8100+i)
    Sys.sleep(i*0.1)
    th.draw <-LV.rprior(LV.hyperparams)
    m<-try(mif2(Daphnia_pomp,
               Nmif=NMIF,
               start=th.draw, # we will initialize
               rw.sd=LV.rw.sd,
               Np=NP,
               cooling.type='geometric',
               cooling.fraction= 0.3,
               max.fail=200,
               transform=TRUE
    ))
    list(pomp=m,start=th.draw)
  }
#setwd("~/Google_Drive/Bytheffect_MS_Code/Models/2017_Analyses/Daphnia_pulse/")
setwd("~/R") ##set to working directory
source("Daphnia_1.R")
    m.out <- rbind(
      pf.lik = sapply(mpar,function(x){
        if(class(x$pomp)=="mif2d.pomp") logLik(x$pomp) else NA
      }),
      sapply(mpar,function(x) {
        if(class(x$pomp)=="mif2d.pomp") coef(x$pomp) else rep(NA,length(coef(LV_restegg)))
      }),
      sapply(mpar,function(x)x$start)
    )
    save(m.out,mpar,file="out.rda")
toc <- Sys.time()
print(toc-tic)
print(m.out[1,])
#stopCluster(cl) ##turn on for flux

                                                         
load("out.rda")
max(as.vector(m.out[1,]), na.rm = TRUE)
p<-which.max(m.out[1,])

plist<-lapply(mpar,function(x)x[[1]])
run<-plist[[p]]
fitparamsnull<-run@params

pf<-replicate(10,pfilter(Daphnia_pomp,params=fitparamsnull,Np=1000))
ll<-sapply(pf,logLik)
mean(ll)
#[1] -210.2033
sd(ll)
#[1] 0.9810041
aic<-2*14-2*log(-mean(ll))
aic
#[1] 17.30385
max(as.vector(m.out[1,]), na.rm = TRUE)
#[1] -210.1886



                                                            

