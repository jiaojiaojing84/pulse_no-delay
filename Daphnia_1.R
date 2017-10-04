#Daphnia model with pulse each year that occurs one week prior to first observation of Daphnia
#setwd("~/Google_Drive/Bytheffect_MS_Code/Models/2017_Analyses/Daphnia_1wk_pulse/")
require(pomp); require(lubridate)
Daphnia_rprocess <- Csnippet("
                            double births;
                            double deaths;
                            double betaz;
                            double dw;

                            if (error_count != 0.0) return;
                            
                            betaz = exp(dot_product(4,&seas_1,&log_betaz1));

                            dw = rnorm(0,sqrt(dt));  // white noise for prey birth

                            births = (betaz*(1-eta*byth)*dt)*prey*(1-prey/cc) + prey*sdbeta*dw + pulse*rlnorm(meanp,sigmap);	// prey births
                            deaths = prey*(alpha*byth + mu)*dt;	// prey deaths due to predation

                            prey += (births - deaths)*nocoup;
                            noise += dw;

                            // check for violations of positivity constraints
                            // nonzero error_count variable signals violation
                            if (prey <= 0.0)  {
                            prey = 0;
                            }
                            ")

params.init <- c(taudem=0.19809275,tauenv=0.34300230,
                 log.betaz1=-7.95806589,log.betaz2=-2.77616985,log.betaz3=-1.40980823,log.betaz4=-0.77433155,
                 alpha=0,mu=0.05921878,eta=0.32815929,meanp=-3.05757131,sigmap=1.73135194,
                 cc=35.46637919,sdbeta=0.26921510,
                 prey.0=0.05281949)

zoop<-read.table("GLERL_M110_Zoop-Chla-Temp-Secchi_1994-2012.txt",header=TRUE)
zoop$day<-as.numeric(1+(as.Date(zoop$Date, format="%m/%d/%y")-rep(as.Date("1994-01-01"),length(zoop$Date)))) 
zoop$cases<-complete.cases(zoop$D.mendotae);zoopcomplete<-subset(zoop,cases==TRUE)

require(timeSeries)
Byth<-rep(c(rep(0,100),rep(NA,265)),19); 
Byth[zoopcomplete$day]<-zoopcomplete$Bythotrephes; Byth[1]<-0;Byth[6935]<-0
Byth<-interpNA(Byth, method = "linear")
#plot(Byth,type="l")
nocoup<-rep(c(0,rep(1,364)),19)


#figure out first date of daphnia observation each year
#need to determine first observation of Daphnia each year
zoopcomplete$year<-year(as.Date(zoopcomplete$Date, format="%m/%d/%y"))
surveyyears<-c(1994:2003,2007:2012)
firstdaph<-rep(0,16)
firstdaphjday<-rep(0,16)

firstdaphdens<-rep(0,16)

for (i in 1:16) {
  zoopyear<-subset(zoopcomplete,year==surveyyears[i])
  zoopyearno0<-subset(zoopyear,D.mendotae>0)
  firstdaph[i]<-zoopyearno0$day[which.min(zoopyearno0$JulianDay)]
  firstdaphjday[i]<-zoopyearno0$JulianDay[which.min(zoopyearno0$JulianDay)]
  firstdaphdens[i]<-zoopyearno0$D.mendotae[which.min(zoopyearno0$JulianDay)]
}

daphstart<-firstdaph

#create vector with 1s at time of first observation
pulse<-rep(0,6935)
pulse[daphstart]<-1

#need to generate subset of Daphnia observations excluding early zero observations

#first, need ot determine how many zero observations each year
#initzerocount<-rep(0,16)

#for (i in 1:16) {
 #zoopyear<-subset(zoopcomplete,year==surveyyears[i])
 #zoopyearpreobs<-subset(zoopyear,JulianDay<firstdaphjday[i])
 #initzerocount[i]<-length(zoopyearpreobs$D.mendotae)
#}
#sum(initzerocount) #41 observations
#zerodays<-rep(0,41)

zoopsubset<-read.table("GLERL_M110_Zoop-Chla-Temp-Secchi_1994-2012_init0sremoved.txt",header=TRUE)
zoopsubset$day<-as.numeric(1+(as.Date(zoopsubset$Date, format="%m/%d/%y")-rep(as.Date("1994-01-01"),length(zoopsubset$Date)))) 


#for(i in 1:16) {
#index<-daphstart[i] + 365*(i-1)
#pulse[index]<-1  
#}

covar.dt <- 1
t0 <- 0
nbasis <- 4
tcovar <- seq(from=t0,to=6934,by=covar.dt)
yr <- 1:365


covartable <- data.frame(
  time=tcovar,
  seas=periodic.bspline.basis(tcovar,nbasis=nbasis,degree=3,period=365),
  byth=Byth,
  nocoup,
  pulse
)

Daphnia_rmeasure <- Csnippet("
                                  double tau, daphnia;
                                  tau=sqrt(taudem*taudem*prey + tauenv*tauenv*prey*prey);
                                  if ((error_count > 0.0) || (!(R_FINITE(prey)))) {
                                  D_mendotae = R_NaReal;
                                  } else {
                                  daphnia = rnorm(prey,tau);
                                  if (daphnia<=0) {
                                  D_mendotae=0;
                                  }
                                  else {
                                  D_mendotae=daphnia;
                                  }
                                  }
                                  ")

Daphnia_dmeasure <- Csnippet("
                                  double tau, tol = 1.0e-18;
                                  tau=sqrt(taudem*taudem*prey + tauenv*tauenv*prey*prey);
                                  double f = 0.0;
                                  if ((error_count > 0.0) || (!(R_FINITE(prey)))) {
                                  lik = tol;
                                  } else {
                                  if (D_mendotae==0) {
                                  f += pnorm(0,prey,tau,1,1)+tol;
                                  }
                                  else {
                                  f += dnorm(D_mendotae, prey, tau, 1)+tol;
                                  }
                                  lik = (give_log) ? f : exp(f);
                                  }
                                  ")


Daphnia_untrans <- Csnippet("
                       Ttaudem = log(taudem);
                       Ttauenv = log(tauenv);
                       Tcc = log(cc);
                       Tsdbeta = log(sdbeta);
                       Tmu = log(mu);
                       Teta = log(eta);
                       Tsigmap = log(sigmap);
                       ")

Daphnia_trans <- Csnippet("
                     Ttaudem = exp(taudem);
                     Ttauenv = exp(tauenv);
                     Tcc = exp(cc);
                     Tsdbeta = exp(sdbeta);
                     Tmu = exp(mu);
                     Teta = exp(eta);
                     Tsigmap = exp(sigmap);
                     ")

Daphnia_pomp <- pomp(
  data=subset(zoopsubset,select=c("D.mendotae","day")),
  times="day",
  t0=0,
  params=params.init,
  rprocess = euler.sim(step.fun = Daphnia_rprocess, delta.t=1),
  rmeasure= Daphnia_rmeasure,
  dmeasure = Daphnia_dmeasure,
  covar=covartable,
  tcovar="time",
  obsnames = c("D.mendotae"),
  zeronames = c("error_count"),
  statenames = c("prey","noise","error_count"),
  paramnames = c("taudem","tauenv","log.betaz1","cc","alpha","mu","eta","sdbeta","meanp","sigmap","prey.0"),
  covarnames = c("seas.1","byth","nocoup","pulse"),
  all.state.names=c("prey","noise", "error_count"),
  comp.names=c("prey"),
  comp.ic.names=c("prey.0"),
  fromEstimationScale=Daphnia_trans,
  toEstimationScale=Daphnia_untrans,
  initializer = function (params, t0, comp.ic.names, comp.names, all.state.names, ...) {
    states <- numeric(length(all.state.names))
    names(states) <- all.state.names
    frac <- params[comp.ic.names]
    states[comp.names] <- frac
    states
  }
)

#Daphnia_sim<-simulate(Daphnia_pomp)
#plot(Daphnia_sim)
pf<-pfilter(Daphnia_pomp,params=params.init,Np=2000)
pf@loglik
