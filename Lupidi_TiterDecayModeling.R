## Lupidi et al. (data modeling
# For questions, please contact Katharine Owers Bonner (owersk@gmail.com) & Peter Diggle (p.diggle@lancaster.ac.uk)

# Data source: Lupidi, R., et al. (1991). "Serological follow-up of patients involved in a 
# localized outbreak of leptospirosis." J Clin Microbiol 29(4): 805-809.

# We treat serovar as 3 levels initially, then switch to 2 contrasts when we add individual effect
# We also test for an initial rise in antibody titer, then remove it when it isn't significant


set.seed(8925)


setwd("~/LupidiModeling")
d=read.csv("LupidiDataFormatted.csv")
d$Patient=as.factor(d$Patient)

logR=log(c(1,2,6.4,20,64,200,640,2000))

################################################################################################
################################################################################################

## Model: Baseline ####

################################################################################################
################################################################################################


logL_bl<-function(k,mu,sigmasq) {
  
  pk=rep(NA,length(k))
 
   for (i in 1:length(k)){
    pk[i]=ifelse(k[i]==0, (pnorm(logR[1],mu,sigmasq)),
                 (pnorm(logR[k[i]+1],mu,sigmasq)-pnorm(logR[k[i]],mu,sigmasq)))
  }
  
  sum(log(pk))
}

my.fn_bl<-function(par) {
  mu<-par[1]; sigmasq<-par[2]
  result<-(-logL_bl(d$Dilution,mu,sigmasq))
  print(round(c(mu,sigmasq,result),3))
  result
}

initial.guess_bl<-c(2.11,1.99) # Sample Mean and Standard Deviation
result_bl<-optim(initial.guess_bl,my.fn_bl) 
result_bl 


################################################################################################
################################################################################################

## Model: Time  ####

################################################################################################
################################################################################################


logL_time<-function(k,Time,alpha,beta,sigmasq){
    
  pk=rep(NA,length(k))
  
  for (i in 1:length(k)){
    mu=alpha+beta*Time[i]
    pk[i]=ifelse(k[i]==0, (pnorm(logR[1],mu,sigmasq)),
                 (pnorm(logR[k[i]+1],mu,sigmasq)-pnorm(logR[k[i]],mu,sigmasq)))
  }
   
  sum(log(pk))
}

my.fn_time<-function(par) {
  alpha<-par[1]; beta<-par[2];sigmasq<-par[3]
  result<-(-logL_time(d$Dilution,d$Time,alpha,beta,sigmasq))
  print(round(c(alpha,beta,sigmasq,result),3))
  result
}

initial.guess_time<-c(result_bl$par[1],0,result_bl$par[2]) 
result_time<-optim(initial.guess_time,my.fn_time)  
result_time

# compare
1-pchisq(2*(-result_time$value+result_bl$value),1) 


################################################################################################
################################################################################################

## Model: Time & log(time) ####
# note:The log(time) term allows for an initial titer rise following infection

################################################################################################
################################################################################################


logL_timeLog<-function(k,Time,alpha,beta,eta,sigmasq){
  pk=rep(NA,length(k))
  
  for (i in 1:length(k)){
    mu=alpha+beta*Time[i]+eta*log(Time[i])
    pk[i]=ifelse(k[i]==0, (pnorm(logR[1],mu,sigmasq)),
                 (pnorm(logR[k[i]+1],mu,sigmasq)-pnorm(logR[k[i]],mu,sigmasq)))
  }
  
  sum(log(pk))
}

my.fn_timeLog<-function(par) {
  alpha<-par[1]; beta<-par[2];eta<-par[3];sigmasq<-par[4]
  result<-(-logL_timeLog(d$Dilution,d$Time,alpha,beta,eta,sigmasq))
  print(round(c(alpha,beta,eta,sigmasq,result),3))
  result
}

initial.guess_timeLog<-c(result_time$par[1],result_time$par[2],0,result_time$par[3]) 
result_timeLog<-optim(initial.guess_timeLog,my.fn_timeLog) 
result_timeLog

# compare
1-pchisq(2*(-result_timeLog$value+result_time$value),1) # This is NOT a significant improvement, so we remove it



################################################################################################
################################################################################################

## Model: Time & Serovar  ####

#NOTE: Serovar has 3 levels here

################################################################################################
################################################################################################



logL_timeSV3<-function(k,Time,SV,alphaAUS,alphaBRAT,alphaLORA,beta,sigmasq) {
  
  pk=rep(NA,length(k))
  for (i in 1:length(k)){
    alpha=ifelse(SV[i]=="australis",alphaAUS,ifelse(SV[i]=="bratislava",alphaBRAT,alphaLORA))
    mu=alpha+beta*Time[i]
    pk[i]=ifelse(k[i]==0, (pnorm(logR[1],mu,sigmasq)),
                 (pnorm(logR[k[i]+1],mu,sigmasq)-pnorm(logR[k[i]],mu,sigmasq)))
  }
  sum(log(pk))
}

my.fn_timeSV3<-function(par) {
  alphaAUS<-par[1];alphaBRAT<-par[2];alphaLORA<-par[3];beta<-par[4];sigmasq<-par[5]
  result<-(-logL_timeSV3(d$Dilution,d$Time,d$SV,alphaAUS,alphaBRAT,alphaLORA,beta,sigmasq))
  print(round(c(alphaAUS,alphaBRAT,alphaLORA,beta,sigmasq,result),3))
  result
}

initial.guess_timeSV3<-c(result_time$par[1],result_time$par[1],result_time$par[1],
                         result_time$par[2], result_time$par[3])
result_timeSV3<-optim(initial.guess_timeSV3,my.fn_timeSV3) 
result_timeSV3 

# compare
1-pchisq(2*(-result_timeSV3$value+result_time$value),2) 

##########################################################################################################
##########################################################################################################
##########################################################################################################

## FINAL MODEL: Time & SV & 18 individual coefficients ##

# NOTE: Re-parameterized serovar to be the contrasts between serovars (2 levels instead of 3) so each of
    # the 18 individuals has an intercept (instead of being a contast with a "baseline" person)
    # The "baseline" concept makes more sense for a serovar (Australis in this case)

##########################################################################################################
##########################################################################################################
##########################################################################################################

logL_TimeSV18<-function(k,Patient,Time,SV,alphaB,alphaL,beta,sigmasq,
          P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18) {
  
  pk=rep(NA,length(k))
  
  for (i in 1:length(k)){
    
    alpha=ifelse(SV[i]=="australis",0,ifelse(SV[i]=="bratislava",alphaB,alphaL))
    
    PatientCoef=ifelse(Patient[i]==1,P1,ifelse(Patient[i]==2,P2,ifelse(Patient[i]==3,P3,ifelse(Patient[i]==4,P4,
           ifelse(Patient[i]==5,P5,ifelse(Patient[i]==6,P6,ifelse(Patient[i]==7,P7,ifelse(Patient[i]==8,P8,
           ifelse(Patient[i]==9,P9,ifelse(Patient[i]==10,P10,ifelse(Patient[i]==11,P11,ifelse(Patient[i]==12,P12,
           ifelse(Patient[i]==13,P13,ifelse(Patient[i]==14,P14,ifelse(Patient[i]==15,P15,ifelse(Patient[i]==16,P16,
           ifelse(Patient[i]==17,P17,P18)))))))))))))))))
   
    mu=alpha+beta*Time[i]+PatientCoef
    
    pk[i]=ifelse(k[i]==0, (pnorm(logR[1],mu,sigmasq)),
                 (pnorm(logR[k[i]+1],mu,sigmasq)-pnorm(logR[k[i]],mu,sigmasq)))
  }
  
  sum(log(pk))
}

my.fn_TimeSV18<-function(par) {
  alphaB<-par[1];alphaL<-par[2];beta<-par[3];sigmasq<-par[4];
  P1<-par[5];P2<-par[6];P3<-par[7];P4<-par[8];P5<-par[9];P6<-par[10];P7<-par[11];P8<-par[12];P9<-par[13];
  P10<-par[14];P11<-par[15];P12<-par[16];P13<-par[17];P14<-par[18];P15<-par[19];P16<-par[20];
  P17<-par[21];P18<-par[22];
  
  result<-(-logL_TimeSV18(d$Dilution,d$Patient,d$Time,d$SV,alphaB,alphaL,beta,sigmasq,
                     P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18))
  
  print(round(c(alphaB,alphaL,beta,sigmasq,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,
                P14,P15,P16,P17,P18,result),3))
  result
}


initial.guess_TimeSV18<-c((result_timeSV3$par[2]-result_timeSV3$par[1]), # alphaB
                          (result_timeSV3$par[3]-result_timeSV3$par[1]), # alphaL
                          result_timeSV3$par[4],result_timeSV3$par[5], # beta, sigmasq
                          rep(result_timeSV3$par[1],18)) #18FE

result_TimeSV18<-optim(initial.guess_TimeSV18,my.fn_TimeSV18,
                         control=list(maxit=100000),hessian=T)
result_TimeSV18

# Compare
1-pchisq(2*(-result_TimeSV18$value+result_timeSV3$value),17)


##########################################################################################################
##########################################################################################################
##########################################################################################################

##########################################################################################################
##########################################################################################################

# Explore & Save the output 

##########################################################################################################
##########################################################################################################

#### Exploring the Hessian matrix
H<-result_TimeSV18$hessian
iH<-solve(H) # The variance matrix

sm<-iH[1:4,1:4] # Only the interesting parameters (2 alpha contrasts, beta, sigmasq)
round(sm,8) 

setwd("~/Results")
write.csv(H, "2018_03_08-HessianMatrix_Time-Serovar-18IndivCoefs.csv",row.names=F)

#### Correlation matrix
rmat<-sm/sqrt(outer(diag(sm),diag(sm),"*"))  
rmat 


## Appendix Table 1 
names=c("Baseline","Time","Time +log(Time)", "Time + Serovar", "Time + Serovar + Individual")
LogL=c(result_bl$value,result_time$value,result_timeLog$value,result_timeSV3$value,result_TimeSV18$value)
Deviance=c(NA, 
           2*(result_time$value-result_bl$value),
           2*(result_timeLog$value-result_time$value),
           2*(result_timeSV3$value-result_time$value),
           2*(result_TimeSV18$value-result_timeSV3$value))
df=c(NA,1,1,2,17)
sig=c(NA, 
          1-pchisq(2*(-result_time$value+result_bl$value),1),
          1-pchisq(2*(-result_timeLog$value+result_time$value),1),
          1-pchisq(2*(-result_timeSV3$value+result_time$value),2),
          1-pchisq(2*(-result_TimeSV18$value+result_timeSV3$value),17)
          )
t1=cbind(names,LogL,Deviance,df,sig)
#write.csv(t1,"2018_03_08_AppendixTable1_DevianceTable.csv",row.names=F)

#### Appendix Table 2 (And general info)
SEmod=sqrt(diag(iH))
Variable<-c("alphaB","alphaL","beta","sigmasq","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12",
            "P13","P14","P15","P16","P17","P18")
CSE=cbind(Variable,result_TimeSV18$par,SEmod)
colnames(CSE)=c("Parameter","Estimate","Standard Error")
#write.csv(CSE,paste0(Sys.Date(),"-TimeSV18-Coefs&StErrs.csv"),row.names=F)


##########################################################################################################
##########################################################################################################

# Calculate and plot Standardized Residuals

##########################################################################################################
##########################################################################################################

n=length(d[,1])
PlotTitle="Standardized Residuals, Time + serovar + Individual"

ResidualPlot_TimeSV18<-function(Dilution,Patient,Time,SV,
                                alphaB,alphaL,beta,sigmasq,
                                P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18){ 
  pk=matrix(NA,nrow=n,ncol=8)
  mu_hat=rep(NA,n)
  sigmasq_hat=rep(NA,n)
  
  for (i in 1:n){
    alpha=ifelse(SV[i]=="australis",0,ifelse(SV[i]=="bratislava",alphaB,alphaL))
    PatientCoef=ifelse(Patient[i]==1,P1,ifelse(Patient[i]==2,P2,ifelse(Patient[i]==3,P3,ifelse(Patient[i]==4,P4,
       ifelse(Patient[i]==5,P5,ifelse(Patient[i]==6,P6,ifelse(Patient[i]==7,P7,ifelse(Patient[i]==8,P8,
       ifelse(Patient[i]==9,P9,ifelse(Patient[i]==10,P10,ifelse(Patient[i]==11,P11,ifelse(Patient[i]==12,P12,
       ifelse(Patient[i]==13,P13,ifelse(Patient[i]==14,P14,ifelse(Patient[i]==15,P15,ifelse(Patient[i]==16,P16,
       ifelse(Patient[i]==17,P17,P18)))))))))))))))))
    mu=alpha+beta*Time[i]+PatientCoef
    
    for (bin in 0:7){ 
      pk[i,bin+1]=ifelse(bin==0, (pnorm(logR[1],mu,sigmasq)),
                         (pnorm(logR[bin+1],mu,sigmasq)-pnorm(logR[bin],mu,sigmasq)))
    }}
  
  for (point in 1:n){
    mu_hat[point]=pk[point,1]*0+pk[point,2]*1+pk[point,3]*2+pk[point,4]*3+pk[point,5]*4+pk[point,6]*5+
      pk[point,7]*6+pk[point,8]*7

    sigmasq_hat[point]=pk[point,1]*0+pk[point,2]*1+pk[point,3]*2^2+pk[point,4]*3^2+pk[point,5]*4^2+
      pk[point,6]*5^2+pk[point,7]*6^2+pk[point,8]*7^2-mu_hat[point]^2
  }

  pk<-cbind(pk,mu_hat,Dilution,sigmasq_hat,((Dilution-mu_hat)/sqrt(sigmasq_hat)))
  plot(jitter(mu_hat, factor=0.6),jitter(pk[,12],factor=75),
       cex=0.5,xlab="Fitted Values",ylab="Standardized Residuals") #main=PlotTitle,
  return(pk)
}


TimeSV18<-ResidualPlot_TimeSV18(d$Dilution,d$Patient,d$Time,d$SV,
                                result_TimeSV18$par[1],result_TimeSV18$par[2],
                                result_TimeSV18$par[3],result_TimeSV18$par[4],result_TimeSV18$par[5],
                                result_TimeSV18$par[6],result_TimeSV18$par[7],result_TimeSV18$par[8],
                                result_TimeSV18$par[9],result_TimeSV18$par[10],result_TimeSV18$par[11],
                                result_TimeSV18$par[12],result_TimeSV18$par[13],result_TimeSV18$par[14],
                                result_TimeSV18$par[15],result_TimeSV18$par[16],result_TimeSV18$par[17],
                                result_TimeSV18$par[18],result_TimeSV18$par[19],result_TimeSV18$par[20],
                                result_TimeSV18$par[21],result_TimeSV18$par[22]) 


colnames(TimeSV18)=c("Lik0","Lik1","Lik2","Lik3","Lik4","Lik5","Lik6","Lik7","mu_hat",
                     "K","sigmasq_hat","StdResid")
res=cbind(d[,c(1:3,7,10)],TimeSV18)
write.csv(res,"2018_03_08-StdResiduals_Time+SV+18FE.csv",row.names=F)

# Appendix Figure
png("2018_03_08-AppFig_Residuals vs fitted and by indiv.png",width=960)
par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(res$mu_hat,res$StdResid,xlab="Fitted Values",ylab="Standardized Residual")

plot(jitter(as.numeric(res$Patient)),res$StdResid,col=as.factor(res$SV),pch=19,xlab="Individual",
     ylab="Standardized Residual")
dev.off()
