## Example script for applying titer decay rate to cohort data
# August 2020
# For questions, please contact Katharine Owers Bonner (owersk@gmail.com) & Peter Diggle (p.diggle@lancaster.ac.uk)


# For this example, we assume two-fold dilutions (1:50, 1:100, 1:200...)


#############################################################################
# Import estimated titer decay rate (R)

# use for this example:
decayrate <- 0.926

#############################################################################
# in a long data.frame (here d), with one row per interval per individual, include columns for 
  # ID
  # StartDilution (K1) (in the dilution format, not titer, i.e. a titer of 1:50 is 1)
  # EndDilution (K2) (in the dilution format, not titer, i.e. a titer of 1:50 is 1)
  # Interval between samples (in months)



ID<- rep(1:5,each=2)
StartDilution <- c(1,4,0,0,2,3,0,5,1,1)
EndDilution <-   c(1,3,1,0,1,4,0,3,0,2)
Interval <-rep(6,10)
 
d<-data.frame(ID,StartDilution,EndDilution,Interval)

# Apply decay to data

# Generate empty vectors to fill with imputed & decayed values
StartImputed=NULL
EndDecayed=NULL
EndRounded=NULL
Infect_Decayed=NULL
Infect_Decayed_Type=NULL

for (i in 1:length(d$ID)){
  
  ## Impute a start value (within its dilution)
  StartImputed[i]=d$StartDilution[i]+runif(1)
  
  ## Apply decay rate to the imputed titer at the start of the time interval
  EndDecayed[i]= StartImputed[i]*decayrate^d$Interval[i] 
  
  ## Assign the decayed titer a dilution to match what happens with MAT, where you censor it 
      # into a bin based on the threshold of 50% agglutination
  EndRounded[i]=floor(EndDecayed[i])
  
  ## Define Infections
  Infect_Decayed[i]=ifelse(EndRounded[i]==0 & d$EndDilution[i] >=1,1, # seroconversion
                           ifelse(EndRounded[i]>=1 & d$EndDilution[i]>=(2+EndRounded[i]),1,0)) #four-fold rise
  
  ## Classify infections
  Infect_Decayed_Type[i]=ifelse(Infect_Decayed[i]==0,NA,
                                ifelse(EndRounded[i]==0,"SC","FFT"))
} 

d<-cbind(d,StartImputed,EndDecayed,EndRounded, Infect_Decayed,Infect_Decayed_Type)  
