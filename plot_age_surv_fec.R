library(readxl)
wd="C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT"

## Set working directory and initialize  ----
setwd("C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT/Data")
file.remove("cogediv.txt",showWarnings = FALSE)
file.remove("cogediv_sensibility.txt",showWarnings = FALSE)
allinfo=0
lifetime<-as.data.frame(read_excel("GENETIC_DIVERSITY_DATA.xlsx",sheet="agene"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)

## 
lifetime<-as.data.frame(read_excel("GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
line=c(line,nrow(lifetime)+1)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)

species_c_f=data.frame(Sp=rep(NA,16),
                       c=rep(0,16),
                       f_l=rep(NA,16),
                       beta=rep(0,16)
)

for (i in 1:16){
  lifetime_tmp=lifetime[line[i]:(line[i+1]-1),]
  
  #Surv
  Surv_c=exp(-((((lifetime_tmp$Length)/(lifetime$Linf[line[i]]))^(-1.5))*lifetime$K[line[i]]))
  Surv_f=exp(-((((lifetime_tmp$Length_F)/(lifetime$Linf_F[line[i]]))^(-1.5))*lifetime$K_F[line[i]]))
  Surv_m=exp(-((((lifetime_tmp$Length_M)/(lifetime$Linf_M[line[i]]))^(-1.5))*lifetime$K_M[line[i]]))
  
  #Cum surv
  Cum_surv_c=c(NA,1)
  Cum_surv_c=c(1)
  for (j in 2:length(Surv_c)){
    Cum_surv_c[j]=Cum_surv_c[j-1]*Surv_c[j-1]
  }
  
  Cum_surv_f=c(NA,1)
  Cum_surv_f=c(1)
  for (j in 2:length(Surv_f)){
    Cum_surv_f[j]=Cum_surv_f[j-1]*Surv_f[j-1]
  }
  
  Cum_surv_m=c(NA,1)
  Cum_surv_m=c(1)
  for (j in 2:length(Surv_m)){
    Cum_surv_m[j]=Cum_surv_m[j-1]*Surv_m[j-1]
  }

  # Fecundity
  # Combined
  if (lifetime$Fecundity_length_relationship[line[i]]=="Linear"){
    fec_inter=lifetime$Alpha[line[i]]+(lifetime_tmp$Length*lifetime$Beta[line[i]]) 
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  } else if (lifetime$Fecundity_length_relationship[line[i]]=="Power"){
    fec_inter=lifetime$Alpha[line[i]]*(lifetime_tmp$Length^(lifetime$Beta[line[i]]))   
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  } else if (lifetime$Fecundity_length_relationship[line[i]]=="Exponential"){
    fec_inter=lifetime$Alpha[line[i]]*exp(lifetime_tmp$Length*lifetime$Beta[line[i]])  
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  }
  
  fec_inter[1]=0
  if (lifetime$Maturity[line[i]]>1){
    for (j in 1:lifetime$Maturity[line[i]]){
      fec_inter[j+1]=0
    }
  }
  
  fec_inter_c=fec_inter
  
  # Female
  if (lifetime$Fecundity_length_relationship[line[i]]=="Linear") {
    fec_inter=lifetime$Alpha[line[i]]+(lifetime_tmp$Length_F*lifetime$Beta[line[i]]) 
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  } else if (lifetime$Fecundity_length_relationship[line[i]]=="Power") {
    fec_inter=lifetime$Alpha[line[i]]*(lifetime_tmp$Length_F^(lifetime$Beta[line[i]]))   
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  } else if (lifetime$Fecundity_length_relationship[line[i]]=="Exponential") {
    fec_inter=lifetime$Alpha[line[i]]*exp(lifetime_tmp$Length_F*lifetime$Beta[line[i]])  
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  }
  fec_inter[1]=0
  if (lifetime$Maturity_F[line[i]]>1){
    for (j in 1:lifetime$Maturity_F[line[i]]){
      fec_inter[j+1]=0
    }
  }
  
  fec_inter_f=fec_inter
  
  # Male
  if (lifetime$Fecundity_length_relationship[line[i]]=="Linear") {
    fec_inter=lifetime$Alpha[line[i]]+(lifetime_tmp$Length_M*lifetime$Beta[line[i]]) 
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  } else if (lifetime$Fecundity_length_relationship[line[i]]=="Power") {
    fec_inter=lifetime$Alpha[line[i]]*(lifetime_tmp$Length_M^(lifetime$Beta[line[i]]))   
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  } else if (lifetime$Fecundity_length_relationship[line[i]]=="Exponential") {
    fec_inter=lifetime$Alpha[line[i]]*exp(lifetime_tmp$Length_M*lifetime$Beta[line[i]])  
    fec_inter=fec_inter/(max(fec_inter,na.rm=T))
  }
  fec_inter[1]=0
  if (lifetime$Maturity_M[line[i]]>1){
    for (j in 1:lifetime$Maturity_M[line[i]]){
      fec_inter[j+1]=0
    }
  }
  
  fec_inter_m=fec_inter
  
  age_specific=data.frame(Age=seq(0,nrow(lifetime_tmp)-1),
                          Surv=Surv_c,
                          Surv_f=Surv_f,
                          Surv_m=Surv_m,
                          Cum_surv_c=Cum_surv_c,
                          Cum_surv_f=Cum_surv_f,
                          Cum_surv_m=Cum_surv_m,
                          fec_c=fec_inter_c,
                          fec_f=fec_inter_f,
                          fec_m=fec_inter_m
  )
  
  # Fit surv
  #ss <- nls(Cum_surv_c ~ exp(-((Age/b)^(c))),
            ss <- nls(Surv ~ 1-(1- exp(((Age/b)^(c))-(((Age+1)/b)^(c)))),
            start = list(b=1,c = 1),
            lower=c(0.005,0.005), upper=c(100,20),
            control=nls.control(maxiter=1000),
            algorithm="port",
            data=age_specific)
  age_specific$est=c(fitted(ss))
  summary(ss)
  
  #plot(age_specific$Age,age_specific$Surv,
  #☻     type="l",
  ###     main="Fit age-specific survival",
  #     xlab="Age",
  #     ylab="Cum. surv.",
  #     ylim=c(0,1))
  
  #lines(age_specific$Age,age_specific$est,col="blue",lty=2)
  #readline("")
  species_c_f[i,]=c(lifetime_tmp$Species_code[1],
                    coefficients(ss)[2],
                    lifetime$Fecundity_length_relationship[line[i]],
                    lifetime$Beta[line[i]])
  
  # Plot
  
  p<-ggplot(age_specific,aes(x=Age,y=Cum_surv_c))+
    geom_line(aes(y=Cum_surv_f),size=2,col="red",alpha=0.5)+
    geom_line(aes(y=Cum_surv_m),size=2,col="blue",alpha=0.5)+
    geom_line(size=2,alpha=0.75)+
    geom_line(aes(y=fec_m),linetype='dashed',size=2,col="red",alpha=0.5)+
    geom_line(aes(y=fec_f),linetype='dashed',size=2,col="blue",alpha=0.5)+
    geom_line(aes(y=fec_c),linetype='dashed',size=2,alpha=0.75)+
    theme_classic()+
    scale_y_continuous(
      "Probability of survival until age", 
      sec.axis = sec_axis(~ . * 1, name = "Age-specific fecundity")
    )+
    #ylim(c(0,1))+
    labs(#title = "Montagu's blenny",
           title = lifetime_tmp$Latin[1])+
    scale_x_continuous("Age",breaks = seq(0,nrow(age_specific)-1))+
    theme(
      plot.title = element_text(hjust = 0.5,face='italic')
      #plot.subtitle = element_text(hjust = 0.5,face='italic')
    )
  
  pdf(paste(wd,"/figures/Surv_fec_",lifetime_tmp$Species_code[1],".pdf",sep=""),width=5.5,height=3.33)
  print(p)
  dev.off()
}

round(age_specific[-1,],2)
round(age_specific[-1,],3)

save(species_c_f,file="species_c_f.Rdata")
