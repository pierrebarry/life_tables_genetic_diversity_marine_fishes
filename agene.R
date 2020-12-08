################################
### This script produces AgeNe 
###############################

## Load packages ----
library(betareg)
library(INLA)
library(readxl)
library(readtext)
library(stringr)
library(ggthemes)

## Load functions ----
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 
Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

## Set working directory and initialize  ----
setwd("C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT/Data/agene")
file.remove("cogediv.txt",showWarnings = FALSE)
file.remove("cogediv_sensibility.txt",showWarnings = FALSE)
allinfo=0
lifetime<-as.data.frame(read_excel("agene.xlsx",sheet="RawData"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)

## Construct and run agene ----

set_model=matrix(nrow=32,ncol=5,rep(0,5*32))
set_model_first_column=c(rep(0,16),rep(1,16))
set_model_second_column=rep(c(rep(0,8),rep(1,8)),2)
set_model_third_column=rep(c(rep(0,4),rep(1,4)),4)
set_model_four_column=rep(c(rep(0,2),rep(1,2)),8)  
set_model_five_column=rep(c(rep(0,1),rep(1,1)),16)   
set_model=as.data.frame(cbind(set_model_five_column,
                              set_model_four_column,
                              set_model_third_column,
                              set_model_second_column,
                              set_model_first_column))
colnames(set_model)=c("AgeMat","Surv","Fec","Sex","Hermaphro")

n_output=length(species)
agene_output=vector('list',5)


for (ele in 1:5){
  agene_output[[ele]]=data.frame(Species_code=rep(NA,n_output),
                                 Latin=rep(NA,n_output),
                                 Vernacular=rep(NA,n_output)
  )
  agene_output[[ele]]$Species_code=species
  agene_output[[ele]]$Latin=latin
  agene_output[[ele]]$Vernacular=vernacular
  
  for (j in 1:32){
    agene_output[[ele]]=cbind(agene_output[[ele]],rep(0,16))
  }
  
  colnames(agene_output[[ele]])=c("Species_code","Latin","Vernacular",paste("Output",seq(1,32),sep=""))
}


for (ll in 1:16){
  
  lifetime<-as.data.frame(read_excel("agene.xlsx",sheet="RawData"))
  
  # Compute loop
  species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
  latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
  vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
  
  line=which(is.na(lifetime$Species_code)==FALSE)
  lifetime$Length_F=as.numeric(lifetime$Length_F)
  lifetime$Length_M=as.numeric(lifetime$Length_M)
  lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
  lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
  
  #print(set_model[ll,])
  file.remove("cogediv_sensibility.txt",showWarnings=FALSE)
  # 0 = no difference
  diff_ageatmaturity=set_model$AgeMat[ll] #OK
  diff_age_specific_survival=set_model$Surv[ll] #OK
  diff_age_specific_fecundity=set_model$Fec[ll] # OK
  diff_sex_specific_vitalrate=set_model$Sex[ll] # OK set to female
  diff_hermaphro=set_model$Hermaphro[ll] # OK
  
  # CONSTRUCT DATA FOR AGENE ------
  
  # Set hermaphrodism Y/N
  if (diff_hermaphro==0){
    for (i in 1:length(species)){
      lifetime$AgeNe_model[line[i]]="TwoSex"
    }
  }
  
  # Set Age at Maturity Y/N
  if (diff_ageatmaturity==0){
    for (i in 1:length(species)){
      lifetime$Maturity_F[line[i]]=1
      lifetime$Maturity_M[line[i]]=1
      lifetime$Maturity[line[i]]=1
    }
  }
  
  # Set Sex-specific vital rates Y/N
  if (diff_sex_specific_vitalrate==0){
    lifetime$Max_age_F=lifetime$Max_age
    lifetime$Max_age_M=lifetime$Max_age
    
    lifetime$Length_F=lifetime$Length
    lifetime$Length_M=lifetime$Length
    
    lifetime$Linf_F=lifetime$Linf
    lifetime$Linf_M=lifetime$Linf
    
    lifetime$K_F=lifetime$K
    lifetime$K_M=lifetime$K
    
    lifetime$t0_F=lifetime$t0
    lifetime$t0_M=lifetime$t0
    
    lifetime$Maturity_F=lifetime$Maturity
    lifetime$Maturity_M=lifetime$Maturity

  }
  
  surviv_plot=vector('list',length(species))
  fecun_plot=vector('list',length(species))
  for (i in 1:length(species)){
    surviv_plot[[i]]=c(surviv_plot[[i]],1)
  }
  for (i in 1:length(species)){
    
    
    #First line
    firstline=paste(lifetime$Latin[line[i]]," - Length ref : ",
                    lifetime$Length_ref[line[i]]," - Fecundity ref : ",
                    lifetime$Fecundity_ref[line[i]])
    if (file.exists("cogediv_sensibility.txt")==FALSE){
      write(firstline,file="cogediv_sensibility.txt",append="FALSE")
    } else {
      write(firstline,file="cogediv_sensibility.txt",append="TRUE")
    }
    
    
    # Second line
    if (lifetime$AgeNe_model[line[i]]=="TwoSex"){
      secondline=c(lifetime$Max_age[line[i]],
                   "1000",
                   lifetime$Initial_SexRatio[line[i]])
      write(secondline,file="cogediv_sensibility.txt",append="TRUE",ncolumns = 3)
    } else {
      secondline=c(lifetime$Max_age[line[i]],
                   "1000",
                   lifetime$First_sexe[line[i]])
      write(secondline,file="cogediv_sensibility.txt",append="TRUE",ncolumns = 3)
    }
    
    
    # Other line
    if (lifetime$AgeNe_model[line[i]]=="TwoSex" | lifetime$AgeNe_model[line[i]]=="SexChange"){
      for (j in 1:lifetime$Max_age[line[i]]){
        
        if (diff_age_specific_survival==1){
          if (j<lifetime$Max_age_F[line[i]]){
            female_survival=exp(-((((lifetime$Length_F[line[i]+j-1])/(lifetime$Linf_F[line[i]]))^(-1.5))*lifetime$K_F[line[i]]))
          } else {
            female_survival=0
          }
        } else {
          const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[i]]-1))
          if (j<lifetime$Max_age_F[line[i]] & j>1){
            female_survival=const_surv
          } else if (j>1){
            female_survival=0
          } else {
            female_survival=1
          }
        }
        
        if (is.na(female_survival)==TRUE){
          female_survival=0
        }
        
        if (diff_age_specific_fecundity==1){
          if (j>=lifetime$Maturity_F[line[i]]){
            if (lifetime$Fecundity_length_relationship[line[i]]=="Power"){
              female_fecundity=lifetime$Alpha[line[i]]*(lifetime$Length_F[line[i]+j-1]^(lifetime$Beta[line[i]]))      
            } else if (lifetime$Fecundity_length_relationship[line[i]]=="Exponential") {
              female_fecundity=lifetime$Alpha[line[i]]*exp(lifetime$Length_F[line[i]+j-1]*lifetime$Beta[line[i]])      
            } else if (lifetime$Fecundity_length_relationship[line[i]]=="Linear") {
              female_fecundity=lifetime$Alpha[line[i]]+(lifetime$Length_F[line[i]+j-1]*lifetime$Beta[line[i]])      
            }      
          } else {
            female_fecundity=0
          }
        } else {
          if (j>=lifetime$Maturity_F[line[i]]){
            female_fecundity=1
          } else {
            female_fecundity=0
          }
          
        }
        

        
        if (is.na(female_fecundity)==FALSE & female_fecundity==0){
          female_fecundity=0.01
        }
        
        if (is.na(female_fecundity)==TRUE){
          female_fecundity=0
        }
        
        
        PF_Female=lifetime$PF_F[line[i]]
        
        if (diff_age_specific_survival==1){
          if (j<as.numeric(lifetime$Max_age_M[line[i]])){
            male_survival=exp(-((((lifetime$Length_M[line[i]+j-1])/(lifetime$Linf_M[line[i]]))^(-1.5))*lifetime$K_M[line[i]]))
          } else {
            male_survival=0
          }
        } else {
          const_surv=exp((log(0.01))/(as.numeric(lifetime$Max_age_M[line[i]])-1))
          if (j==1){
            male_survival=1
          } else if (j<as.numeric(lifetime$Max_age_M[line[i]])){
            male_survival=const_surv
          } else {
            male_survival=0
          }
        }
        
        if (is.na(male_survival)==TRUE){
          male_survival=0
        }
        
        if (diff_age_specific_fecundity==1){
          if (j>=lifetime$Maturity_M[line[i]]){
            if (lifetime$Fecundity_length_relationship[line[i]]=="Power"){
              male_fecundity=lifetime$Alpha[line[i]]*(lifetime$Length_M[line[i]+j-1]^(lifetime$Beta[line[i]]))      
            } else if (lifetime$Fecundity_length_relationship[line[i]]=="Exponential") {
              male_fecundity=lifetime$Alpha[line[i]]*exp(lifetime$Length_M[line[i]+j-1]*lifetime$Beta[line[i]])      
            } else if (lifetime$Fecundity_length_relationship[line[i]]=="Linear") {
              male_fecundity=lifetime$Alpha[line[i]]+(lifetime$Length_M[line[i]+j-1]*lifetime$Beta[line[i]])      
            }      
          } else {
            male_fecundity=0
          }
          if (is.na(male_fecundity)==TRUE){
            male_fecundity=0
          }
        } else {
          if (j>=lifetime$Maturity_M[line[i]]){
            male_fecundity=1
          } else {
            male_fecundity=0
          }
        }
        

        
        if (is.na(male_fecundity)==FALSE & male_fecundity==0){
          male_fecundity=0.01
        }
        
        if (is.na(male_fecundity)==TRUE){
          male_fecundity=0
        }
        
        PF_Male=lifetime$PF_M[line[i]]    
        
        if (lifetime$AgeNe_model[line[i]]=="TwoSex"){
          otherline=c(j,
                      round(female_survival,2),
                      round(female_fecundity,2),
                      PF_Female,
                      round(male_survival,2),
                      round(male_fecundity,2),
                      PF_Male)
        } else {
          sex_ratio=lifetime$SexRatio[line[i]+j-1]
          
          #if (female_survival==0){
          #  sex_ratio=1
          #} else if (male_survival==0){
          #  sex_ratio=0
          #}
          
          otherline=c(j,
                      round(female_survival,2),
                      round(female_fecundity,2),
                      PF_Female,
                      round(male_survival,2),
                      round(male_fecundity,2),
                      PF_Male,
                      sex_ratio)
          
        }
        
        surviv_plot[[i]]=c(surviv_plot[[i]],round(female_survival,2)*surviv_plot[[i]][length(surviv_plot[[i]])])
        #surviv_plot[[i]]=c(surviv_plot[[i]],round(female_survival,2))
        fecun_plot[[i]]=c(fecun_plot[[i]],round(female_fecundity,2))
        
        write(otherline,file="cogediv_sensibility.txt",append="TRUE",ncolumns=7)
        
      }
    } else {
      for (j in 1:lifetime$Max_age[line[i]]){
        
        if (diff_age_specific_survival==1){
          if (j<lifetime$Max_age_F[line[i]]){
            female_survival=exp(-((((lifetime$Length_F[line[i]+j-1])/(lifetime$Linf_F[line[i]]))^(-1.5))*lifetime$K_F[line[i]]))
          } else {
            female_survival=0
          }
        } else {
          const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[i]]-1))
          if (j<lifetime$Max_age_F[line[i]] & j>1){
            female_survival=const_surv
          } else if (j>1){
            female_survival=0
          } else {
            female_survival=1
          }
        }
        
        if (is.na(female_survival)==TRUE){
          female_survival=0
        }
        
        
        if (diff_age_specific_fecundity==1){
          if (j>=lifetime$Maturity_F[line[i]]){
            if (lifetime$Fecundity_length_relationship[line[i]]=="Power"){
              female_fecundity=lifetime$Alpha[line[i]]*(lifetime$Length_F[line[i]+j-1]^(lifetime$Beta[line[i]]))      
            } else if (lifetime$Fecundity_length_relationship[line[i]]=="Exponential") {
              female_fecundity=lifetime$Alpha[line[i]]*exp(lifetime$Length_F[line[i]+j-1]*lifetime$Beta[line[i]])      
            } else if (lifetime$Fecundity_length_relationship[line[i]]=="Linear") {
              female_fecundity=lifetime$Alpha[line[i]]+(lifetime$Length_F[line[i]+j-1]*lifetime$Beta[line[i]])      
            }      
          } else {
            female_fecundity=0
          }
        } else {
          female_fecundity=1
        }
        
        if (is.na(female_fecundity)==TRUE){
          female_fecundity=0
        }
        
        
        PF_Female=lifetime$PF_F[line[i]]
        
        if (diff_age_specific_survival==1){
          if (j<as.numeric(lifetime$Max_age_M[line[i]])){
            male_survival=exp(-((((lifetime$Length_M[line[i]+j-1])/(lifetime$Linf_M[line[i]]))^(-1.5))*lifetime$K_M[line[i]]))
          } else {
            male_survival=0
          }
        } else {
          const_surv=exp((log(0.01))/(as.numeric(lifetime$Max_age_M[line[i]])-1))
          if (j==1){
            male_survival=1
          } else if (j<as.numeric(lifetime$Max_age_M[line[i]])){
            male_survival=const_surv
          } else {
            male_survival=0
          }
        }
        
        
        if (is.na(male_survival)==TRUE){
          male_survival=0
        }
        
        if (diff_age_specific_fecundity==1){
          
          if (j>=lifetime$Maturity_F[line[i]]){
            if (lifetime$Fecundity_length_relationship[line[i]]=="Power"){
              male_fecundity=lifetime$Alpha[line[i]]*(lifetime$Length_M[line[i]+j-1]^(lifetime$Beta[line[i]]))      
            } else if (lifetime$Fecundity_length_relationship[line[i]]=="Exponential") {
              male_fecundity=lifetime$Alpha[line[i]]*exp(lifetime$Length_M[line[i]+j-1]*lifetime$Beta[line[i]])      
            } else if (lifetime$Fecundity_length_relationship[line[i]]=="Linear") {
              male_fecundity=lifetime$Alpha[line[i]]+(lifetime$Length_M[line[i]+j-1]*lifetime$Beta[line[i]])      
            }      
          } else {
            male_fecundity=0
          }
          if (is.na(male_fecundity)==TRUE){
            male_fecundity=0
          }
        } else {
          male_fecundity=1
        }
        
        if (is.na(male_fecundity)==TRUE){
          male_fecundity=0
        }
        
        PF_Male=lifetime$PF_M[line[i]]    
        
        sex_ratio=lifetime$SexRatio[line[i]+j-1]
        
        otherline=c(j,
                    round(female_survival,2),
                    round(female_fecundity,2),
                    PF_Female,
                    round(male_survival,2),
                    round(male_fecundity,2),
                    PF_Male,
                    sex_ratio)
        surviv_plot[[i]]=c(surviv_plot[[i]],round(female_survival,2)*surviv_plot[[i]][length(surviv_plot[[i]])])
        #surviv_plot[[i]]=c(surviv_plot[[i]],round(female_survival,2))
        fecun_plot[[i]]=c(fecun_plot[[i]],round(female_fecundity,2))
        write(otherline,file="cogediv_sensibility.txt",append="TRUE",ncolumns=8)
        
      }
    }
    
  }
  
  #readline(print(set_model[ll,]))
  # RUN AGENE -------
  system( shQuote( "AgeNe.exe"),input=c("cogediv_sensibility.txt","output_sensibility.txt"))

  ## EXTRACT INFORMATIONS ----
  agene<-readtext("output_sensibility.txt")
  agene<-strsplit(agene$text,"\n")[[1]]
  
  position=c()
  for (i in 1:n_output){
    position[i]=grep(latin[i],agene)
  }
  position=c(position,length(agene))
  
  for (i in 1:n_output){
    
    inf=position[i]
    sup=position[i+1]
    
    if (length(intersect(grep("overall",agene[inf:sup]),grep("Vk",agene[inf:sup])))>0){
      
      intersect=intersect(grep("overall",agene[inf:sup]),grep("Vk",agene[inf:sup]))
      agene_output[[1]][i,ll+3]=as.numeric(numextract(agene[inf:sup][intersect]))
      intersect=intersect(grep("female",agene[inf:sup]),grep("Vk",agene[inf:sup]))
      agene_output[[2]][i,ll+3]=as.numeric(numextract(agene[inf:sup][intersect]))
      intersect=intersect+1
      agene_output[[3]][i,ll+3]=as.numeric(numextract(agene[inf:sup][intersect]))
      intersect=intersect(grep("Ne",agene[inf:sup]),grep("Total N",agene[inf:sup]))
      agene_output[[4]][i,ll+3]=as.numeric(Numextract(agene[inf:sup][intersect]))[2]
      intersect=intersect(grep("Ne",agene[inf:sup]),grep("Adult N",agene[inf:sup]))
      agene_output[[5]][i,ll+3]=as.numeric(Numextract(agene[inf:sup][intersect]))[2]
      
    } else {
      
      intersect=grep("Vk",agene[inf:sup])
      agene_output[[1]][i,ll+3]=as.numeric(numextract(agene[inf:sup][intersect]))
      intersect=intersect(grep("Ne",agene[inf:sup]),grep("Total N",agene[inf:sup]))
      agene_output[[4]][i,ll+3]=as.numeric(Numextract(agene[inf:sup][intersect]))
      intersect=intersect(grep("Ne",agene[inf:sup]),grep("Adult N",agene[inf:sup]))
      agene_output[[5]][i,ll+3]=as.numeric(Numextract(agene[inf:sup][intersect]))
      
    }
    
  }
}

save(agene_output,file="agene_output.Rdata")
#load(file="C:/Users/pierr/Documents/PROJETS/COGEDIV/ARTICLE/Genetic_diversity_LHT/Data/agene_output.Rdata")















##### Analysis ------
wd="C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT"
setwd(wd)
load(file="Data/agene/agene_output.Rdata")

model_fit_Vk=as.data.frame(matrix(nrow=32,ncol=23,rep(0,32*23)))
colnames(model_fit_Vk)=c("slope-gaussian-freq","sd-gaussian-freq","pvalue-gaussian-freq","R-squared-gaussian-freq","AIC-gaussian-freq",
                         "slope-gaussian-bayesian","sd-gaussian-bayesian","DIC-gaussian-bayesian",
                         "slope-beta-freq","sd-beta-freq","pvalue-beta-freq","R-squared-beta-freq","AIC-beta-freq",
                         "slope-beta-bayesian","sd-beta-bayesian","DIC-beta-bayesian",
                         "slope-gamma-freq","sd-gamma-freq","pvalue-gamma-freq","AIC-gamma-freq",
                         "slope-gamma-bayesian","sd-gamma-bayesian","DIC-gamma-bayesian")

model_fit_NeN=as.data.frame(matrix(nrow=32,ncol=23,rep(0,32*23)))
colnames(model_fit_NeN)=c("slope-gaussian-freq","sd-gaussian-freq","pvalue-gaussian-freq","R-squared-gaussian-freq","AIC-gaussian-freq",
                         "slope-gaussian-bayesian","sd-gaussian-bayesian","DIC-gaussian-bayesian",
                         "slope-beta-freq","sd-beta-freq","pvalue-beta-freq","R-squared-beta-freq","AIC-beta-freq",
                         "slope-beta-bayesian","sd-beta-bayesian","DIC-beta-bayesian",
                         "slope-gamma-freq","sd-gamma-freq","pvalue-gamma-freq","AIC-gamma-freq",
                         "slope-gamma-bayesian","sd-gamma-bayesian","DIC-gamma-bayesian")

Vk=list(model_fit_Vk,model_fit_Vk)
NeN=list(model_fit_NeN,model_fit_NeN)

for (i in 1:5){
  agene_output[[i]]=agene_output[[i]][order(agene_output[[i]]$Species_code),]
}

for (ll in 1:16){
  ### Vk

  data_plot=data.frame(y=div/100,
                       x=agene_output[[1]][,ll+3],
                       x2=lfh$Parental_Care)
  
  # Gaussian frequentist - all data
  model<-lm(y~x,data=data_plot)
  Vk[[1]][ll,]$`slope-gaussian-freq`=model$coefficients[2]
  Vk[[1]][ll,]$`sd-gaussian-freq`=coefficients(summary(model))[2,2]
  Vk[[1]][ll,]$`pvalue-gaussian-freq`=coefficients(summary(model))[2,4]
  Vk[[1]][ll,]$`R-squared-gaussian-freq`=summary(model)$r.squared
  Vk[[1]][ll,]$`AIC-gaussian-freq`=AIC(model)
  # Gaussian frequentist - no parental care 
  model<-lm(y~x,data=data_plot[data_plot$x2=="No",])
  Vk[[2]][ll,]$`slope-gaussian-freq`=model$coefficients[2]
  Vk[[2]][ll,]$`sd-gaussian-freq`=coefficients(summary(model))[2,2]
  Vk[[2]][ll,]$`pvalue-gaussian-freq`=coefficients(summary(model))[2,4]
  Vk[[2]][ll,]$`R-squared-gaussian-freq`=summary(model)$r.squared
  Vk[[2]][ll,]$`AIC-gaussian-freq`=AIC(model)
  
  # Gaussian bayesian - all data
  model <- inla(y~x, family='gaussian',
                data=data_plot,
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  Vk[[1]][ll,]$`slope-gaussian-bayesian`=model$summary.fixed[2,1]
  Vk[[1]][ll,]$`sd-gaussian-bayesian`=model$summary.fixed[2,2]
  Vk[[1]][ll,]$`DIC-gaussian-bayesian`=model$dic$dic
  # Gaussian bayesian - no parental care
  model <- inla(y~x, family='gaussian',
                data=data_plot[data_plot$x2=="No",],
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  Vk[[2]][ll,]$`slope-gaussian-bayesian`=model$summary.fixed[2,1]
  Vk[[2]][ll,]$`sd-gaussian-bayesian`=model$summary.fixed[2,2]
  Vk[[2]][ll,]$`DIC-gaussian-bayesian`=model$dic$dic
  
  # Beta frequentist - all data
  model<-betareg(y~x,data=data_plot)
  Vk[[1]][ll,]$`slope-beta-freq`=model$coefficients$mean[2]
  Vk[[1]][ll,]$`sd-beta-freq`=coefficients(summary(model))$mean[2,2]
  Vk[[1]][ll,]$`pvalue-beta-freq`=coefficients(summary(model))$mean[2,4]
  Vk[[1]][ll,]$`R-squared-beta-freq`=model$pseudo.r.squared
  Vk[[1]][ll,]$`AIC-beta-freq`=AIC(model)
  # Beta frequentist - no parental care 
  model<-betareg(y~x,data=data_plot[data_plot$x2=="No",])
  Vk[[2]][ll,]$`slope-beta-freq`=model$coefficients$mean[2]
  Vk[[2]][ll,]$`sd-beta-freq`=coefficients(summary(model))$mean[2,2]
  Vk[[2]][ll,]$`pvalue-beta-freq`=coefficients(summary(model))$mean[2,4]
  Vk[[2]][ll,]$`R-squared-beta-freq`=model$pseudo.r.squared
  Vk[[2]][ll,]$`AIC-beta-freq`=AIC(model)
  
  # Beta bayesian - all data
  model <- inla(y~x, family='beta',
                data=data_plot,
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  Vk[[1]][ll,]$`slope-beta-bayesian`=model$summary.fixed[2,1]
  Vk[[1]][ll,]$`sd-beta-bayesian`=model$summary.fixed[2,2]
  Vk[[1]][ll,]$`DIC-beta-bayesian`=model$dic$dic
  # Beta bayesian - no parental care
  model <- inla(y~x, family='beta',
                data=data_plot[data_plot$x2=="No",],
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  Vk[[2]][ll,]$`slope-beta-bayesian`=model$summary.fixed[2,1]
  Vk[[2]][ll,]$`sd-beta-bayesian`=model$summary.fixed[2,2]
  Vk[[2]][ll,]$`DIC-beta-bayesian`=model$dic$dic
  
  # Gamma frequentist - all data
  model<-glm(y~x,family=Gamma(link='log'),data=data_plot)
  Vk[[1]][ll,]$`slope-gamma-freq`=coefficients(summary(model))[2,1]
  Vk[[1]][ll,]$`sd-gamma-freq`=coefficients(summary(model))[2,2]
  Vk[[1]][ll,]$`pvalue-gamma-freq`=coefficients(summary(model))[2,4]
  Vk[[1]][ll,]$`AIC-beta-freq`=AIC(model)
  # Gamma frequentist - no parental care 
  model<-glm(y~x,family=Gamma(link='log'),data=data_plot[data_plot$x2=="No",])
  Vk[[2]][ll,]$`slope-gamma-freq`=coefficients(summary(model))[2,1]
  Vk[[2]][ll,]$`sd-gamma-freq`=coefficients(summary(model))[2,2]
  Vk[[2]][ll,]$`pvalue-gamma-freq`=coefficients(summary(model))[2,4]
  Vk[[2]][ll,]$`AIC-beta-freq`=AIC(model)
  
  # Gamma bayesian - all data
  model <- inla(y~x, family='gamma',
                data=data_plot,
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  Vk[[1]][ll,]$`slope-gamma-bayesian`=model$summary.fixed[2,1]
  Vk[[1]][ll,]$`sd-gamma-bayesian`=model$summary.fixed[2,2]
  Vk[[1]][ll,]$`DIC-gamma-bayesian`=model$dic$dic
  # Gamma bayesian - no parental care
  model <- inla(y~x, family='gamma',
                data=data_plot[data_plot$x2=="No",],
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  Vk[[2]][ll,]$`slope-gamma-bayesian`=model$summary.fixed[2,1]
  Vk[[2]][ll,]$`sd-gamma-bayesian`=model$summary.fixed[2,2]
  Vk[[2]][ll,]$`DIC-gamma-bayesian`=model$dic$dic
  
  ### Ne/N
  
  data_plot=data.frame(y=div/100,
                       x=agene_output[[4]][,ll+3],
                       x2=lfh$Parental_Care)
  
  # Gaussian frequentist - all data
  model<-lm(y~x,data=data_plot)
  NeN[[1]][ll,]$`slope-gaussian-freq`=model$coefficients[2]
  NeN[[1]][ll,]$`sd-gaussian-freq`=coefficients(summary(model))[2,2]
  NeN[[1]][ll,]$`pvalue-gaussian-freq`=coefficients(summary(model))[2,4]
  NeN[[1]][ll,]$`R-squared-gaussian-freq`=summary(model)$r.squared
  NeN[[1]][ll,]$`AIC-gaussian-freq`=AIC(model)
  # Gaussian frequentist - no parental care 
  model<-lm(y~x,data=data_plot[data_plot$x2=="No",])
  NeN[[2]][ll,]$`slope-gaussian-freq`=model$coefficients[2]
  NeN[[2]][ll,]$`sd-gaussian-freq`=coefficients(summary(model))[2,2]
  NeN[[2]][ll,]$`pvalue-gaussian-freq`=coefficients(summary(model))[2,4]
  NeN[[2]][ll,]$`R-squared-gaussian-freq`=summary(model)$r.squared
  NeN[[2]][ll,]$`AIC-gaussian-freq`=AIC(model)
  
  # Gaussian bayesian - all data
  model <- inla(y~x, family='gaussian',
                data=data_plot,
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  NeN[[1]][ll,]$`slope-gaussian-bayesian`=model$summary.fixed[2,1]
  NeN[[1]][ll,]$`sd-gaussian-bayesian`=model$summary.fixed[2,2]
  NeN[[1]][ll,]$`DIC-gaussian-bayesian`=model$dic$dic
  # Gaussian bayesian - no parental care
  model <- inla(y~x, family='gaussian',
                data=data_plot[data_plot$x2=="No",],
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  NeN[[2]][ll,]$`slope-gaussian-bayesian`=model$summary.fixed[2,1]
  NeN[[2]][ll,]$`sd-gaussian-bayesian`=model$summary.fixed[2,2]
  NeN[[2]][ll,]$`DIC-gaussian-bayesian`=model$dic$dic
  
  # Beta frequentist - all data
  model<-betareg(y~x,data=data_plot)
  NeN[[1]][ll,]$`slope-beta-freq`=model$coefficients$mean[2]
  NeN[[1]][ll,]$`sd-beta-freq`=coefficients(summary(model))$mean[2,2]
  NeN[[1]][ll,]$`pvalue-beta-freq`=coefficients(summary(model))$mean[2,4]
  NeN[[1]][ll,]$`R-squared-beta-freq`=model$pseudo.r.squared
  NeN[[1]][ll,]$`AIC-beta-freq`=AIC(model)
  # Beta frequentist - no parental care 
  model<-betareg(y~x,data=data_plot[data_plot$x2=="No",])
  NeN[[2]][ll,]$`slope-beta-freq`=model$coefficients$mean[2]
  NeN[[2]][ll,]$`sd-beta-freq`=coefficients(summary(model))$mean[2,2]
  NeN[[2]][ll,]$`pvalue-beta-freq`=coefficients(summary(model))$mean[2,4]
  NeN[[2]][ll,]$`R-squared-beta-freq`=model$pseudo.r.squared
  NeN[[2]][ll,]$`AIC-beta-freq`=AIC(model)
  
  # Beta bayesian - all data
  model <- inla(y~x, family='beta',
                data=data_plot,
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  NeN[[1]][ll,]$`slope-beta-bayesian`=model$summary.fixed[2,1]
  NeN[[1]][ll,]$`sd-beta-bayesian`=model$summary.fixed[2,2]
  NeN[[1]][ll,]$`DIC-beta-bayesian`=model$dic$dic
  # Beta bayesian - no parental care
  model <- inla(y~x, family='beta',
                data=data_plot[data_plot$x2=="No",],
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  NeN[[2]][ll,]$`slope-beta-bayesian`=model$summary.fixed[2,1]
  NeN[[2]][ll,]$`sd-beta-bayesian`=model$summary.fixed[2,2]
  NeN[[2]][ll,]$`DIC-beta-bayesian`=model$dic$dic
  
  # Gamma frequentist - all data
  model<-glm(y~x,family=Gamma(link='log'),data=data_plot)
  NeN[[1]][ll,]$`slope-gamma-freq`=coefficients(summary(model))[2,1]
  NeN[[1]][ll,]$`sd-gamma-freq`=coefficients(summary(model))[2,2]
  NeN[[1]][ll,]$`pvalue-gamma-freq`=coefficients(summary(model))[2,4]
  NeN[[1]][ll,]$`AIC-beta-freq`=AIC(model)
  # Gamma frequentist - no parental care 
  model<-glm(y~x,family=Gamma(link='log'),data=data_plot[data_plot$x2=="No",])
  NeN[[2]][ll,]$`slope-gamma-freq`=coefficients(summary(model))[2,1]
  NeN[[2]][ll,]$`sd-gamma-freq`=coefficients(summary(model))[2,2]
  NeN[[2]][ll,]$`pvalue-gamma-freq`=coefficients(summary(model))[2,4]
  NeN[[2]][ll,]$`AIC-beta-freq`=AIC(model)
  
  # Gamma bayesian - all data
  model <- inla(y~x, family='gamma',
                data=data_plot,
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  NeN[[1]][ll,]$`slope-gamma-bayesian`=model$summary.fixed[2,1]
  NeN[[1]][ll,]$`sd-gamma-bayesian`=model$summary.fixed[2,2]
  NeN[[1]][ll,]$`DIC-gamma-bayesian`=model$dic$dic
  # Gamma bayesian - no parental care
  model <- inla(y~x, family='gamma',
                data=data_plot[data_plot$x2=="No",],
                #control.family=list(link='logit'),
                control.predictor=list(link=1, compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
  NeN[[2]][ll,]$`slope-gamma-bayesian`=model$summary.fixed[2,1]
  NeN[[2]][ll,]$`sd-gamma-bayesian`=model$summary.fixed[2,2]
  NeN[[2]][ll,]$`DIC-gamma-bayesian`=model$dic$dic
  
}

#### Fit model -----
for (ll in 1:16){
 
  data_plot=cbind(data_plot,model$summary.fitted.values)
  p<-ggplot(data_plot, 
            aes(x = x, 
                y = y,
                label=lfh$Species_plot,col=x2)) +
    geom_line(aes(x=data_plot$x,
                  y=data_plot$mean),
              col="black")+
    geom_ribbon(aes(ymin=data_plot$`0.025quant`,
                    ymax=data_plot$`0.975quant`),
                fill='dodgerblue2',
                alpha=0.25,
                col="white") +
    xlab("Variance in reproductive success")+
    ylab("Heterozygosity (%)")+
    geom_rangeframe()+
    theme_classic()+
    geom_text_repel()+
    theme(legend.position = "none")+
    geom_point(size=2,
               pch=19
               #aes(col=x2)
    ) +
    #scale_color_manual(values=brewer.pal(2,"Spectral")[c(1,3)])
    scale_colour_viridis_d(end=0.7)
  
  p<- p + annotate(geom="text", 
                   x=0.90*max(data_plot$x), 
                   y=0.95*max(data_plot$y),
                   label=paste("p-value : ",round(coefficients(summary(reg))[2,4],3),sep=""),
                   color="black")
  p<- p + annotate(geom="text", 
                   x=0.90*max(data_plot$x), 
                   y=0.90*max(data_plot$y),
                   label=paste("R² : ",round(summary(reg)$r.squared,3),sep=""),
                   color="black")
  
  print(p)
 
  
}
work_complete()

wd="C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT/images"
setwd(wd)

sim=16
compare_lfh_genet(x=agene_output[[4]][,sim+3],
                  y=div,
                  likelihood='gaussian',
                  x_legend="Ne/N",
                  export_image="agene16",
                  lfh=lfh)









print(l<-summary(lm(compare_het_stat$ratio_mean~compare_het_stat$ratio_truehet)))
r.x <- lm(compare_het_stat$ratio_mean ~ compare_het_stat$ratio_truehet)
r1 <- lm(compare_het_stat$ratio_mean ~1 + offset(compare_het_stat$ratio_truehet))
anova(r.x, r1)




delta_agemat=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==1 &
        set_model$AgeMat[i]==0 & 
        set_model$Surv[j]==set_model$Surv[i] & 
        set_model$Fec[j]==set_model$Fec[i] & 
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_agemat=c(delta_agemat,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_surv=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] & 
        set_model$Surv[j]==1 &
        set_model$Surv[i]==0 & 
        set_model$Fec[j]==set_model$Fec[i] & 
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_surv=c(delta_surv,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_fec=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] & 
        set_model$Surv[j]==set_model$Surv[i] & 
        set_model$Fec[j]==1 &
        set_model$Fec[i]==0 & 
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_fec=c(delta_fec,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_sex=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] & 
        set_model$Surv[j]==set_model$Surv[i] & 
        set_model$Fec[j]==set_model$Fec[i] & 
        set_model$Sex[j]==1 &
        set_model$Sex[i]==0 &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_sex=c(delta_sex,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_hermaphro=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] & 
        set_model$Surv[j]==set_model$Surv[i] & 
        set_model$Fec[j]==set_model$Fec[i] & 
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==1 &
        set_model$Hermaphro[i]==0){
      delta_hermaphro=c(delta_hermaphro,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_agematsurv=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==1 &
        set_model$AgeMat[i]==0 & 
        set_model$Surv[j]==1 &
        set_model$Surv[i]==0 & 
        set_model$Fec[j]==set_model$Fec[i] & 
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_agematsurv=c(delta_agematsurv,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_agematfec=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==1 &
        set_model$AgeMat[i]==0 & 
        set_model$Surv[j]==set_model$Surv[i] &
        set_model$Fec[j]==1 &
        set_model$Fec[i]==0 & 
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_agematfec=c(delta_agematfec,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_agematsex=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==1 &
        set_model$AgeMat[i]==0 & 
        set_model$Surv[j]==set_model$Surv[i] &
        set_model$Fec[j]==set_model$Fec[i] &
        set_model$Sex[j]==1 &
        set_model$Sex[i]==0 &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_agematsex=c(delta_agematsex,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_agemathermaphro=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==1 &
        set_model$AgeMat[i]==0 & 
        set_model$Surv[j]==set_model$Surv[i] &
        set_model$Fec[j]==set_model$Fec[i] &
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==1 & 
        set_model$Hermaphro[i]==0){
      delta_agemathermaphro=c(delta_agemathermaphro,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_survfec=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] &
        set_model$Surv[j]==1 &
        set_model$Surv[i]==0 &
        set_model$Fec[j]==1 &
        set_model$Fec[i]==0 &
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_survfec=c(delta_survfec,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_survsex=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] &
        set_model$Surv[j]==1 &
        set_model$Surv[i]==0 &
        set_model$Fec[j]==set_model$Fec[i] &
        set_model$Sex[j]==1 &
        set_model$Sex[i]==0 &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_survsex=c(delta_survsex,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_survhermahpro=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] &
        set_model$Surv[j]==1 &
        set_model$Surv[i]==0 &
        set_model$Fec[j]==set_model$Fec[i] &
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==1 & 
        set_model$Hermaphro[i]==0){
      delta_survhermahpro=c(delta_survhermahpro,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_fecsex=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] &
        set_model$Surv[j]==set_model$Surv[i] &
        set_model$Fec[j]==1 &
        set_model$Fec[i]==0 &
        set_model$Sex[j]==1 &
        set_model$Sex[i]==0 &
        set_model$Hermaphro[j]==set_model$Hermaphro[i]){
      delta_fecsex=c(delta_fecsex,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_fechermahpro=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] &
        set_model$Surv[j]==set_model$Surv[i] &
        set_model$Fec[j]==1 &
        set_model$Fec[i]==0 &
        set_model$Sex[j]==set_model$Sex[i] &
        set_model$Hermaphro[j]==1 &
        set_model$Hermaphro[i]==0){
      delta_fechermahpro=c(delta_fechermahpro,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}
delta_sexhermahpro=c()
for (i in 1:(nrow(set_model)-1)){
  for (j in (i+1):nrow(set_model)){
    if (set_model$AgeMat[j]==set_model$AgeMat[i] &
        set_model$Surv[j]==set_model$Surv[i] &
        set_model$Fec[j]==set_model$Fec[i] &
        set_model$Sex[j]==1 &
        set_model$Sex[i]==0 &
        set_model$Hermaphro[j]==1 &
        set_model$Hermaphro[i]==0){
      delta_sexhermahpro=c(delta_sexhermahpro,model_fit_Vk$DIC[j]-model_fit_Vk$DIC[i])
    }
  }
}

data_delta=data.frame(value=c(delta_agemat,
                              delta_surv,
                              delta_fec,
                              delta_sex,
                              delta_hermaphro,
                              delta_agematsurv,
                              delta_agematfec,
                              delta_agematsex,
                              delta_agemathermaphro,
                              delta_survfec,
                              delta_survsex,
                              delta_survhermahpro,
                              delta_fecsex,
                              delta_fechermahpro,
                              delta_sexhermahpro),
                      fact=c(rep("delta_agemat",16),
                             rep("delta_surv",16),
                             rep("delta_fec",16),
                             rep("delta_sex",16),
                             rep("delta_hermaphro",16),
                             rep("delta_agematsurv",8),
                             rep("delta_agematfec",8),
                             rep("delta_agematsex",8),
                             rep("delta_agemathermaphro",8),
                             rep("delta_survfec",8),
                             rep("delta_survsex",8),
                             rep("delta_survhermahpro",8),
                             rep("delta_fecsex",8),
                             rep("delta_fechermahpro",8),
                             rep("delta_sexhermahpro",8)))
ggplot(data_delta,aes(y=value,x=fact))+
  geom_boxplot()+
  geom_jitter()


data_Vk=cbind(model_fit_Vk[1:16,],set_model[1:16,])
colnames(data_Vk)=c("slope","pvalue","rsquared","AIC","BIC","DIC","WAIC","MLIK","AgeMat","Surv","Fec","Sex","Hermaphro")
summary(lm(pvalue~Fec+AgeMat+Surv+Sex,data=data_Vk))
data_NeN=cbind(model_fit_NeN[1:16,],set_model[1:16,])
colnames(data_NeN)=c("slope","pvalue","rsquared","AIC","BIC","DIC","WAIC","MLIK","AgeMat","Surv","Fec","Sex","Hermaphro")
summary(lm(BIC~Fec+AgeMat+Surv+Sex,data=data_NeN))


#------------
p <- ggplot(fastp_process, aes(x=SPECIES, y=fastp_process[,i]),color=color) + 
  geom_boxplot(outlier.shape = NA,alpha=0)+
  geom_violin(color="black",alpha=0) +
  geom_jitter(position=position_jitter(0.2),color=color,pch=redo,size=0.00001,alpha=0)+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  #geom_errorbar(aes(ymin=Summary_GenomeScope[,i-1], ymax=Summary_GenomeScope[,i+1]), width=.1,position=position_dodge(width=0.5)) +
  ylab(colnames(fastp_process)[i])+
  ggtitle(paste(analogy_legend[i]))+
  #geom_text_repel(sub<-subset(fastp_process,fastp_process[,i]<quantile(fastp_process[,i],probs=c(0.025)) | 
  #                              fastp_process[,i]>quantile(fastp_process[,i],probs=c(0.975))), 
  #                mapping=aes(x=SPECIES, y=sub[,i],label=SAMPLE),color="black")+
  scale_x_discrete(name="Species",
                   labels=labels_species)+
  ylab(analogy_legend[i])+
  theme_tufte()+
  ylim(c(min(fastp_process[,i]),
         max(fastp_process[,i])+max(fastp_process[,i])*(gross*1)))+
  geom_rangeframe()+
  geom_point_interactive(aes(x=SPECIES, y=fastp_process[,i],
                             tooltip=SAMPLE,
                             data_id=Num_batch), size = 1.5,position="jitter",
                         color=color,pch=redo,alpha=0.75) +
  scale_fill_discrete(name="Location",
                      breaks=c("Li", "Mu", "Fa","Ga"),
                      labels=c("Gulf of Lion","Murcia","Faro","Gulf of Gascogne")) +
  theme(legend.position="top")
p<-ggMarginal(p, type="histogram",margins='y',alpha=0.25,col="black",fill="yellow",cex=0.5)

# Sensibility test: only lifespan ----
for (i in 1:length(species)){
  
  
  #First line
  firstline=paste(lifetime$Latin[line[i]],"Lifespan-test")
  if (file.exists("cogediv_lifespan.txt")==FALSE){
    write(firstline,file="cogediv_lifespan.txt",append="FALSE")
  } else {
    write(firstline,file="cogediv_lifespan.txt",append="TRUE")
  }
  
  
  # Second line
  secondline=c(lifetime$Max_age[line[i]],
               "1000",
               0.5)
  write(secondline,file="cogediv_lifespan.txt",append="TRUE",ncolumns = 3)
  
  
  # Other line
  for (j in 1:lifetime$Max_age[line[i]]){
    
    const_surv=exp((log(0.01))/(lifetime$Max_age[line[i]]-1))
    if (j<lifetime$Max_age_F[line[i]] & j>1){
      female_survival=const_surv
    } else if (j>1){
      female_survival=0
    } else {
      female_survival=1
    }
    
    female_fecundity=1
    
    PF_Female=1
    
    const_surv=exp((log(0.01))/(lifetime$Max_age[line[i]]-1))
    if (j<lifetime$Max_age_F[line[i]] & j>1){
      male_survival=const_surv
    } else if (j>1){
      male_survival=0
    } else {
      male_survival=1
    }
    
    male_fecundity=1
    
    PF_Male=1
    
    otherline=c(j,
                round(female_survival,2),
                round(female_fecundity,2),
                PF_Female,
                round(male_survival,2),
                round(male_fecundity,2),
                PF_Male)
    write(otherline,file="cogediv_lifespan.txt",append="TRUE",ncolumns=7)
    
  }
  
  
  
  
}
system( shQuote( "AgeNe.exe"),input=c("cogediv_lifespan.txt","output_lifespan.txt"))
agene_output_lifespan=data.frame(Species_code=rep(NA,n_output),
                                 Latin=rep(NA,n_output),
                                 Vernacular=rep(NA,n_output),
                                 Vk=rep(0,n_output),
                                 Vk_female=rep(0,n_output),
                                 Vk_male=rep(0,n_output),
                                 Ne_N=rep(0,n_output),
                                 Ne_AdultN=rep(0,n_output)
)

agene_output_lifespan$Species_code=species
agene_output_lifespan$Latin=latin
agene_output_lifespan$Vernacular=vernacular

agene<-readtext("C:/Users/pierr/Desktop/output_lifespan.txt")
agene<-strsplit(agene$text,"\n")[[1]]

position=c()
for (i in 1:n_output){
  position[i]=grep(latin[i],agene)
}
position=c(position,length(agene))

for (i in 1:n_output){
  
  inf=position[i]
  sup=position[i+1]
  
  if (length(intersect(grep("overall",agene[inf:sup]),grep("Vk",agene[inf:sup])))>0){
    
    intersect=intersect(grep("overall",agene[inf:sup]),grep("Vk",agene[inf:sup]))
    agene_output_lifespan$Vk[i]=as.numeric(numextract(agene[inf:sup][intersect]))
    intersect=intersect(grep("female",agene[inf:sup]),grep("Vk",agene[inf:sup]))
    agene_output_lifespan$Vk_female[i]=as.numeric(numextract(agene[inf:sup][intersect]))
    intersect=intersect+1
    agene_output_lifespan$Vk_male[i]=as.numeric(numextract(agene[inf:sup][intersect]))
    intersect=intersect(grep("Ne",agene[inf:sup]),grep("Total N",agene[inf:sup]))
    agene_output_lifespan$Ne_N[i]=as.numeric(Numextract(agene[inf:sup][intersect]))[2]
    intersect=intersect(grep("Ne",agene[inf:sup]),grep("Adult N",agene[inf:sup]))
    agene_output_lifespan$Ne_AdultN[i]=as.numeric(Numextract(agene[inf:sup][intersect]))[2]
    
  } else {
    
    intersect=grep("Vk",agene[inf:sup])
    agene_output_lifespan$Vk[i]=as.numeric(numextract(agene[inf:sup][intersect]))
    intersect=intersect(grep("Ne",agene[inf:sup]),grep("Total N",agene[inf:sup]))
    agene_output_lifespan$Ne_N[i]=as.numeric(Numextract(agene[inf:sup][intersect]))
    intersect=intersect(grep("Ne",agene[inf:sup]),grep("Adult N",agene[inf:sup]))
    agene_output_lifespan$Ne_AdultN[i]=as.numeric(Numextract(agene[inf:sup][intersect]))
    
  }
  
}

print(agene_output_lifespan)
setwd("C:/Users/pierr/Documents/PROJETS/COGEDIV/BIOINFO/AgeNe")
save(agene_output_lifespan,file="agene_output_lifespan.Rdata")


#link to compare lfh_genetic
# Sensibility test : only increaseing survival rates ----
for (i in 1:length(species)){
  
  
  #First line
  firstline=paste(lifetime$Latin[line[i]],"Lifespan-survival-test")
  if (file.exists("cogediv_lifespan_survival.txt")==FALSE){
    write(firstline,file="cogediv_lifespan_survival.txt",append="FALSE")
  } else {
    write(firstline,file="cogediv_lifespan_survival.txt",append="TRUE")
  }
  
  
  # Second line
  secondline=c(lifetime$Max_age[line[i]],
               "1000",
               0.5)
  write(secondline,file="cogediv_lifespan_survival.txt",append="TRUE",ncolumns = 3)
  
  
  # Other line
  for (j in 1:lifetime$Max_age[line[i]]){
    
    if (j<lifetime$Max_age_F[line[i]]){
      female_survival=exp(-((((lifetime$Length_F[line[i]+j-1])/(lifetime$Linf_F[line[i]]))^(-1.5))*lifetime$K_F[line[i]]))
    } else {
      female_survival=0
    }
    
    female_fecundity=1
    PF_Female=1
    
    if (j<lifetime$Max_age_M[line[i]]){
      male_survival=exp(-((((lifetime$Length_M[line[i]+j-1])/(lifetime$Linf_M[line[i]]))^(-1.5))*lifetime$K_M[line[i]]))
    } else {
      male_survival=0
    }
    
    male_fecundity=1
    PF_Male=1    
    
    
    otherline=c(j,
                round(female_survival,2),
                round(female_fecundity,2),
                PF_Female,
                round(male_survival,2),
                round(male_fecundity,2),
                PF_Male)
    
    write(otherline,file="cogediv_lifespan_survival.txt",append="TRUE",ncolumns=7)
    
  }
  
  
}
n_output=length(species)
system( shQuote( "AgeNe.exe"),input=c("cogediv_lifespan_survival.txt","output_lifespan_survival.txt"))
agene_output_lifespan_survival=data.frame(Species_code=rep(NA,n_output),
                                          Latin=rep(NA,n_output),
                                          Vernacular=rep(NA,n_output),
                                          Vk=rep(0,n_output),
                                          Vk_female=rep(0,n_output),
                                          Vk_male=rep(0,n_output),
                                          Ne_N=rep(0,n_output),
                                          Ne_AdultN=rep(0,n_output)
)

agene_output_lifespan_survival$Species_code=species
agene_output_lifespan_survival$Latin=latin
agene_output_lifespan_survival$Vernacular=vernacular

agene<-readtext("C:/Users/pierr/Desktop/output_lifespan_survival.txt")
agene<-strsplit(agene$text,"\n")[[1]]

position=c()
for (i in 1:n_output){
  position[i]=grep(latin[i],agene)
}
position=c(position,length(agene))

for (i in 1:n_output){
  
  inf=position[i]
  sup=position[i+1]
  
  if (length(intersect(grep("overall",agene[inf:sup]),grep("Vk",agene[inf:sup])))>0){
    
    intersect=intersect(grep("overall",agene[inf:sup]),grep("Vk",agene[inf:sup]))
    agene_output_lifespan_survival$Vk[i]=as.numeric(numextract(agene[inf:sup][intersect]))
    intersect=intersect(grep("female",agene[inf:sup]),grep("Vk",agene[inf:sup]))
    agene_output_lifespan_survival$Vk_female[i]=as.numeric(numextract(agene[inf:sup][intersect]))
    intersect=intersect+1
    agene_output_lifespan_survival$Vk_male[i]=as.numeric(numextract(agene[inf:sup][intersect]))
    intersect=intersect(grep("Ne",agene[inf:sup]),grep("Total N",agene[inf:sup]))
    agene_output_lifespan_survival$Ne_N[i]=as.numeric(Numextract(agene[inf:sup][intersect]))[2]
    intersect=intersect(grep("Ne",agene[inf:sup]),grep("Adult N",agene[inf:sup]))
    agene_output_lifespan_survival$Ne_AdultN[i]=as.numeric(Numextract(agene[inf:sup][intersect]))[2]
    
  } else {
    
    intersect=grep("Vk",agene[inf:sup])
    agene_output_lifespan_survival$Vk[i]=as.numeric(numextract(agene[inf:sup][intersect]))
    intersect=intersect(grep("Ne",agene[inf:sup]),grep("Total N",agene[inf:sup]))
    agene_output_lifespan_survival$Ne_N[i]=as.numeric(Numextract(agene[inf:sup][intersect]))
    intersect=intersect(grep("Ne",agene[inf:sup]),grep("Adult N",agene[inf:sup]))
    agene_output_lifespan_survival$Ne_AdultN[i]=as.numeric(Numextract(agene[inf:sup][intersect]))
    
  }
  
}

print(agene_output_lifespan_survival)
setwd("C:/Users/pierr/Documents/PROJETS/COGEDIV/BIOINFO/AgeNe")
save(agene_output_lifespan_survival,file="agene_output_lifespan_survival.Rdata")


#ANALYSE -----
load(file="C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/concatenate_species.Rdata")
cogediv_species=concatenate_species()

het=rep(NA,21)

for (i in 1:16){
  for (j in 1:21){
    if (agene_output_lifespan_survival$Species_code[i]==cogediv_species$Species_code[j]){
      het[j]=cogediv_species$Heterozygosity_median[j]
    }
    
  }
}
het=het[is.na(het)==FALSE]
plot(agene_output_lifespan_survival$Vk_male,het)

for (i in 1:length(species)){
  for (j in 1:length(fecun_plot[[i]])){
    fecun_plot[[i]][j]=fecun_plot[[i]][j]/fecun_plot[[i]][length(fecun_plot[[i]])]
  }
}

colfunc<-colorRampPalette(c("royalblue","springgreen","yellow","red"))
het_color=floor(rescale(het,c(1,100000)))


plot(seq(0,21,length.out=100),seq(0,1,length.out=100),
     type="n",
     xlab="",
     ylab="",
     log=c("x"),
     xlim=c(1,20))
# Plot survival curves
for (i in 1:length(species)){
  lines(surviv_plot[[i]],
        lwd=4,
        col=colfunc(100000-1)[het_color[i]])
}
plot(seq(0,21,length.out=100),seq(0,1,length.out=100),
     type="n",
     xlab="",
     ylab="",
     log=c("x"),
     xlim=c(1,20))
for (i in 1:length(species)){
  lines(fecun_plot[[i]],
        lwd=4,col=colfunc(100000-1)[het_color[i]])
}

# Does Vk is linked to heterozygosity ?
data_plot=data.frame(species=lfh$Species_code,genetic=genetic$Heterozygosity_median,
                     genetic_sd=genetic[,gen+1],
                     trait=lfh$Ne_N*lfh$Trophic_Level,
                     trait2=lfh$Parental_Care)

data_plot=data_plot[data_plot$trait2=="No",]
data_plot=data_plot[is.na(data_plot$genetic)==FALSE &
                      is.na(data_plot$trait)==FALSE,]

p<-ggplot(data_plot, aes(x = trait, y = genetic,label=species)) +
  geom_point(size=2,pch=19) +
  #geom_errorbar(aes(ymin = genetic[,gen]-genetic[,gen+1], ymax = genetic[,gen]+genetic[,gen+1]),position='dodge',width=max(lfh[,trait],na.rm=TRUE)/80)+
  #xlab("Variance in reproductive success")+ylab("Heterozygosity")+
  ggtitle(paste("Correlation between ","Variance in reproductive success"," and ","Heterozygosity",sep=""))+
  #theme_tufte()+
  #geom_rangeframe()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  geom_text_repel()


model<-betareg(genetic~trait,data=data_plot)
print(summary(model))

p<-p+
  geom_line(aes(y = predict(model, data_plot)), 
            color="red", linetype = "solid",size=0.5)

p<-p+
  geom_text(x=max(data_plot$trait,na.rm=TRUE)*0.9, y=max(data_plot$genetic,na.rm=TRUE),label=paste("p-value = ",round(coef(summary(model))$mean[2,4],4),sep=""))+
  geom_text(x=max(data_plot$trait,na.rm=TRUE)*0.9, y=max(data_plot$genetic,na.rm=TRUE)*0.95,label=paste("Pseudo R² = ",round(model$pseudo.r.squared,3),sep=""))


print(p)

#○ Does Vk is linked to heterozygosity ?
data_plot=data.frame(species=lfh$Species_code,genetic=genetic$Heterozygosity_median,
                     genetic_sd=genetic[,gen+1],
                     trait=lfh$Ne_N)

data_plot=data_plot[is.na(data_plot$genetic)==FALSE &
                      is.na(data_plot$trait)==FALSE,]

p<-ggplot(data_plot, aes(x = trait, y = genetic,label=species)) +
  geom_point(size=2,pch=19) +
  #geom_errorbar(aes(ymin = genetic[,gen]-genetic[,gen+1], ymax = genetic[,gen]+genetic[,gen+1]),position='dodge',width=max(lfh[,trait],na.rm=TRUE)/80)+
  xlab("Ne/N")+ylab("Heterozygosity")+
  ggtitle(paste("Correlation between ","Ne/N"," and ","Heterozygosity",sep=""))+
  #theme_tufte()+
  #geom_rangeframe()+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  geom_text_repel()


model<-betareg(genetic~trait,data=data_plot)
print(summary(model))

p<-p+
  geom_line(aes(y = predict(model, data_plot)), 
            color="red", linetype = "solid",size=0.5)

p<-p+
  geom_text(x=max(data_plot$trait,na.rm=TRUE)*0.9, y=max(data_plot$genetic,na.rm=TRUE),label=paste("p-value = ",round(coef(summary(model))$mean[2,4],4),sep=""))+
  geom_text(x=max(data_plot$trait,na.rm=TRUE)*0.9, y=max(data_plot$genetic,na.rm=TRUE)*0.95,label=paste("Pseudo R² = ",round(model$pseudo.r.squared,3),sep=""))


print(p)

plot(lfh$Ne_N,
     genetic$Heterozygosity_median)

