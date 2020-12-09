#------------------------------------------------------------------------#
#                                                                        #
#                             Run AgeNe                                  #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages ----
library(readxl)
library(readtext)
library(stringr)
# Load functions ----
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 
Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

# Set working directory and initialize  ----
setwd("Data/agene")
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
# Construct and run agene ----
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

save_new_output=0
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
        fecun_plot[[i]]=c(fecun_plot[[i]],round(female_fecundity,2))
        write(otherline,file="cogediv_sensibility.txt",append="TRUE",ncolumns=8)
        
      }
    }
    
  }
  
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

# Save new output (if wanted) ----
if (save_new_output==1){
  save(agene_output,file="agene_output.Rdata")
}




