#------------------------------------------------------------------------#
#                                                                        #
#                        Get SLiM input                                  #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages ----
library(readxl)
library(readtext)
# Get life tables data ----
species=c("Spilc",
          "Hgutt",
          "Mmerl",
          "Scabr",
          "Dlabr",
          "Msurm",
          "Lmorm",
          "Dpunt",
          "Peryt",
          "Cjuli",
          "Ssard",
          "Cgale",
          "Scine",
          "Lbude",
          "Styph",
          "Scant")
max_fec=100

lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
line=c(line,nrow(lifetime)+1)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
line=line[-17]


# Get input for SLiM simulations ----------
## 1. K = 500 | Surv = / | Fec = / | Sex = / | AgeMat = / ---------
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input1/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input1/testdiv_",species[sp],".txt",sep=""))
  }
}

lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input1",sep=""))
for (sp in 1:length(species)){
  lifetime$Maturity_F[line[sp]]=1
  lifetime$Maturity_M[line[sp]]=1
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv
  surv_grep=grep(c("Surv"),divforSlim)
  female_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[sp]]))
  
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j==(lifetime$Max_age_F[line[sp]]+2)){
      female_survival=paste(female_survival,1,")",sep="")
    } else {
      female_survival=paste(female_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv",',female_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<=lifetime$Max_age_F[line[sp]] & j>lifetime$Maturity_F[line[sp]]){
      female_fecundity=paste(female_fecundity,max_fec,",",sep="")
      
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)) {
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
      
    } else {
      female_fecundity=paste(female_fecundity,0,",",sep="")
      
    }
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim//Output1/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output1/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output1/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}







## 2. K = 500 | Surv = Yes | Fec = /| Sex = / | AgeMat = / ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input2/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input2/testdiv_",species[sp],".txt",sep=""))
  }
}

lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input2",sep=""))
for (sp in 1:length(species)){
  lifetime$Maturity_F[line[sp]]=1
  lifetime$Maturity_M[line[sp]]=1
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv
  surv_grep=grep(c("Surv"),divforSlim)
  female_survival="c("
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j<(lifetime$Max_age_F[line[sp]]+2)){
      female_survival=paste(female_survival,round(1-(exp(-((((lifetime$Length_F[line[sp]+j-1])/(lifetime$Linf_F[line[sp]]))^(-1.5))*lifetime$K_F[line[sp]]))),2),",",sep="")
    } else {
      female_survival=paste(female_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv",',female_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<=lifetime$Max_age_F[line[sp]] & j>lifetime$Maturity_F[line[sp]]){
      female_fecundity=paste(female_fecundity,max_fec,",",sep="")
      
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)) {
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
      
    } else {
      female_fecundity=paste(female_fecundity,0,",",sep="")
      
    }
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output2/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output2/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output2/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}











## 3. K = 500 | Surv = / | Fec = Yes | Sex = / | AgeMat = / ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input3/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input3/testdiv_",species[sp],".txt",sep=""))
  }
}

lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input3",sep=""))
for (sp in 1:length(species)){
  lifetime$Maturity_F[line[sp]]=1
  lifetime$Maturity_M[line[sp]]=1
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv
  surv_grep=grep(c("Surv"),divforSlim)
  female_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[sp]]))
  
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j==(lifetime$Max_age_F[line[sp]]+2)){
      female_survival=paste(female_survival,1,")",sep="")
    } else {
      female_survival=paste(female_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv",',female_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j>=(lifetime$Maturity_F[line[sp]]+1) & j<(lifetime$Max_age_F[line[sp]]+1)){
      if (lifetime$Fecundity_length_relationship[line[sp]]=="Power"){
        fec_inter=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+j-1]^(lifetime$Beta[line[sp]]))   
        maxfec=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+lifetime$Max_age_F[line[sp]]]^(lifetime$Beta[line[sp]]))  
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Exponential") {
        fec_inter=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]])  
        maxfec=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])    
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Linear") {
        fec_inter=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]]) 
        maxfec=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])   
      }  
      if (fec_inter<0){
        female_fecundity=paste(female_fecundity,0.001,",",sep="")
      } else {
        female_fecundity=paste(female_fecundity,round(fec_inter*max_fec/maxfec,2),",",sep="")
      }
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
    } else {
      fec_inter=0
      female_fecundity=paste(female_fecundity,0,",",sep="")
    }
    
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output3/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output3/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output3/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}














## 4. K = 500 | Surv = Yes | Fec = Yes | Sex = / | AgeMat = / ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input4/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input4/testdiv_",species[sp],".txt",sep=""))
  }
}

lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input4",sep=""))
for (sp in 1:length(species)){
  lifetime$Maturity_F[line[sp]]=1
  lifetime$Maturity_M[line[sp]]=1
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv
  surv_grep=grep(c("Surv"),divforSlim)
  female_survival="c("
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j<(lifetime$Max_age_F[line[sp]]+2)){
      female_survival=paste(female_survival,round(1-(exp(-((((lifetime$Length_F[line[sp]+j-1])/(lifetime$Linf_F[line[sp]]))^(-1.5))*lifetime$K_F[line[sp]]))),2),",",sep="")
    } else {
      female_survival=paste(female_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv",',female_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j>=(lifetime$Maturity_F[line[sp]]+1) & j<(lifetime$Max_age_F[line[sp]]+1)){
      if (lifetime$Fecundity_length_relationship[line[sp]]=="Power"){
        fec_inter=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+j-1]^(lifetime$Beta[line[sp]]))   
        maxfec=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+lifetime$Max_age_F[line[sp]]]^(lifetime$Beta[line[sp]]))  
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Exponential") {
        fec_inter=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]])  
        maxfec=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])    
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Linear") {
        fec_inter=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]]) 
        maxfec=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])   
      }  
      if (fec_inter<0){
        female_fecundity=paste(female_fecundity,0.001,",",sep="")
      } else {
        female_fecundity=paste(female_fecundity,round(fec_inter*max_fec/maxfec,2),",",sep="")
      }
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
    } else {
      fec_inter=0
      female_fecundity=paste(female_fecundity,0,",",sep="")
    }
    
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output4/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output4/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output4/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}



## 5. K = 500 | Surv = / | Fec = / | AgeMat = Yes | Sex = / ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input5/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input5/testdiv_",species[sp],".txt",sep=""))
  }
}

lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input5",sep=""))
for (sp in 1:length(species)){
  lifetime$Maturity_F[line[sp]]=1
  lifetime$Maturity_M[line[sp]]=1
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv_female
  surv_grep=grep(c("Surv_female"),divforSlim)
  female_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[sp]]))
  
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j==(lifetime$Max_age_F[line[sp]]+2)){
      female_survival=paste(female_survival,1,")",sep="")
    } else {
      female_survival=paste(female_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_female",',female_survival,');',sep="")
  # Surv_male
  surv_grep=grep(c("Surv_male"),divforSlim)
  male_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_M[line[sp]]))
  
  for (j in 2:(lifetime$Max_age_M[line[sp]]+2)){
    if (j==(lifetime$Max_age_M[line[sp]]+2)){
      male_survival=paste(male_survival,1,")",sep="")
    } else {
      male_survival=paste(male_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_male",',male_survival,');',sep="")
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<=lifetime$Max_age_F[line[sp]] & j>lifetime$Maturity_F[line[sp]]){
      female_fecundity=paste(female_fecundity,max_fec,",",sep="")
      
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)) {
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
      
    } else {
      female_fecundity=paste(female_fecundity,0,",",sep="")
      
    }
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat_female"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Female",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #AgeMat_male
  agemat_grep=grep(c("agemat_male"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Male",',lifetime$Maturity_M[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  if (species[sp]=="Peryt" | species[sp]=="Cjuli"){
    divforSlim[starter]="p1.individuals.age = rdunif(K, min=0, max=size(Surv_female)-1);"
  } else {
    divforSlim[starter]="p1.individuals.age = rdunif(K, min=0, max=size(Surv_male)-1);"
  }
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output5/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output5/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output5/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}


















## 6. K = 500 | Surv = Yes | Fec = / | Sex = Yes | AgeMat = / ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input6/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input6/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input6",sep=""))
for (sp in 1:length(species)){
  lifetime$Maturity_F[line[sp]]=1
  lifetime$Maturity_M[line[sp]]=1
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv_female
  surv_grep=grep(c("Surv_female"),divforSlim)
  female_survival="c("
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j<(lifetime$Max_age_F[line[sp]]+2)){
      female_survival=paste(female_survival,round(1-(exp(-((((lifetime$Length_F[line[sp]+j-1])/(lifetime$Linf_F[line[sp]]))^(-1.5))*lifetime$K_F[line[sp]]))),2),",",sep="")
    } else {
      female_survival=paste(female_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_female",',female_survival,');',sep="")
  # Surv_male
  surv_grep=grep(c("Surv_male"),divforSlim)
  male_survival="c("
  for (j in 2:(lifetime$Max_age_M[line[sp]]+2)){
    if (j<(lifetime$Max_age_M[line[sp]]+2)){
      male_survival=paste(male_survival,round(1-(exp(-((((lifetime$Length_M[line[sp]+j-1])/(lifetime$Linf_M[line[sp]]))^(-1.5))*lifetime$K_M[line[sp]]))),2),",",sep="")
    } else {
      male_survival=paste(male_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_male",',male_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<=lifetime$Max_age_F[line[sp]] & j>lifetime$Maturity_F[line[sp]]){
      female_fecundity=paste(female_fecundity,max_fec,",",sep="")
      
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)) {
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
      
    } else {
      female_fecundity=paste(female_fecundity,0,",",sep="")
      
    }
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat_female"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Female",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #AgeMat_male
  agemat_grep=grep(c("agemat_male"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Male",',lifetime$Maturity_M[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  if (species[sp]=="Peryt" | species[sp]=="Cjuli"){
    divforSlim[starter]="p1.individuals.age = rdunif(K, min=0, max=size(Surv_female)-1);"
  } else {
    divforSlim[starter]="p1.individuals.age = rdunif(K, min=0, max=size(Surv_male)-1);"
  }
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output6/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output6/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output6/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}













## 7. K = 500 | Surv = / | Fec = Yes | Sex = Sex | AgeMat = / ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input7/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input7/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input7",sep=""))
for (sp in 1:length(species)){
  
  lifetime$Maturity_F[line[sp]]=1
  lifetime$Maturity_M[line[sp]]=1
  
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv_female
  surv_grep=grep(c("Surv_female"),divforSlim)
  female_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[sp]]))
  
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j==(lifetime$Max_age_F[line[sp]]+2)){
      female_survival=paste(female_survival,1,")",sep="")
    } else {
      female_survival=paste(female_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_female",',female_survival,');',sep="")
  # Surv_male
  surv_grep=grep(c("Surv_male"),divforSlim)
  male_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_M[line[sp]]))
  
  for (j in 2:(lifetime$Max_age_M[line[sp]]+2)){
    if (j==(lifetime$Max_age_M[line[sp]]+2)){
      male_survival=paste(male_survival,1,")",sep="")
    } else {
      male_survival=paste(male_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_male",',male_survival,');',sep="")
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j>=(lifetime$Maturity_F[line[sp]]+1) & j<(lifetime$Max_age_F[line[sp]]+1)){
      if (lifetime$Fecundity_length_relationship[line[sp]]=="Power"){
        fec_inter=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+j-1]^(lifetime$Beta[line[sp]]))   
        maxfec=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+lifetime$Max_age_F[line[sp]]]^(lifetime$Beta[line[sp]]))  
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Exponential") {
        fec_inter=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]])  
        maxfec=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])    
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Linear") {
        fec_inter=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]]) 
        maxfec=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])   
      }  
      if (fec_inter<0){
        female_fecundity=paste(female_fecundity,0.001,",",sep="")
      } else {
        female_fecundity=paste(female_fecundity,round(fec_inter*max_fec/maxfec,2),",",sep="")
      }
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
    } else {
      fec_inter=0
      female_fecundity=paste(female_fecundity,0,",",sep="")
    }
    
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat_female"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Female",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #AgeMat_male
  agemat_grep=grep(c("agemat_male"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Male",',lifetime$Maturity_M[line[sp]],');',sep="")
  
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  if (species[sp]=="Peryt" | species[sp]=="Cjuli"){
    divforSlim[starter]="p1.individuals.age = rdunif(K, min=0, max=size(Surv_female)-1);"
  } else {
    divforSlim[starter]="p1.individuals.age = rdunif(K, min=0, max=size(Surv_male)-1);"
  }
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output7/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output7/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output7/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}












## 8. K = 500 | Surv = Yes | Fec = Yes | Sex = Yes | AgeMat = / ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input8/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input8/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input8",sep=""))
for (sp in 1:length(species)){
  
  lifetime$Maturity_F[line[sp]]=1
  lifetime$Maturity_M[line[sp]]=1
  
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv_female
  surv_grep=grep(c("Surv_female"),divforSlim)
  female_survival="c("
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j<(lifetime$Max_age_F[line[sp]])+2){
      female_survival=paste(female_survival,round(1-(exp(-((((lifetime$Length_F[line[sp]+j-1])/(lifetime$Linf_F[line[sp]]))^(-1.5))*lifetime$K_F[line[sp]]))),2),",",sep="")
    } else {
      female_survival=paste(female_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_female",',female_survival,');',sep="")
  # Surv_male
  surv_grep=grep(c("Surv_male"),divforSlim)
  male_survival="c("
  for (j in 2:(lifetime$Max_age_M[line[sp]]+2)){
    if (j<(lifetime$Max_age_M[line[sp]]+2)){
      male_survival=paste(male_survival,round(1-(exp(-((((lifetime$Length_M[line[sp]+j-1])/(lifetime$Linf_M[line[sp]]))^(-1.5))*lifetime$K_M[line[sp]]))),2),",",sep="")
    } else {
      male_survival=paste(male_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_male",',male_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j>=(lifetime$Maturity_F[line[sp]]+1) & j<(lifetime$Max_age_F[line[sp]]+1)){
      if (lifetime$Fecundity_length_relationship[line[sp]]=="Power"){
        fec_inter=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+j-1]^(lifetime$Beta[line[sp]]))   
        maxfec=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+lifetime$Max_age_F[line[sp]]]^(lifetime$Beta[line[sp]]))  
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Exponential") {
        fec_inter=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]])  
        maxfec=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])    
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Linear") {
        fec_inter=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]]) 
        maxfec=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])   
      }  
      if (fec_inter<0){
        female_fecundity=paste(female_fecundity,0.001,",",sep="")
      } else {
        female_fecundity=paste(female_fecundity,round(fec_inter*max_fec/maxfec,2),",",sep="")
      }
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
    } else {
      fec_inter=0
      female_fecundity=paste(female_fecundity,0,",",sep="")
    }
    
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat_female"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Female",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #AgeMat_male
  agemat_grep=grep(c("agemat_male"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Male",',lifetime$Maturity_M[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  if (species[sp]=="Peryt" | species[sp]=="Cjuli"){
    divforSlim[starter]="p1.individuals.age = rdunif(K, min=0, max=size(Surv_female)-1);"
  } else {
    divforSlim[starter]="p1.individuals.age = rdunif(K, min=0, max=size(Surv_male)-1);"
  }
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output8/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output8/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output8/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}






## 9. K = 500 | Surv = / | Fec = / | Sex = / | AgeMat = Yes ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input9/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input9/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/SLiM.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input9",sep=""))
for (sp in 1:length(species)){
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv
  surv_grep=grep(c("Surv"),divforSlim)
  female_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[sp]]))
  
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_survival=paste(female_survival,1,")",sep="")
    } else {
      female_survival=paste(female_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv",',female_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<=lifetime$Max_age_F[line[sp]] & j>lifetime$Maturity_F[line[sp]]){
      female_fecundity=paste(female_fecundity,max_fec,",",sep="")
      
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)) {
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
      
    } else {
      female_fecundity=paste(female_fecundity,0,",",sep="")
      
    }
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv)-1);",sep="")
  #lifetime$Maturity_F[line[sp]]-1
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output9/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output9/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output9/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}





## 10. K = 500 | Surv = Yes | Fec = /| Sex = / | AgeMat = Yes ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input10/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input10/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/SLiM.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input10",sep=""))
for (sp in 1:length(species)){
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv
  surv_grep=grep(c("Surv"),divforSlim)
  female_survival="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<(lifetime$Max_age_F[line[sp]]+1)){
      female_survival=paste(female_survival,round(1-(exp(-((((lifetime$Length_F[line[sp]+j-1])/(lifetime$Linf_F[line[sp]]))^(-1.5))*lifetime$K_F[line[sp]]))),2),",",sep="")
    } else {
      female_survival=paste(female_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv",',female_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<=lifetime$Max_age_F[line[sp]] & j>lifetime$Maturity_F[line[sp]]){
      female_fecundity=paste(female_fecundity,max_fec,",",sep="")
      
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)) {
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
      
    } else {
      female_fecundity=paste(female_fecundity,0,",",sep="")
      
    }
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv)-1);",sep="")
  #lifetime$Maturity_F[line[sp]]-1
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output10/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output10/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output10/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}









## 11. K = 500 | Surv = / | Fec = Yes | Sex = / | AgeMat = Yes ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input11/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input11/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/SLiM.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input11",sep=""))
for (sp in 1:length(species)){
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv
  surv_grep=grep(c("Surv"),divforSlim)
  female_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[sp]]))
  
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_survival=paste(female_survival,1,")",sep="")
    } else {
      female_survival=paste(female_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv",',female_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j>=(lifetime$Maturity_F[line[sp]]+1) & j<(lifetime$Max_age_F[line[sp]]+1)){
      if (lifetime$Fecundity_length_relationship[line[sp]]=="Power"){
        fec_inter=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+j-1]^(lifetime$Beta[line[sp]]))   
        maxfec=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+lifetime$Max_age_F[line[sp]]]^(lifetime$Beta[line[sp]]))  
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Exponential") {
        fec_inter=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]])  
        maxfec=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])    
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Linear") {
        fec_inter=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]]) 
        maxfec=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])   
      }  
      female_fecundity=paste(female_fecundity,round(fec_inter*max_fec/maxfec,2),",",sep="")
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
    } else {
      fec_inter=0
      female_fecundity=paste(female_fecundity,0,",",sep="")
    }
    
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv)-1);",sep="")
  #lifetime$Maturity_F[line[sp]]-1
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output11/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output11/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output11/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}










## 12. K = 500 | Surv = Yes | Fec = Yes | Sex = / | AgeMat = Yes ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input12/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input12/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/SLiM.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input12",sep=""))
for (sp in 1:length(species)){
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv
  surv_grep=grep(c("Surv"),divforSlim)
  female_survival="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<(lifetime$Max_age_F[line[sp]]+1)){
      female_survival=paste(female_survival,round(1-(exp(-((((lifetime$Length_F[line[sp]+j-1])/(lifetime$Linf_F[line[sp]]))^(-1.5))*lifetime$K_F[line[sp]]))),2),",",sep="")
    } else {
      female_survival=paste(female_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv",',female_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j>=(lifetime$Maturity_F[line[sp]]+1) & j<(lifetime$Max_age_F[line[sp]]+1)){
      if (lifetime$Fecundity_length_relationship[line[sp]]=="Power"){
        fec_inter=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+j-1]^(lifetime$Beta[line[sp]]))   
        maxfec=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+lifetime$Max_age_F[line[sp]]]^(lifetime$Beta[line[sp]]))  
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Exponential") {
        fec_inter=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]])  
        maxfec=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])    
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Linear") {
        fec_inter=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]]) 
        maxfec=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])   
      }  
      female_fecundity=paste(female_fecundity,round(fec_inter*max_fec/maxfec,2),",",sep="")
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
    } else {
      fec_inter=0
      female_fecundity=paste(female_fecundity,0,",",sep="")
    }
    
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv)-1);",sep="")
  #lifetime$Maturity_F[line[sp]]-1
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output12/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output12/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output12/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}











## 13. K = 500 | Surv = / | Fec = / | AgeMat = Yes | Sex = Yes ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input13/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input13/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/SLiM.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input13",sep=""))
for (sp in 1:length(species)){
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv_female
  surv_grep=grep(c("Surv_female"),divforSlim)
  female_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[sp]]))
  
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_survival=paste(female_survival,1,")",sep="")
    } else {
      female_survival=paste(female_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_female",',female_survival,');',sep="")
  # Surv_male
  surv_grep=grep(c("Surv_male"),divforSlim)
  male_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_M[line[sp]]))
  
  for (j in 1:(lifetime$Max_age_M[line[sp]]+1)){
    if (j==(lifetime$Max_age_M[line[sp]]+1)){
      male_survival=paste(male_survival,1,")",sep="")
    } else {
      male_survival=paste(male_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_male",',male_survival,');',sep="")
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<=lifetime$Max_age_F[line[sp]] & j>lifetime$Maturity_F[line[sp]]){
      female_fecundity=paste(female_fecundity,max_fec,",",sep="")
      
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)) {
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
      
    } else {
      female_fecundity=paste(female_fecundity,0,",",sep="")
      
    }
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat_female"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Female",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #AgeMat_male
  agemat_grep=grep(c("agemat_male"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Male",',lifetime$Maturity_M[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  if (species[sp]=="Peryt" | species[sp]=="Cjuli"){
    divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv_female)-1);",sep="")
  } else {
    divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv_male)-1);",sep="")
  }
  #lifetime$Maturity_F[line[sp]]-1
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output13/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output13/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output13/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}














## 14. K = 500 | Surv = Yes | Fec = / | Sex = Yes | AgeMat = Yes ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input14/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input14/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/SLiM.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input14",sep=""))
for (sp in 1:length(species)){
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv_female
  surv_grep=grep(c("Surv_female"),divforSlim)
  female_survival="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<(lifetime$Max_age_F[line[sp]]+1)){
      female_survival=paste(female_survival,round(1-(exp(-((((lifetime$Length_F[line[sp]+j-1])/(lifetime$Linf_F[line[sp]]))^(-1.5))*lifetime$K_F[line[sp]]))),2),",",sep="")
    } else {
      female_survival=paste(female_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_female",',female_survival,');',sep="")
  # Surv_male
  surv_grep=grep(c("Surv_male"),divforSlim)
  male_survival="c("
  for (j in 1:(lifetime$Max_age_M[line[sp]]+1)){
    if (j<(lifetime$Max_age_M[line[sp]]+1)){
      male_survival=paste(male_survival,round(1-(exp(-((((lifetime$Length_M[line[sp]+j-1])/(lifetime$Linf_M[line[sp]]))^(-1.5))*lifetime$K_M[line[sp]]))),2),",",sep="")
    } else {
      male_survival=paste(male_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_male",',male_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j<=lifetime$Max_age_F[line[sp]] & j>lifetime$Maturity_F[line[sp]]){
      female_fecundity=paste(female_fecundity,max_fec,",",sep="")
      
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)) {
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
      
    } else {
      female_fecundity=paste(female_fecundity,0,",",sep="")
      
    }
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat_female"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Female",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #AgeMat_male
  agemat_grep=grep(c("agemat_male"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Male",',lifetime$Maturity_M[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  if (species[sp]=="Peryt" | species[sp]=="Cjuli"){
    divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv_female)-1);",sep="")
  } else {
    divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv_male)-1);",sep="")
  }
  #lifetime$Maturity_F[line[sp]]-1
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output14/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output14/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output14/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}








## 15. K = 500 | Surv = / | Fec = Yes | Sex = Sex | AgeMat = Yes ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input15/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input15/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/SLiM.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
setwd(paste(wd,"/Data/forward_slim/Input15",sep=""))
for (sp in 1:length(species)){
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv_female
  surv_grep=grep(c("Surv_female"),divforSlim)
  female_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_F[line[sp]]))
  
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_survival=paste(female_survival,1,")",sep="")
    } else {
      female_survival=paste(female_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_female",',female_survival,');',sep="")
  # Surv_male
  surv_grep=grep(c("Surv_male"),divforSlim)
  male_survival="c("
  const_surv=exp((log(0.01))/(lifetime$Max_age_M[line[sp]]))
  
  for (j in 1:(lifetime$Max_age_M[line[sp]]+1)){
    if (j==(lifetime$Max_age_M[line[sp]]+1)){
      male_survival=paste(male_survival,1,")",sep="")
    } else {
      male_survival=paste(male_survival,round(const_surv,2),",",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_male",',male_survival,');',sep="")
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j>=(lifetime$Maturity_F[line[sp]]+1) & j<(lifetime$Max_age_F[line[sp]]+1)){
      if (lifetime$Fecundity_length_relationship[line[sp]]=="Power"){
        fec_inter=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+j-1]^(lifetime$Beta[line[sp]]))   
        maxfec=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+lifetime$Max_age_F[line[sp]]]^(lifetime$Beta[line[sp]]))  
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Exponential") {
        fec_inter=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]])  
        maxfec=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])    
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Linear") {
        fec_inter=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]]) 
        maxfec=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])   
      }  
      female_fecundity=paste(female_fecundity,round(fec_inter*max_fec/maxfec,2),",",sep="")
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
    } else {
      fec_inter=0
      female_fecundity=paste(female_fecundity,0,",",sep="")
    }
    
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat_female"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Female",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #AgeMat_male
  agemat_grep=grep(c("agemat_male"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Male",',lifetime$Maturity_M[line[sp]],');',sep="")
  
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  if (species[sp]=="Peryt" | species[sp]=="Cjuli"){
    divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv_female)-1);",sep="")
  } else {
    divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv_male)-1);",sep="")
  }
  #lifetime$Maturity_F[line[sp]]-1
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output15/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output15/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output15/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}








## 16. K = 500 | Surv = Yes | Fec = Yes | Sex = Yes | AgeMat = Yes ---------
setwd(wd)
for (sp in 1:16){
  if (file.exists(paste("Data/forward_slim/Input16/testdiv_",species[sp],".txt",sep=""))==T){
    file.remove(paste("Data/forward_slim/Input16/testdiv_",species[sp],".txt",sep=""))
  }
}
lifetime<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="SLiM"))
species=lifetime[is.na(lifetime$Species_code)==FALSE,]$Species_code
latin=lifetime[is.na(lifetime$Species_code)==FALSE,]$Latin
vernacular=lifetime[is.na(lifetime$Species_code)==FALSE,]$Vernacular
line=which(is.na(lifetime$Species_code)==FALSE)
lifetime$Length=as.numeric(lifetime$Length)
lifetime$Length_F=as.numeric(lifetime$Length_F)
lifetime$Length_M=as.numeric(lifetime$Length_M)
lifetime$Max_age=as.numeric(lifetime$Max_age)
lifetime$Max_age_F=as.numeric(lifetime$Max_age_F)
lifetime$Max_age_M=as.numeric(lifetime$Max_age_M)
## CONSTRUCT SLIM MODELS
setwd(paste(wd,"/Data/forward_slim/Input16",sep=""))
for (sp in 1:length(species)){
  divforSlim<-readtext("model.txt")
  divforSlim<-strsplit(divforSlim$text,"\n")[[1]]
  # Surv_female
  surv_grep=grep(c("Surv_female"),divforSlim)
  female_survival="c("
  for (j in 2:(lifetime$Max_age_F[line[sp]]+2)){
    if (j<(lifetime$Max_age_F[line[sp]]+2)){
      female_survival=paste(female_survival,round(1-(exp(-((((lifetime$Length_F[line[sp]+j-1])/(lifetime$Linf_F[line[sp]]))^(-1.5))*lifetime$K_F[line[sp]]))),2),",",sep="")
    } else {
      female_survival=paste(female_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_female",',female_survival,');',sep="")
  # Surv_male
  surv_grep=grep(c("Surv_male"),divforSlim)
  male_survival="c("
  for (j in 2:(lifetime$Max_age_M[line[sp]]+2)){
    if (j<(lifetime$Max_age_M[line[sp]]+2)){
      male_survival=paste(male_survival,round(1-(exp(-((((lifetime$Length_M[line[sp]+j-1])/(lifetime$Linf_M[line[sp]]))^(-1.5))*lifetime$K_M[line[sp]]))),2),",",sep="")
    } else {
      male_survival=paste(male_survival,1,")",sep="")
    }
  }
  
  divforSlim[surv_grep[1]]=paste('defineConstant("Surv_male",',male_survival,');',sep="")
  
  #Fec
  fec_grep=grep(c("fec"),divforSlim)
  female_fecundity="c("
  for (j in 1:(lifetime$Max_age_F[line[sp]]+1)){
    if (j>=(lifetime$Maturity_F[line[sp]]+1) & j<(lifetime$Max_age_F[line[sp]]+1)){
      if (lifetime$Fecundity_length_relationship[line[sp]]=="Power"){
        fec_inter=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+j-1]^(lifetime$Beta[line[sp]]))   
        maxfec=lifetime$Alpha[line[sp]]*(lifetime$Length_F[line[sp]+lifetime$Max_age_F[line[sp]]]^(lifetime$Beta[line[sp]]))  
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Exponential") {
        fec_inter=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]])  
        maxfec=lifetime$Alpha[line[sp]]*exp(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])    
      } else if (lifetime$Fecundity_length_relationship[line[sp]]=="Linear") {
        fec_inter=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+j-1]*lifetime$Beta[line[sp]]) 
        maxfec=lifetime$Alpha[line[sp]]+(lifetime$Length_F[line[sp]+lifetime$Max_age[line[sp]]]*lifetime$Beta[line[sp]])   
      }  
      female_fecundity=paste(female_fecundity,round(fec_inter*max_fec/maxfec,2),",",sep="")
    } else if (j==(lifetime$Max_age_F[line[sp]]+1)){
      female_fecundity=paste(female_fecundity,max_fec,")",sep="")
    } else {
      fec_inter=0
      female_fecundity=paste(female_fecundity,0,",",sep="")
    }
    
    
  }
  divforSlim[fec_grep[1]]=paste(' defineConstant("f",',female_fecundity,');',sep="")
  
  #AgeMat_Female
  agemat_grep=grep(c("agemat_female"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Female",',lifetime$Maturity_F[line[sp]],');',sep="")
  
  #AgeMat_male
  agemat_grep=grep(c("agemat_male"),divforSlim)
  divforSlim[agemat_grep[1]]=paste('defineConstant("AgeMat_Male",',lifetime$Maturity_M[line[sp]],');',sep="")
  
  #Starter sex-specific
  starter=grep(c("starter"),divforSlim)
  if (species[sp]=="Peryt" | species[sp]=="Cjuli"){
    divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv_female)-1);",sep="")
  } else {
    divforSlim[starter]=paste("p1.individuals.age = rdunif(K, min=0, max=size(Surv_male)-1);",sep="")
  }
  #lifetime$Maturity_F[line[sp]]-1
  
  #Writefile
  writefile_grep=grep(c("writefile"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output16/div_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_popadult"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output16/popadult_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  #Writefile
  writefile_grep=grep(c("writefile_poptotal"),divforSlim)
  divforSlim[writefile_grep[1]]=paste(' writeFile("/shared/projects/abc_fish/forward_slim/Output16/poptotal_',species[sp],'_"+K+"_"+iter+".txt", line, append=T);',sep="")
  
  for (i in 1:length(divforSlim)){
    cat(divforSlim[i],file=paste("testdiv_",species[sp],".txt",sep=""),append=T,sep="\n")
  }
}

