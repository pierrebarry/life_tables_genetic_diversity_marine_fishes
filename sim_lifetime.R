#------------------------------------------------------------------------#
#                                                                        #
#               Run life tables simulations with AgeNe                   #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages ----
library(emdbook)
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
range01 <- function(x,ff){((1-ff)*((x-min(x))/(max(x)-min(x))))+ff}
#https://www.symbolab.com/solver/non-linear-system-of-equations-calculator/solve%20for%20b%2C%20ln%5Cleft(0.01%5Cright)%20%2B%5Cleft(%5Cfrac%7B20%7D%7Bb%7D%5Cright)%5E%7B1.4%7D%3D0

# Preliminary : examples of age-survival and fecundity ----
# Age-survival model
c=lseq(0.1,30,length.out=50)
L=10
plot_true=1
for (j in 1:length(c)){
  x=seq(0,L)
  b=-L/(-4.60517^(1/c[j]))
  Sx=c(1,rep(NA,length(x)-1))
  for (i in 2:length(x)){
    Sx[i]=exp(-((x[i]/b)^(c[j])))
  }

  Mx=c(1,rep(NA,length(x)-1))
  for (i in 2:length(x)){
    Mx[i]=1-(1-exp((x[i-1]/b)^(c[j])-(x[i]/b)^(c[j])))
  }
  
  if (plot_true==1){
    if (j==1){
      plot(x[-1],Mx[-1],type="l",ylim=c(0,1),main="Age-specific mortality",
           xlab="Age",ylab="Probability")
    } else {
      lines(x[-1],Mx[-1],type="l")
    }
    
    print(paste("Probability of survival until age ",x,":",Sx,sep=""))
    readline(c[j])
  }
  

}

for (j in 1:length(c)){
  x=seq(0,L)
  b=-L/(-4.60517^(1/c[j]))
  Sx=c(1,rep(NA,length(x)-1))
  for (i in 2:length(x)){
    Sx[i]=exp(-((x[i]/b)^(c[j])))
  }
  
  Mx=c(1,rep(NA,length(x)-1))
  for (i in 2:length(x)){
    Mx[i]=1-(1-exp((x[i-1]/b)^(c[j])-(x[i]/b)^(c[j])))
  }
  
  if (plot_true==1){
    
    if (j==1){
      plot(x[-1],Sx[-1],type="l",ylim=c(0,1),
           main="Cumulative survival",xlab="Age",ylab="Probability")
    } else {
      lines(x[-1],Sx[-1],type="l")
    }
    print(paste("Probability of survival until age ",x,":",Sx,sep=""))
    readline(c[j])
  }
  
  
}

# Age-fecundity model
Lifespan=c(10)
i=1
Maturity=c(1)
# Linear : F = a*L + b
f=seq(0,100,length.out=50)
for (u in 1:50){
  
  age=seq(1,Lifespan[i])
  fec=rep(0,Lifespan[i])
  for (j in (Maturity[i]):(Lifespan[i])){
    fec[j]=f[u]*(j-Maturity[i])+100
  }

  fec[Maturity[i]:Lifespan[i]]=(fec[Maturity[i]:Lifespan[i]]*f[u]/max(fec[Maturity[i]:Lifespan[i]]))
  
  if (u==1){
    plot(age,fec,type="l",ylim=c(0,100))
  } else {
    lines(fec)
  }
  readline(f[u])
}

# Exponential : f = a*exp(L*b)
f=seq(-1,1,length.out=50)
for (u in 1:length(f)){
  
  age=seq(1,Lifespan[i])
  fec=rep(0,Lifespan[i])
  for (j in Maturity[i]:Lifespan[i]){
    fec[j]=1*exp((j-Maturity[i]+1)*f[u])
  }
  fec[Maturity[i]:Lifespan[i]]=fec[Maturity[i]:Lifespan[i]]/max(fec[Maturity[i]:Lifespan[i]])
  
  if (u==1){
    plot(age,fec,type="l")
  } else {
    lines(fec)
  }
  readline(f[u])
}

# Power-law : f = a*L^(b)
f=c(seq(-5,0,length.out=25),seq(0,5,length.out=25))

for (u in 1:length(f)){
  
  age=seq(1,Lifespan[i])
  fec=rep(0,Lifespan[i])
  for (j in Maturity[i]:Lifespan[i]){
    fec[j]=1*((j-Maturity[i]+1)^(f[u]))
  }
  
  if (u==1){
    plot(age,fec,type="l",ylim=c(0,10000))
  } else {
    lines(fec)
  }
  
  readline(f[u])
  
}


# Polynomial function
f=seq(-1,1,length.out=50)
for (u in 1:length(f)){
  cc=1
  d=1
  age=seq(1,Lifespan[i])
  fec=rep(0,Lifespan[i])
  for (j in 1:Lifespan[i]){
    fec[j]=cc*(j-Maturity[i])*((Maturity[i]+Lifespan[i]-j)^(2))
  }

  fec[Maturity[i]:Lifespan[i]]=fec[Maturity[i]:Lifespan[i]]+abs(min(fec[Maturity[i]:Lifespan[i]]))
  fec[Maturity[i]:Lifespan[i]]=fec[Maturity[i]:Lifespan[i]]/max(fec[Maturity[i]:Lifespan[i]])
  if (Maturity[i]>=2){
    for (j in 1:(Maturity[i]-1)){
      fec[j]=0
    }
  }
  range01 <- function(x,ff){((1-ff)*((x-min(x))/(max(x)-min(x))))+ff}
  
  fec[Maturity[i]:Lifespan[i]]=range01(fec[Maturity[i]:Lifespan[i]],f[u])
  
  if (u==1){
    plot(age,fec,type="l",ylim=c(0,1))
  } else {
    lines(fec)
  }
  readline(f[u])
  
}

# Built simulated life tables and run AgeNe ----
setwd("Data/sim_lifetime")
lifetime<-as.data.frame(read_excel("C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT/Data/agene/agene.xlsx",sheet="RawData"))
Lifespan=lifetime$Max_age[is.na(lifetime$Max_age)==FALSE]
Maturity=lifetime$Maturity[is.na(lifetime$Maturity)==FALSE]
red=0.01

# c values
c=lseq(0.1,30,length.out=50)

for (fecundity_model in seq(1,5)){
  remember=1
  par(mfrow=c(1,1))
  
  if (fecundity_model==1){
    f=c(1)
  }
  
  if (fecundity_model==2){
    
    f=seq(-1,1,length.out=50)
    
    for (i in 1:length(Lifespan)){
      for (u in 1:length(f)){
        
        age=seq(1,Lifespan[i])
        fec=rep(0,Lifespan[i])
        if (f[u]<0){
          for (j in 1:Lifespan[i]){
            fec[j]=((-f[u]-1)/Lifespan[i])*age[j]+1 
          }
        } else {
          for (j in 1:Lifespan[i]){
            fec[j]=((1-f[u])/Lifespan[i])*age[j]+f[u] 
          }        
        }
        
        
        if (Maturity[i]>=2){
          for (j in 1:(Maturity[i]-1)){
            fec[j]=0
          }
        }
        
        
        if (u==1){
          plot(c(0,age),
               c(0,fec),
               type="l",
               ylim=c(0,1),
               xlim=c(0,Lifespan[i]))
        } else {
          lines(c(0,age),
                c(0,fec))
        }
        
      }
      
    }
    
  }
  
  if (fecundity_model==3){
    f=seq(-1,1,length.out=50)
    
    for (i in 1:length(Lifespan)){
      for (u in 1:length(f)){
        
        age=seq(1,Lifespan[i])
        fec=rep(0,Lifespan[i])
        for (j in Maturity[i]:Lifespan[i]){
          fec[j]=1*exp((j-Maturity[i]+1)*f[u])
        }
        
        if (u==1){
          plot(age,fec,type="l",ylim=c(0,1))
        } else {
          lines(fec)
        }
      }
    }
    
  }
  
  if (fecundity_model==4){
    f=seq(-5,5,length.out=50)
    
    for (i in 1:length(Lifespan)){
      for (u in 1:length(f)){
        
        age=seq(1,Lifespan[i])
        fec=rep(0,Lifespan[i])
        for (j in Maturity[i]:Lifespan[i]){
          fec[j]=1*((j-Maturity[i]+1)^(f[u]))
        }
        
        if (u==1){
          plot(age,fec,type="l",ylim=c(0,1))
        } else {
          lines(fec)
        }
      }
    }
  }
  
  if (fecundity_model==5){
    
    f=seq(-1,1,length.out=50)
    
    for (i in 1:length(Lifespan)){
      for (u in 1:length(f)){
        
        cc=1
        d=1
        age=seq(1,Lifespan[i])
        fec=rep(0,Lifespan[i])
        for (j in 1:Lifespan[i]){
          fec[j]=cc*(j-Maturity[i])*((Maturity[i]+Lifespan[i]-j)^(2))
        }
        
        
        if (Maturity[i]>=2){
          for (j in 1:(Maturity[i]-1)){
            fec[j]=0
          }
        }
        
        fec=fec/(max(fec))
        
        if (u==1){
          plot(age,fec,type="l",ylim=c(0,1))
        } else {
          lines(fec)
        }
      }
    }
  }
  
  species=paste("Sp",seq(1,16),sep="")
  sim_lifetime=vector('list',5)
  for (k in 1:length(sim_lifetime)){
    sim_lifetime[[k]]=data.frame(Species=paste("Sp",seq(1,16),sep=""))
  }
  
  
  for (k in 1:length(sim_lifetime)){
    for (u in 1:(length(c)*length(f))){
      sim_lifetime[[k]]=cbind(sim_lifetime[[k]],rep(0,length(species)))
    }
  }
  
  col=c()
  for (j in 1:length(f)){
    for (i in 1:length(c)){
      col=c(col,paste("f",f[j],"c",c[i],sep=""))
    }
  }
  
  for (k in 1:length(sim_lifetime)){
    colnames(sim_lifetime[[k]])=c("Species",col)
  }
  
  for (i in 1:length(species)){
    for (k in 1:length(c)){
      
      x=seq(0,Lifespan[i])
      b=-Lifespan[i]/(-4.60517^(1/c[k]))
      Sx=c(1,rep(NA,length(x)-1))
      for (j in 2:length(x)){
        Sx[j]=exp(-((x[j]/b)^(c[k])))
      }
      
      Mx=c(1,rep(NA,length(x)-1))
      for (j in 2:length(x)){
        Mx[j]=1-(1-exp((x[j-1]/b)^(c[k])-(x[j]/b)^(c[k])))
      }
      Mx=Mx
      x=x
    }
    
  }
  
  for (u in 1:length(f)){
    for (k in 1:length(c)){
      file.remove("cogediv_sensibility.txt",showWarnings = FALSE)
      species=paste("Sp",seq(1,16),sep="")
      for (i in 1:length(species)){
        
        x=seq(0,Lifespan[i])
        b=-Lifespan[i]/(-4.60517^(1/c[k]))
        Sx=c(1,rep(NA,length(x)-1))
        for (j in 2:length(x)){
          Sx[j]=exp(-((x[j]/b)^(c[k])))
        }
        
        Mx=c(1,rep(NA,length(x)-1))
        for (j in 2:length(x)){
          Mx[j]=1-(1-exp((x[j-1]/b)^(c[k])-(x[j]/b)^(c[k])))
        }
        Mx=Mx
        
        
        # f = constant
        if (fecundity_model==1){
          fec=rep(0,Lifespan[i])
          for (j in Maturity[i]:Lifespan[i]){
            fec[j]=1
          }
          
          if (Maturity[i]>=2){
            for (j in 1:(Maturity[i]-1)){
              fec[j]=0
            }
          }
        }
        
        ## Compute f if linear
        if (fecundity_model==2){
          
          
          age=seq(1,Lifespan[i])
          fec=rep(0,Lifespan[i])
          if (f[u]<0){
            for (j in 1:Lifespan[i]){
              fec[j]=((-f[u]-1)/Lifespan[i])*age[j]+1 
            }
          } else {
            for (j in 1:Lifespan[i]){
              fec[j]=((1-f[u])/Lifespan[i])*age[j]+f[u] 
            }        
          }
          
          
          if (Maturity[i]>=2){
            for (j in 1:(Maturity[i]-1)){
              fec[j]=0
            }
          }
          
        } 
        
        if (fecundity_model==3){
          
          age=seq(1,Lifespan[i])
          fec=rep(0,Lifespan[i])
          for (j in Maturity[i]:Lifespan[i]){
            fec[j]=1*exp((j-Maturity[i]+1)*f[u])
          }
          
        }
        
        if (fecundity_model==4){
          
          age=seq(1,Lifespan[i])
          fec=rep(0,Lifespan[i])
          for (j in Maturity[i]:Lifespan[i]){
            fec[j]=1*((j-Maturity[i]+1)^(f[u]))
          }
          
        }
        
        if (fecundity_model==5){
          
          cc=1
          d=1
          age=seq(1,Lifespan[i])
          fec=rep(0,Lifespan[i])
          for (j in 1:Lifespan[i]){
            fec[j]=cc*(j-Maturity[i])*((Maturity[i]+Lifespan[i]-j)^(2))
          }
          
          
          if (Maturity[i]>=2){
            for (j in 1:(Maturity[i]-1)){
              fec[j]=0
            }
          }
          
          fec=fec/(max(fec))
          
        }
        
        #First line
        firstline=paste(species[i]," ; c = ",c[k]," ; f = ",1,sep="")
        if (file.exists("cogediv_sensibility.txt")==FALSE){
          write(firstline,file="cogediv_sensibility.txt",append="FALSE")
        } else {
          write(firstline,file="cogediv_sensibility.txt",append="TRUE")
        }
        
        # Second line
        secondline=c(Lifespan[i],
                     "1000",
                     0.5)
        write(secondline,file="cogediv_sensibility.txt",append="TRUE",ncolumns = 3)
        
        # Other line
        for (j in 1:Lifespan[i]){
          
          if (j<Maturity[i]){
            fec[j]=0
          }
          
          if (j==Lifespan[i]){
            Mx[j]=0
          }
          otherline=c(j,
                      round(Mx[j],2),
                      round(fec[j],2),
                      1,
                      round(Mx[j],2),
                      round(fec[j],2),
                      1)
          
          write(otherline,file="cogediv_sensibility.txt",append="TRUE",ncolumns=7)
          
        }
        
      }
      
      system( shQuote( "AgeNe.exe"),input=c("cogediv_sensibility.txt","output_sensibility.txt"))
      
      ## EXTRACT INFORMATIONS ----
      agene<-readtext("output_sensibility.txt")
      agene<-strsplit(agene$text,"\n")[[1]]
      
      position=c()
      for (j in 1:length(species)){
        position[j]=grep(species[j],agene)[1]
      }
      position=c(position,length(agene))
      
      remember=remember+1
      for (i in 1:length(species)){
        
        inf=position[i]
        sup=position[i+1]
        
        if (length(intersect(grep("overall",agene[inf:sup]),grep("Vk",agene[inf:sup])))>0){
          
          intersect=intersect(grep("overall",agene[inf:sup]),grep("Vk",agene[inf:sup]))
          sim_lifetime[[1]][i,remember]=as.numeric(numextract(agene[inf:sup][intersect]))
          intersect=intersect(grep("female",agene[inf:sup]),grep("Vk",agene[inf:sup]))
          sim_lifetime[[2]][i,remember]=as.numeric(numextract(agene[inf:sup][intersect]))
          intersect=intersect+1
          sim_lifetime[[3]][i,remember]=as.numeric(numextract(agene[inf:sup][intersect]))
          intersect=intersect(grep("Ne",agene[inf:sup]),grep("Total N",agene[inf:sup]))
          sim_lifetime[[4]][i,remember]=as.numeric(Numextract(agene[inf:sup][intersect]))[2]
          intersect=intersect(grep("Ne",agene[inf:sup]),grep("Adult N",agene[inf:sup]))
          sim_lifetime[[5]][i,remember]=as.numeric(Numextract(agene[inf:sup][intersect]))[2]
          
        } else {
          
          intersect=grep("Vk",agene[inf:sup])
          sim_lifetime[[1]][i,remember]=as.numeric(numextract(agene[inf:sup][intersect]))
          intersect=intersect(grep("Ne",agene[inf:sup]),grep("Total N",agene[inf:sup]))
          sim_lifetime[[4]][i,remember]=as.numeric(Numextract(agene[inf:sup][intersect]))
          intersect=intersect(grep("Ne",agene[inf:sup]),grep("Adult N",agene[inf:sup]))
          sim_lifetime[[5]][i,remember]=as.numeric(Numextract(agene[inf:sup][intersect]))
          
        }
        
      }
      
    }
  }
  
  # Constant
  if (fecundity_model==1){
    simulation_lifetime_constant=sim_lifetime
    save(simulation_lifetime_constant,file=paste("simulation_lifetime_constant_",red,".Rdata",sep=""))
  }
  
  # Linear
  if (fecundity_model==2){
    simulation_lifetime_linear=sim_lifetime
    save(simulation_lifetime_linear,file=paste("simulation_lifetime_linear_",red,".Rdata",sep=""))
  }
  
  # Exponential
  if (fecundity_model==3){
    simulation_lifetime_exp=sim_lifetime
    save(simulation_lifetime_exp,file=paste("simulation_lifetime_exp_",red,"_noscale.Rdata",sep=""))
  }
  
  # Power
  if (fecundity_model==4){
    simulation_lifetime_power=sim_lifetime
    save(simulation_lifetime_power,file=paste("simulation_lifetime_power_",red,"_noscale.Rdata",sep=""))
  }
  
  # Polynomial
  if (fecundity_model==5){
    simulation_lifetime_poly=sim_lifetime
    save(simulation_lifetime_poly,file=paste("simulation_lifetime_poly_",red,"_first.Rdata",sep=""))
  }
  
}



