#------------------------------------------------------------------------#
#                                                                        #
#                     Analyse simulated life tables                      #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages ----
library(emdbook)
library(ggplot2)
library(ggpubr)
library(readxl)
library(ggsci)
library(wesanderson)
# Load functions ----
Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
} 
# Plot Figure 4 - Example ----
# Survival
c=lseq(0.1,30,length.out=25)
L=10
data_prelm=data.frame(Age=c(NA),
                      c=c(NA),
                      Sx=c(NA))
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
  
  data_dd=cbind(seq(0,L),
                rep(c[j],L+1),
                Sx)
  colnames(data_dd)=colnames(data_prelm)
  data_prelm=rbind(data_prelm,
                   data_dd)
}
data_prelm=data_prelm[-1,]

size=c(rep(1,11),
       rep(0.1,11*11),
       rep(1,11),
       rep(0.1,11*11),
       rep(1,11))


cc=wes_palette("Zissou1", 3, type = "continuous")
col=c(rep(cc[1],11),
      rep("gray",11*11),
      rep(cc[2],11),
      rep("gray",11*11),
      rep(cc[3],11))

data_prelm$c=factor(data_prelm$c)
sim_1<-ggplot(data_prelm,aes(x=Age,y=Sx))+
  geom_line(aes(fill=c),size=size,col=col)+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("")+
  scale_x_continuous(breaks=seq(0,10,by=1),
                     labels=as.character(seq(0,10,by=1)))+
  xlab("")+
  ylab("Probability of survival")+
  annotate("text", x = 0.25, y = 0.15, label = "Type III",size=3,col=cc[1])+
  annotate("text", x = 4.25, y = 0.5, label = "Type II",size=3,col=cc[2])+
  annotate("text", x = 9.25, y = 0.95, label = "Type I",size=3,col=cc[3])

# Fecundity
f=seq(0,1,length.out=25)
L=10
data_prelm=data.frame(Age=c(NA),
                      f=c(NA),
                      Fec=c(NA))
for (j in 1:length(f)){
  x=seq(0,L)

  fec=c()
  for (i in 0:L){
    fec[i+1]=1*exp(f[j]*i)
  }
  
  
  if (f[j]<0){
    fec=fec-1
  }
  
  data_dd=cbind(seq(0,L),
                rep(f[j],L+1),
                fec)
  colnames(data_dd)=colnames(data_prelm)
  data_prelm=rbind(data_prelm,
                   data_dd)
}
data_prelm=data_prelm[-1,]

data_prelm$f=factor(data_prelm$f)

size=c(rep(0.1,11*12),
       rep(0.1,11),
       rep(0.1,4*11),
       rep(0.1,11),
       rep(0.1,11*7))

sim_2.1<-ggplot(data_prelm,aes(x=Age,y=Fec))+
  geom_line(aes(fill=f),color="gray")+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("")+
  scale_y_log10()+
  scale_x_continuous(breaks=seq(0,10,by=1),
                     labels=as.character(seq(0,10,by=1)))+
  ylab("")

f=seq(-1,-0.04166667,length.out=24)
L=10
data_prelm=data.frame(Age=c(NA),
                      f=c(NA),
                      Fec=c(NA))
for (j in 1:length(f)){
  x=seq(0,L)
  
  fec=c()
  for (i in 0:L){
    fec[i+1]=1*exp(f[j]*i)
  }
  
  fec=fec/(max(fec))
  
  data_dd=cbind(seq(0,L),
                rep(f[j],L+1),
                fec)
  colnames(data_dd)=colnames(data_prelm)
  data_prelm=rbind(data_prelm,
                   data_dd)
}
data_prelm=data_prelm[-1,]

data_prelm$f=factor(data_prelm$f)

size=c(rep(0.1,11*12),
       rep(0.1,11),
       rep(0.1,4*11),
       rep(0.1,11),
       rep(0.1,11*7))

sim_2.2<-ggplot(data_prelm,aes(x=Age,y=Fec))+
  geom_line(aes(fill=f),color="gray",size=0.1)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_log10()+
  xlab("")+
  scale_x_continuous(breaks=seq(0,10,by=1),
                     labels=as.character(seq(0,10,by=1)))+
  xlab("Age")+
  ylab("")

sim_2<-ggarrange(sim_2.1,sim_2.2,nrow=2)
sim_2<-annotate_figure(sim_2,
                    left=text_grob(expression("Age-specific fecundity"),rot=90,vjust=2),
)

pres<-ggarrange(sim_1,sim_2,nrow=2)

# Analysis ----
lifetime<-as.data.frame(read_excel("Data/agene/agene.xlsx",sheet="RawData"))
Lifespan=lifetime$Max_age[is.na(lifetime$Max_age)==FALSE]
Maturity=lifetime$Maturity[is.na(lifetime$Maturity)==FALSE]
Adult_Lifespan=Lifespan-Maturity
species=paste("Sp",seq(1,16),sep="")
load(file="Data/species_c_f.Rdata")
species_c_f$c=as.numeric(species_c_f$c)
species_c_f$beta=as.numeric(species_c_f$beta)
pc=c(1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1)

breaks_slope=seq(-0.08,0.015,by=0.005)
breaks_pvalue=c(0,0.05,1)
breaks_rsquared=rev(seq(0,1,by=0.05))

## Fec = cst ----
load(file="Data/sim_lifetime/simulation_lifetime_constant_0.01.Rdata")
c=lseq(0.1,30,length.out=50)
f=c(1)
n=length(f)*length(c)
data_sim=data.frame(f=rep(0,n),
                    c=rep(0,n),
                    slope=rep(0,n),
                    pvalue=rep(0,n),
                    rsquared=rep(0,n))

for (i in 2:length(colnames(simulation_lifetime_constant[[4]]))){
  if (substr(colnames(simulation_lifetime_constant[[4]])[i],2,2)=="-"){
    data_sim$f[i-1]=-as.numeric(Numextract(colnames(simulation_lifetime_constant[[4]])[i])[1]) 
  } else {
    data_sim$f[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_constant[[4]])[i])[1])
  }
  data_sim$c[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_constant[[4]])[i])[2])
  data_plot=data.frame(species=species,
                       y=simulation_lifetime_constant[[4]][,i],
                       x=Lifespan)
  m1<-lm(y~x,data=data_plot[pc==0,])
  data_sim$slope[i-1]=summary(m1)$coefficients[2,1]
  data_sim$pvalue[i-1]=summary(m1)$coefficients[2,4]
  data_sim$rsquared[i-1]=summary(m1)$r.squared
  
}

load(file="Data/sim_lifetime/simulation_lifetime_exp_0.01_noscale.Rdata")
c=lseq(0.1,30,length.out=50)
f=seq(-1,1,length.out=50)
n=length(f)*length(c)

data_sim_exp=data.frame(f=rep(0,n),
                        c=rep(0,n),
                        slope=rep(0,n),
                        pvalue=rep(0,n),
                        rsquared=rep(0,n))

for (i in 2:length(colnames(simulation_lifetime_exp[[4]]))){
  if (substr(colnames(simulation_lifetime_exp[[4]])[i],2,2)=="-"){
    data_sim_exp$f[i-1]=-as.numeric(Numextract(colnames(simulation_lifetime_exp[[4]])[i])[1]) 
  } else {
    data_sim_exp$f[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_exp[[4]])[i])[1])
  }
  data_sim_exp$c[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_exp[[4]])[i])[2])
  data_plot=data.frame(species=species,
                       y=simulation_lifetime_exp[[4]][,i],
                       x=Adult_Lifespan)
  m1<-lm(y~x,data=data_plot)
  data_sim_exp$slope[i-1]=summary(m1)$coefficients[2,1]
  data_sim_exp$pvalue[i-1]=summary(m1)$coefficients[2,4]
  data_sim_exp$rsquared[i-1]=summary(m1)$r.squared
  
}

data_sim$fecundity=rep("Constant",nrow(data_sim))
data_sim_exp$fecundity=rep("Exponential",nrow(data_sim_exp))
data_sim=rbind(data_sim,data_sim_exp[data_sim_exp$f=="0.142857142857143",])

data_sim$Animate=c(rep(0,50),seq(1,50))

animate=0
if (animate==1){
  p1<-ggplot(NULL)+
    geom_line(data=data_sim[data_sim$fecundity=="Constant",],aes(x=c,y=slope,size=fecundity),col=rep(wes_palette("Zissou1", 50, type = "continuous"),1))+
    theme_classic()+
    scale_x_log10()+
    ylim(c(-0.075,0.0075))+
    xlab("c")+
    ylab("Slope adult lifespan ~ Ne/N")+
    scale_size_manual(name="Fecundity model",values=c(0.5,1))+
    theme(legend.position = "bottom")
  #transition_reveal(Animate)
  
  p1<-ggplot(NULL)+
    geom_line(data=data_sim,aes(x=c,y=slope,size=fecundity),col=rep(wes_palette("Zissou1", 50, type = "continuous"),2))+
    theme_classic()+
    scale_x_log10()+
    ylim(c(-0.075,0.0075))+
    xlab("")+
    ylab("Slope adult lifespan ~ Ne/N")+
    scale_size_manual(name="Fecundity model",values=c(0.5,1))+
    theme(legend.position = "top")
  # transition_reveal(Animate)
  
  animate(p1, fps = 20, renderer = gifski_renderer(loop = FALSE),
          height = 4, width = 6, units = "in", res = 300)
  
  anim_save("C:/Users/ordinateur/ownCloud/COGEDIV/PRESENTATION/JDD_2020/RMARKDOWN/lifetime_2.gif")
  
} else {
  p1<-ggplot(NULL)+
    geom_line(data=data_sim,aes(x=c,y=slope,size=fecundity),col=rep(wes_palette("Zissou1", 50, type = "continuous"),2))+
    scale_x_log10()+
    xlab("")+
    ylab("Slope adult lifespan ~ Ne/N")+
    scale_size_manual(name="Fecundity model",values=c(0.5,1))+
    #transition_reveal(c)+
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.position = "bottom"
    )
}

p2<-ggplot(NULL)+
  geom_line(data=data_sim,aes(x=c,y=rsquared,size=fecundity),col=rep(wes_palette("Zissou1", 50, type = "continuous"),2))+  
  theme_classic()+
  scale_x_log10()+
  xlab("")+
  ylab("R²")+
  scale_size_manual(name="Fecundity model",values=c(0.5,1))+
  theme(legend.position = "bottom")
p3<-ggplot(NULL)+
  geom_line(data=data_sim,aes(x=c,y=pvalue),col=rep(wes_palette("Zissou1", 50, type = "continuous"),2))+
  theme_classic()+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(yintercept=0.05,col="red")

p4<-ggarrange(p1,p2,nrow=2,common.legend = T,legend="top")

p4<-annotate_figure(p4,
                bottom=text_grob(expression("Age-specific survival variable constant ("*italic(c)*")"),vjust=-1),
)

sim_top=ggarrange(pres,p4,ncol=2,labels=c("A","B"))

## Fec = exponential ----
load(file="Data/sim_lifetime/simulation_lifetime_exp_0.01_noscale.Rdata")
c=lseq(0.1,30,length.out=50)
f=seq(-1,1,length.out=50)
n=length(f)*length(c)

data_sim_exp=data.frame(f=rep(0,n),
                    c=rep(0,n),
                    slope=rep(0,n),
                    pvalue=rep(0,n),
                    rsquared=rep(0,n))

for (i in 2:length(colnames(simulation_lifetime_exp[[4]]))){
  if (substr(colnames(simulation_lifetime_exp[[4]])[i],2,2)=="-"){
    data_sim_exp$f[i-1]=-as.numeric(Numextract(colnames(simulation_lifetime_exp[[4]])[i])[1]) 
  } else {
    data_sim_exp$f[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_exp[[4]])[i])[1])
  }
  data_sim_exp$c[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_exp[[4]])[i])[2])
  data_plot=data.frame(species=species,
                       y=simulation_lifetime_exp[[4]][,i],
                       x=Adult_Lifespan)
  m1<-lm(y~x,data=data_plot)
  data_sim_exp$slope[i-1]=summary(m1)$coefficients[2,1]
  data_sim_exp$pvalue[i-1]=summary(m1)$coefficients[2,4]
  data_sim_exp$rsquared[i-1]=summary(m1)$r.squared
  
}

p3 <- ggplot(NULL) +
  geom_contour_filled(data=data_sim_exp,aes(c,f,z=slope),breaks=breaks_slope)+
  geom_contour(data=data_sim_exp,aes(c,f,z=pvalue),breaks=breaks_pvalue,col="red",linetype='dashed',size=0.5)+
  scale_x_log10()+
  scale_color_jco()+
  theme_classic()+
  xlab(expression("Age-specific survival variable ("*italic(c)*")"))+
  ylab(expression("Age-specific fecundity variable ("*italic(f)*")"))+
  ggtitle("Exponential")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="Slope")+
  theme(legend.position="left")
p3.rsq <- ggplot(NULL) +
  geom_contour_filled(data=data_sim_exp,aes(c,f,z=rsquared),breaks=breaks_rsquared)+
  geom_contour(data=data_sim_exp,aes(c,f,z=pvalue),breaks=breaks_pvalue,col="red",linetype='dashed',size=0.5)+
  scale_x_log10()+
  theme_classic()+
  xlab(expression("Age-specific survival variable ("*italic(c)*")"))+
  ylab(expression("Age-specific fecundity variable ("*italic(f)*")"))+
  ggtitle("Exponential")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="R²")+
  theme(legend.position="right")

sim_bottom=ggarrange(p3,p3.rsq,ncol=2,labels=c("C","D"),label.x=-0.0075)
sim_fig=ggarrange(sim_top,sim_bottom,nrow=2)

pdf(paste(wd,"/figures/sim_lifetime_new.pdf",sep=""),width=15,height=10)
print(sim_fig)
dev.off()

pdf(paste(wd,"/figures/sim_top.pdf",sep=""),width=10,height=5)
print(sim_top)
dev.off()

pdf(paste(wd,"/figures/sim_bottom.pdf",sep=""),width=14,height=5.5)
print(sim_bottom)
dev.off()

## Fec = linear ----
load(file="Data/sim_lifetime/simulation_lifetime_linear_0.01.Rdata")
c=lseq(0.1,30,length.out=50)
f=seq(-1,1,length.out=50)
n=length(f)*length(c)
data_sim=data.frame(f=rep(0,n),
                    c=rep(0,n),
                    slope=rep(0,n),
                    pvalue=rep(0,n),
                    rsquared=rep(0,n))

for (i in 2:length(colnames(simulation_lifetime_linear[[4]]))){
  if (substr(colnames(simulation_lifetime_linear[[4]])[i],2,2)=="-"){
    data_sim$f[i-1]=-as.numeric(Numextract(colnames(simulation_lifetime_linear[[4]])[i])[1]) 
  } else {
    data_sim$f[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_linear[[4]])[i])[1])
  }
  data_sim$c[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_linear[[4]])[i])[2])
  data_plot=data.frame(species=species,
                       y=simulation_lifetime_linear[[4]][,i],
                       x=Adult_Lifespan)
  m1<-lm(y~x,data=data_plot)
  data_sim$slope[i-1]=summary(m1)$coefficients[2,1]
  data_sim$pvalue[i-1]=summary(m1)$coefficients[2,4]
  data_sim$rsquared[i-1]=summary(m1)$r.squared
  
}

p3 <- ggplot(NULL) +
  geom_contour_filled(data=data_sim,aes(c,f,z=slope),breaks=breaks_slope)+
  geom_contour(data=data_sim,aes(c,f,z=pvalue),breaks=breaks_pvalue,col="red",linetype='dashed',size=0.5)+
  scale_x_log10()+
  scale_color_jco()+
  theme_classic()+
  xlab(expression("Age-specific survival variable ("*italic(c)*")"))+
  ylab(expression("Age-specific fecundity variable ("*italic(f)*")"))+
  ggtitle("Linear")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="Slope")+
  theme(legend.position="left")

p3.rsq <- ggplot(NULL) +
  geom_contour_filled(data=data_sim,aes(c,f,z=rsquared),breaks=breaks_rsquared)+
  geom_contour(data=data_sim,aes(c,f,z=pvalue),breaks=breaks_pvalue,col="red",linetype='dashed',size=0.5)+
  scale_x_log10()+
  theme_classic()+
  xlab(expression("Age-specific survival variable ("*italic(c)*")"))+
  ylab(expression("Age-specific fecundity variable ("*italic(f)*")"))+
  ggtitle("Linear")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="R²")+
  theme(legend.position="right")


sim_bottom_linear=ggarrange(p3,p3.rsq,ncol=2)

## Fec = polynomual ----
load(file="Data/sim_lifetime/simulation_lifetime_poly_0.01_first.Rdata")
c=lseq(0.25,30,length.out=50)
f=seq(-1,1,length.out=50)
n=length(f)*length(c)
data_sim=data.frame(f=rep(0,n),
                    c=rep(0,n),
                    slope=rep(0,n),
                    pvalue=rep(0,n),
                    rsquared=rep(0,n))

for (i in 2:length(colnames(simulation_lifetime_poly[[4]]))){
  if (substr(colnames(simulation_lifetime_poly[[4]])[i],2,2)=="-"){
    data_sim$f[i-1]=-as.numeric(Numextract(colnames(simulation_lifetime_poly[[4]])[i])[1]) 
  } else {
    data_sim$f[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_poly[[4]])[i])[1])
  }
  data_sim$c[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_poly[[4]])[i])[2])
  data_plot=data.frame(species=species,
                       y=simulation_lifetime_poly[[4]][,i],
                       x=Adult_Lifespan)
  m1<-lm(y~x,data=data_plot)
  data_sim$slope[i-1]=summary(m1)$coefficients[2,1]
  data_sim$pvalue[i-1]=summary(m1)$coefficients[2,4]
  data_sim$rsquared[i-1]=summary(m1)$r.squared
  
}

p3 <- ggplot(NULL) +
  geom_contour_filled(data=data_sim,aes(c,f,z=slope),breaks=breaks_slope)+
  geom_contour(data=data_sim,aes(c,f,z=pvalue),breaks=breaks_pvalue,col="red",linetype='dashed',size=0.5)+
  scale_x_log10()+
  scale_color_jco()+
  theme_classic()+
  xlab(expression("Age-specific survival variable ("*italic(c)*")"))+
  ylab(expression("Age-specific fecundity variable ("*italic(f)*")"))+
  ggtitle("Polynomial")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="Slope")+
  theme(legend.position="left")

p3.rsq <- ggplot(NULL) +
  geom_contour_filled(data=data_sim,aes(c,f,z=rsquared),breaks=breaks_rsquared)+
  geom_contour(data=data_sim,aes(c,f,z=pvalue),breaks=breaks_pvalue,col="red",linetype='dashed',size=0.5)+
  scale_x_log10()+
  theme_classic()+
  xlab(expression("Age-specific survival variable ("*italic(c)*")"))+
  ylab(expression("Age-specific fecundity variable ("*italic(f)*")"))+
  ggtitle("Polynomial")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="R²")+
  theme(legend.position="right")


sim_bottom_poly=ggarrange(p3,p3.rsq,ncol=2)


## Fec = Power ----
load(file="Data/sim_lifetime/simulation_lifetime_power_0.01_noscale.Rdata")
c=lseq(0.1,30,length.out=50)
f=seq(-5,5,length.out=50)

n=length(f)*length(c)
data_sim=data.frame(f=rep(0,n),c=rep(0,n),slope=rep(0,n),pvalue=rep(0,n),
                    rsquared=rep(0,n))
for (i in 2:length(colnames(simulation_lifetime_power[[4]]))){
  if (substr(colnames(simulation_lifetime_power[[4]])[i],2,2)=="-"){
    data_sim$f[i-1]=-as.numeric(Numextract(colnames(simulation_lifetime_power[[4]])[i])[1]) 
  } else {
    data_sim$f[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_power[[4]])[i])[1])
  }
  data_sim$c[i-1]=as.numeric(Numextract(colnames(simulation_lifetime_power[[4]])[i])[2])
  data_plot=data.frame(species=species,
                       y=simulation_lifetime_power[[4]][,i],
                       x=Lifespan)
  m1<-lm(y~x,data=data_plot)
  data_sim$slope[i-1]=summary(m1)$coefficients[2,1]
  data_sim$pvalue[i-1]=summary(m1)$coefficients[2,4]
  data_sim$rsquared[i-1]=summary(m1)$r.squared
  
}

p3 <- ggplot(NULL) +
  geom_contour_filled(data=data_sim,aes(c,f,z=slope),breaks=breaks_slope)+
  geom_contour(data=data_sim,aes(c,f,z=pvalue),breaks=breaks_pvalue,col="red",linetype='dashed',size=0.5)+
  scale_x_log10()+
  scale_color_jco()+
  theme_classic()+
  xlab(expression("Age-specific survival variable ("*italic(c)*")"))+
  ylab(expression("Age-specific fecundity variable ("*italic(f)*")"))+
  ggtitle("Power")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="Slope")+
  theme(legend.position="left")

p3.rsq <- ggplot(NULL) +
  geom_contour_filled(data=data_sim,aes(c,f,z=rsquared),breaks=breaks_rsquared)+
  geom_contour(data=data_sim,aes(c,f,z=pvalue),breaks=breaks_pvalue,col="red",linetype='dashed',size=0.5)+
  scale_x_log10()+
  theme_classic()+
  xlab(expression("Age-specific survival variable ("*italic(c)*")"))+
  ylab(expression("Age-specific fecundity variable ("*italic(f)*")"))+
  ggtitle("Power")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="R²")+
  theme(legend.position="right")


sim_bottom_power=ggarrange(p3,p3.rsq,ncol=2)

## Fig Suppmat alternative fecundity model ----
pdf(paste(wd,"/figures/sim_lifetime_suppmat.pdf",sep=""),width=20,height=20)
print(ggarrange(sim_bottom_linear,
                sim_bottom_poly,
                sim_bottom_power,
                nrow=3))
dev.off()