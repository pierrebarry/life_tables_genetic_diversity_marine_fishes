#------------------------------------------------------------------------#
#                                                                        #
#               Life-history traits / Genetic diversity                  #
#                                                                        #
#------------------------------------------------------------------------#

# Load packages -----
library(lmtest)
library(ggcorrplot)
library(corrplot)
library(ggiraph)
library(Hmisc)
library(viridis)
library(readxl)
library(ggrepel)
library(ade4)
library(betareg)
library(ggthemes)
time_code=0
# Load data ----
load(file="Data/div.Rdata")
lfh<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="lfh"))
load(file="Data/div.Rdata")
lfh$div=as.vector(div)
lfh$log_fec=log(lfh$Fec)

# Correlation between life history traits and genetic diversity ----

infer_lfh=data.frame(Predictor=c(NA),
                     Dataset=c(NA),
                     pvalue=c(NA),
                     Rsquared=c(NA),
                     Estimate=c(NA),
                     sd_025=c(NA),
                     sd_0975=c(NA))

col=c(5,14,13,10,8,6,9,17,11,12)
subdata=c("Whole data set","No parental care species")
predictor=c("Body Size","Trophic level","Propagule size","Fecundity","Lifespan","Age at first maturity", "Adult lifespan","Log_Fecundity","Hermaphroditism","Parental Care")

for (j in 1:length(subdata)){
  for (i in 1:length(predictor)){

    if (i==10 & j==2){
      
      
    } else {
      
      data_plot=data.frame(species=lfh$Species,
                           y=lfh$div,
                           x=lfh[,col[i]],
                           x2=lfh$Parental_Care)
      
      if (j==1){
        m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
      } else {
        m1<-betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
      }
      infer_lfh[i+((j-1)*length(predictor)),]=c(predictor[i],
                                                subdata[j],
                                                round(coefficients(summary(m1))$mean[2,4],10),
                                                round(m1$pseudo.r.squared,3),
                                                round(coefficients(summary(m1))$mean[2,1],7),
                                                round(coefficients(summary(m1))$mean[2,1]-2*coefficients(summary(m1))$mean[2,2],7),
                                                round(coefficients(summary(m1))$mean[2,1]+2*coefficients(summary(m1))$mean[2,2],7)
      )
      
    }

  }
}

# Sensiblity ----
combi=combn(seq(1,16),11)
data_plot=data.frame(species=lfh$Species,
                     y=lfh$div,
                     x=lfh$Adult_Lifespan)
slope=c()
rsquared=c()

if (time_code==1){
  for (i in 1:ncol(combi)){
    
    data_plot_tmp=data_plot[combi[,i],]
    m1<-betareg(I(y/100)~x,data=data_plot_tmp,link='logit')
    slope=c(slope,round(coefficients(summary(m1))$mean[2,1],7))
    rsquared=c(rsquared,round(m1$pseudo.r.squared,3))
  }
  
  dd=data.frame(x=slope)
  p<-ggplot(dd,aes(x=x)) +
    geom_density(fill="grey",alpha=0.15)+
    geom_vline(xintercept=as.numeric(infer_lfh$Estimate[17]), color="red",
               linetype="dashed")+
    geom_vline(xintercept=quantile(slope,c(0.025)), color="blue")+
    geom_vline(xintercept=quantile(slope,c(0.975)), color="blue")+
    labs(x="Slope between adult lifespan and genetic diversity", 
         y = "Frequency")+
    labs(tag="A")+
    theme_classic()
  print(p)
  
  dd=data.frame(x=rsquared)
  p1<-ggplot(dd,aes(x=x)) +
    geom_density(fill="grey",alpha=0.15)+
    geom_vline(xintercept=as.numeric(infer_lfh$Rsquared[17]), color="red",
               linetype="dashed")+
    geom_vline(xintercept=quantile(rsquared,c(0.025)), color="blue")+
    geom_vline(xintercept=quantile(rsquared,c(0.975)), color="blue")+
    labs(x="Pseudo R² between adult lifespan and genetic diversity", 
         y = "Frequency")+  
    theme_classic()+
    labs(tag="B")
  print(p1)
  
  figure<-ggarrange(p,p1,ncol=1)
  pdf("figures/sensibility_lfh_slope.pdf",width=7.5,height=5)
  print(figure)
  dev.off()
}


# Interaction with hermaphroditism ----
inter_hermaphro=data.frame(Predictor=c(NA),
                     Dataset=c(NA),
                     pvalue=c(NA))

col=c(5,14,13,10,8,6,9,17)
predictor=c("Body Size","Trophic level","Propagule size","Fecundity","Lifespan","Age at first maturity", "Adult lifespan","Log_Fecundity")

for (i in 1:length(predictor)){
  
    data_plot=data.frame(species=lfh$Species,
                         y=lfh$div,
                         x=lfh[,col[i]],
                         x2=lfh$Hermaphrodism)
    
    m1<-betareg(I(y/100)~x,data=data_plot)
    m2 <- betareg(I(y/100)~x*x2,data=data_plot)
    l<-lrtest(m1,m2)
    
    inter_hermaphro[i,]=c(predictor[i],
                                              subdata[j],
                                              l$`Pr(>Chisq)`[2]
            
    )
    
}

# Interaction with parental care ----

inter_parental_care=data.frame(Predictor=c(NA),
                           Dataset=c(NA),
                           pvalue=c(NA))

col=c(5,14,13,10,8,6,9,17)
predictor=c("Body Size","Trophic level","Propagule size","Fecundity","Lifespan","Age at first maturity", "Adult lifespan","Log_Fecundity")

for (i in 1:length(predictor)){
  
  data_plot=data.frame(species=lfh$Species,
                       y=lfh$div,
                       x=lfh[,col[i]],
                       x2=lfh$Parental_Care)
  
  m1<-betareg(I(y/100)~x,data=data_plot)
  m2 <- betareg(I(y/100)~x*x2,data=data_plot)
  l<-lrtest(m1,m2)
  
  inter_parental_care[i,]=c(predictor[i],
                        subdata[j],
                        l$`Pr(>Chisq)`[2]
                        
  )
  
}

# Relation ship between log fecundity and propagule size ----
summary(lm(log(Fec)~Propagule_Size,data=lfh))
# Correlational plot ----
corr_lfh=lfh[,c("div","Body_Size","Trophic_Level","Fec","Propagule_Size","Lifespan","Age_Mat","Adult_Lifespan")]
colnames(corr_lfh)=c("Genetic diversity","Body size","Trophic level","Fecundity","Propagule size","Lifespan","Age at maturity","Adult lifespan")

cor_5 <- rcorr(as.matrix(corr_lfh))
M <- cor_5$r
p_mat <- cor_5$P

pdf(paste(wd,"/figures/corr_plot.pdf",sep=""),width=7.5,height=7.5)
corr_lfh_notebook<-corrplot.mixed(M,
                         upper="ellipse",
                         lower='number',
                         tl.pos='lt',
                         tl.col='black',
                         tl.srt=30,
                         diag='n',
                         p.mat = p_mat,
                         #addCoef.col = "black",
                         sig.level=1,
                         insig='p-value',
                         lower.col=viridis(100)[1:90],
                         upper.col=viridis(100)[1:90])
dev.off()


#corrplot(M, method = "color", col = col(200),
#         type = "upper", order = "hclust", number.cex = .7,
#         addCoef.col = "black", # Add coefficient of correlation
#         tl.col = "black", tl.srt = 90, # Text label color and rotation
#         # Combine with significance
#         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
#         # hide correlation coefficient on the principal diagonal
#         diag = FALSE)

# Plot lifespan ----
data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=lfh$Adult_Lifespan,
                     x2=lfh$Parental_Care)
m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

fit_m1_nopc_tmp=fitted(m1_nopc)
sd_025_m1_no_pc_tmp=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
sd_975_m1_no_pc_tmp=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

fit_m1_pc_tmp=fitted(m1_pc)
sd_025_m1_pc_tmp=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
sd_975_m1_pc_tmp=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

data_plot_tmp=data_plot
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$mean=fitted(m1_nopc)
data_plot_tmp$sd_025=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

data_plot_tmp2=data_plot
data_plot_tmp2=data_plot_tmp2[data_plot_tmp2$x2=="Yes",]
data_plot_tmp2$mean=fitted(m1_pc)
data_plot_tmp2$sd_025=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
data_plot_tmp2$sd_975=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

min_trait=min(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
max_trait=max(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
diff_max_min_gen=max(data_plot$y,na.rm=T)-min(data_plot$y,na.rm=T)
diff_max_min_trait=max_trait-min_trait

italic_species=c("italic('C. galerita')",
                 "italic('C. julis')",
                 "italic('D. labrax')",
                 "italic('D. puntazzo')",
                 "italic('H. guttulatus')",
                 "italic('L. budegassa')",
                 "italic('L. mormyrus')",
                 "italic('M. merluccius')",
                 "italic('M. surmuletus')",
                 "italic('P. erythrinus')",
                 "italic('S. cabrilla')",
                 "italic('S. cantharus')",
                 "italic('S. cinereus')",
                 "italic('S. pilchardus')",
                 "italic('S. sarda')",
                 "italic('S. typhle')")
col=viridis(100)
data_plot_whole=rbind(data_plot,data_plot_tmp)
data_plot_whole$gg=c(rep(1,16),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

data_plot_whole$Animate=seq(1,nrow(data_plot_whole))
data_plot_whole <- data_plot_whole[order(data_plot_whole$gg, data_plot_whole$x),]

p<-ggplot(data_plot_whole,
          aes(x=x,
              y=mean*100,
              label=species,
              group=x2))+
  geom_line(aes(x=x,y=mean*100,color=gg,group=gg),size=1)+
  geom_ribbon(data=data_plot_whole,
              aes(ymin=sd_025,
                  ymax=sd_975,
                  group=gg),
              alpha=0.1) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 shape=x2
             ))+
  geom_rangeframe()+
  theme_classic()+
  geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  scale_x_continuous(breaks=seq(1,21),
                     labels=as.character(seq(1,21)))+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                     labels=as.character(seq(0,1.5,by=0.25)))+
  xlab("Adult lifespan (years)")+
  ylab("Genetic diversity (%)")+
  scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
  scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
  theme(legend.position = c(0.84, 0.74))

x=7
pdf(paste("figures/div_Adult_lifespan.pdf",sep=""),width=1*x,height=(2/3)*x)
print(p)
dev.off()

pdf(paste("figures/div_propagule_size.pdf",sep=""),width=1*x,height=(2/3)*x)
print(p)
dev.off()
## Plot suppmat ----
data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=lfh$Body_Size,
                     x2=lfh$Parental_Care)
m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

fit_m1_nopc_tmp=fitted(m1_nopc)
sd_025_m1_no_pc_tmp=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
sd_975_m1_no_pc_tmp=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

fit_m1_pc_tmp=fitted(m1_pc)
sd_025_m1_pc_tmp=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
sd_975_m1_pc_tmp=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

data_plot_tmp=data_plot
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$mean=fitted(m1_nopc)
data_plot_tmp$sd_025=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

data_plot_tmp2=data_plot
data_plot_tmp2=data_plot_tmp2[data_plot_tmp2$x2=="Yes",]
data_plot_tmp2$mean=fitted(m1_pc)
data_plot_tmp2$sd_025=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
data_plot_tmp2$sd_975=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

min_trait=min(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
max_trait=max(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
diff_max_min_gen=max(data_plot$y,na.rm=T)-min(data_plot$y,na.rm=T)
diff_max_min_trait=max_trait-min_trait

italic_species=c("italic('C. galerita')",
                 "italic('C. julis')",
                 "italic('D. labrax')",
                 "italic('D. puntazzo')",
                 "italic('H. guttulatus')",
                 "italic('L. budegassa')",
                 "italic('L. mormyrus')",
                 "italic('M. merluccius')",
                 "italic('M. surmuletus')",
                 "italic('P. erythrinus')",
                 "italic('S. cabrilla')",
                 "italic('S. cantharus')",
                 "italic('S. cinereus')",
                 "italic('S. pilchardus')",
                 "italic('S. sarda')",
                 "italic('S. typhle')")
col=viridis(100)
data_plot_whole=rbind(data_plot,data_plot_tmp)
data_plot_whole$gg=c(rep(1,16),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

data_plot_whole$Animate=seq(1,nrow(data_plot_whole))
data_plot_whole <- data_plot_whole[order(data_plot_whole$gg, data_plot_whole$x),]

p<-ggplot(data_plot_whole,
          aes(x=x,
              y=mean*100,
              label=species,
              group=x2))+
  geom_line(aes(x=x,y=mean*100,color=gg,group=gg),size=1)+
  geom_ribbon(data=data_plot_whole,
              aes(ymin=sd_025,
                  ymax=sd_975,
                  group=gg),
              alpha=0.1) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 shape=x2
             ))+
  geom_rangeframe()+
  theme_classic()+
  geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  scale_x_continuous(breaks=seq(0,120,by=10),
                     labels=as.character(seq(0,120,by=10)))+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                     labels=as.character(seq(0,1.5,by=0.25)))+
  xlab("Body size (cm)")+
  ylab("Genetic diversity (%)")+
  scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
  scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
  theme(legend.position = c(0.84, 0.74))
p_bodysize=p

data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=lfh$Trophic_Level,
                     x2=lfh$Parental_Care)
m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

fit_m1_nopc_tmp=fitted(m1_nopc)
sd_025_m1_no_pc_tmp=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
sd_975_m1_no_pc_tmp=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

fit_m1_pc_tmp=fitted(m1_pc)
sd_025_m1_pc_tmp=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
sd_975_m1_pc_tmp=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

data_plot_tmp=data_plot
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$mean=fitted(m1_nopc)
data_plot_tmp$sd_025=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

data_plot_tmp2=data_plot
data_plot_tmp2=data_plot_tmp2[data_plot_tmp2$x2=="Yes",]
data_plot_tmp2$mean=fitted(m1_pc)
data_plot_tmp2$sd_025=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
data_plot_tmp2$sd_975=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

min_trait=min(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
max_trait=max(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
diff_max_min_gen=max(data_plot$y,na.rm=T)-min(data_plot$y,na.rm=T)
diff_max_min_trait=max_trait-min_trait

italic_species=c("italic('C. galerita')",
                 "italic('C. julis')",
                 "italic('D. labrax')",
                 "italic('D. puntazzo')",
                 "italic('H. guttulatus')",
                 "italic('L. budegassa')",
                 "italic('L. mormyrus')",
                 "italic('M. merluccius')",
                 "italic('M. surmuletus')",
                 "italic('P. erythrinus')",
                 "italic('S. cabrilla')",
                 "italic('S. cantharus')",
                 "italic('S. cinereus')",
                 "italic('S. pilchardus')",
                 "italic('S. sarda')",
                 "italic('S. typhle')")
col=viridis(100)
data_plot_whole=rbind(data_plot,data_plot_tmp)
data_plot_whole$gg=c(rep(1,16),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

data_plot_whole$Animate=seq(1,nrow(data_plot_whole))
data_plot_whole <- data_plot_whole[order(data_plot_whole$gg, data_plot_whole$x),]

p<-ggplot(data_plot_whole,
          aes(x=x,
              y=mean*100,
              label=species,
              group=x2))+
  geom_line(aes(x=x,y=mean*100,color=gg,group=gg),size=1)+
  geom_ribbon(data=data_plot_whole,
              aes(ymin=sd_025,
                  ymax=sd_975,
                  group=gg),
              alpha=0.1) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 shape=x2
             ))+
  geom_rangeframe()+
  theme_classic()+
  geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  scale_x_continuous(breaks=seq(0,5,by=0.25),
                     labels=as.character(seq(0,5,by=0.25)))+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                     labels=as.character(seq(0,1.5,by=0.25)))+
  xlab("Trophic level")+
  ylab("Genetic diversity (%)")+
  scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
  scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
  theme(legend.position = c(0.84, 0.74))
p_trophiclevel=p

data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=log(lfh$Propagule_Size),
                     x2=lfh$Parental_Care)
m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

fit_m1_nopc_tmp=fitted(m1_nopc)
sd_025_m1_no_pc_tmp=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
sd_975_m1_no_pc_tmp=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

fit_m1_pc_tmp=fitted(m1_pc)
sd_025_m1_pc_tmp=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
sd_975_m1_pc_tmp=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

data_plot_tmp=data_plot
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$mean=fitted(m1_nopc)
data_plot_tmp$sd_025=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

data_plot_tmp2=data_plot
data_plot_tmp2=data_plot_tmp2[data_plot_tmp2$x2=="Yes",]
data_plot_tmp2$mean=fitted(m1_pc)
data_plot_tmp2$sd_025=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
data_plot_tmp2$sd_975=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

min_trait=min(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
max_trait=max(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
diff_max_min_gen=max(data_plot$y,na.rm=T)-min(data_plot$y,na.rm=T)
diff_max_min_trait=max_trait-min_trait

italic_species=c("italic('C. galerita')",
                 "italic('C. julis')",
                 "italic('D. labrax')",
                 "italic('D. puntazzo')",
                 "italic('H. guttulatus')",
                 "italic('L. budegassa')",
                 "italic('L. mormyrus')",
                 "italic('M. merluccius')",
                 "italic('M. surmuletus')",
                 "italic('P. erythrinus')",
                 "italic('S. cabrilla')",
                 "italic('S. cantharus')",
                 "italic('S. cinereus')",
                 "italic('S. pilchardus')",
                 "italic('S. sarda')",
                 "italic('S. typhle')")
col=viridis(100)
data_plot_whole=rbind(data_plot,data_plot_tmp)
data_plot_whole$gg=c(rep(1,16),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

data_plot_whole$Animate=seq(1,nrow(data_plot_whole))
data_plot_whole <- data_plot_whole[order(data_plot_whole$gg, data_plot_whole$x),]

p<-ggplot(data_plot_whole,
          aes(x=x,
              y=mean*100,
              label=species,
              group=x2))+
  geom_line(aes(x=x,y=mean*100,color=gg,group=gg),size=1)+
  geom_ribbon(data=data_plot_whole,
              aes(ymin=sd_025,
                  ymax=sd_975,
                  group=gg),
              alpha=0.1) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 shape=x2
             ))+
  geom_rangeframe()+
  theme_classic()+
  geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  scale_x_continuous(breaks=seq(0,20,by=1),
                     labels=as.character(seq(0,20,by=1)))+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                     labels=as.character(seq(0,1.5,by=0.25)))+
  xlab("log[Propagule size] (mm)")+
  ylab("Genetic diversity (%)")+
  scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
  scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
  theme(legend.position = c(0.84, 0.74))
p_propagule=p

data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=lfh$Fec,
                     x2=lfh$Parental_Care)
m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=c(NA,fitted(m1))
data_plot$sd_025=c(NA,fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2])
data_plot$sd_975=c(NA,fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2])

fit_m1_nopc_tmp=fitted(m1_nopc)
sd_025_m1_no_pc_tmp=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
sd_975_m1_no_pc_tmp=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

fit_m1_pc_tmp=fitted(m1_pc)
sd_025_m1_pc_tmp=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
sd_975_m1_pc_tmp=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

data_plot_tmp=data_plot
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$mean=fitted(m1_nopc)
data_plot_tmp$sd_025=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

data_plot_tmp2=data_plot
data_plot_tmp2=data_plot_tmp2[data_plot_tmp2$x2=="Yes",]
data_plot_tmp2$mean=c(NA,fitted(m1_pc))
data_plot_tmp2$sd_025=c(NA,fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2])
data_plot_tmp2$sd_975=c(NA,fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2])

min_trait=min(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
max_trait=max(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
diff_max_min_gen=max(data_plot$y,na.rm=T)-min(data_plot$y,na.rm=T)
diff_max_min_trait=max_trait-min_trait

italic_species=c("italic('C. galerita')",
                 "italic('C. julis')",
                 "italic('D. labrax')",
                 "italic('D. puntazzo')",
                 "italic('H. guttulatus')",
                 "italic('L. budegassa')",
                 "italic('L. mormyrus')",
                 "italic('M. merluccius')",
                 "italic('M. surmuletus')",
                 "italic('P. erythrinus')",
                 "italic('S. cabrilla')",
                 "italic('S. cantharus')",
                 "italic('S. cinereus')",
                 "italic('S. pilchardus')",
                 "italic('S. sarda')",
                 "italic('S. typhle')")
col=viridis(100)
data_plot_whole=rbind(data_plot,data_plot_tmp)
data_plot_whole$gg=c(rep(1,16),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

data_plot_whole$Animate=seq(1,nrow(data_plot_whole))
data_plot_whole <- data_plot_whole[order(data_plot_whole$gg, data_plot_whole$x),]

p<-ggplot(data_plot_whole,
          aes(x=x,
              y=mean*100,
              label=species,
              group=x2))+
  geom_line(aes(x=x,y=mean*100,color=gg,group=gg),size=1)+
  geom_ribbon(data=data_plot_whole,
              aes(ymin=sd_025,
                  ymax=sd_975,
                  group=gg),
              alpha=0.1) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 shape=x2
             ))+
  geom_rangeframe()+
  theme_classic()+
  geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  scale_x_continuous(breaks=seq(0,50000,by=5000),
                     labels=as.character(seq(0,50000,by=5000)))+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                     labels=as.character(seq(0,1.5,by=0.25)))+
  xlab("Fecundity (number of eggs)")+
  ylab("Genetic diversity (%)")+
  scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
  scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
  theme(legend.position = c(0.84, 0.74))
p_fec=p

data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=lfh$Lifespan,
                     x2=lfh$Parental_Care)
m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

fit_m1_nopc_tmp=fitted(m1_nopc)
sd_025_m1_no_pc_tmp=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
sd_975_m1_no_pc_tmp=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

fit_m1_pc_tmp=fitted(m1_pc)
sd_025_m1_pc_tmp=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
sd_975_m1_pc_tmp=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

data_plot_tmp=data_plot
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$mean=fitted(m1_nopc)
data_plot_tmp$sd_025=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]

data_plot_tmp2=data_plot
data_plot_tmp2=data_plot_tmp2[data_plot_tmp2$x2=="Yes",]
data_plot_tmp2$mean=fitted(m1_pc)
data_plot_tmp2$sd_025=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
data_plot_tmp2$sd_975=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]

min_trait=min(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
max_trait=max(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
diff_max_min_gen=max(data_plot$y,na.rm=T)-min(data_plot$y,na.rm=T)
diff_max_min_trait=max_trait-min_trait

italic_species=c("italic('C. galerita')",
                 "italic('C. julis')",
                 "italic('D. labrax')",
                 "italic('D. puntazzo')",
                 "italic('H. guttulatus')",
                 "italic('L. budegassa')",
                 "italic('L. mormyrus')",
                 "italic('M. merluccius')",
                 "italic('M. surmuletus')",
                 "italic('P. erythrinus')",
                 "italic('S. cabrilla')",
                 "italic('S. cantharus')",
                 "italic('S. cinereus')",
                 "italic('S. pilchardus')",
                 "italic('S. sarda')",
                 "italic('S. typhle')")
col=viridis(100)
data_plot_whole=rbind(data_plot,data_plot_tmp)
data_plot_whole$gg=c(rep(1,16),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

data_plot_whole$Animate=seq(1,nrow(data_plot_whole))
data_plot_whole <- data_plot_whole[order(data_plot_whole$gg, data_plot_whole$x),]

p<-ggplot(data_plot_whole,
          aes(x=x,
              y=mean*100,
              label=species,
              group=x2))+
  geom_line(aes(x=x,y=mean*100,color=gg,group=gg),size=1)+
  geom_ribbon(data=data_plot_whole,
              aes(ymin=sd_025,
                  ymax=sd_975,
                  group=gg),
              alpha=0.1) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 shape=x2
             ))+
  geom_rangeframe()+
  theme_classic()+
  geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  scale_x_continuous(breaks=seq(0,25,by=1),
                     labels=as.character(seq(0,25,by=1)))+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                     labels=as.character(seq(0,1.5,by=0.25)))+
  xlab("Lifespan (years)")+
  ylab("Genetic diversity (%)")+
  scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
  scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
  theme(legend.position = c(0.84, 0.74))
p_lifespan=p

pdf(paste("figures/div_suppmat.pdf",sep=""),width=15,height=17.5)
print(ggarrange(p_bodysize,p_trophiclevel,p_fec,p_propagule,p_lifespan,
                ncol=2,nrow=3,
                labels=c("A","B","C","D","E")))
dev.off()
## Plot for notebook ----

legend_notebook=c("Body size (cm)",
                  "Age at maturity (years)",
                  "Length at maturity (cm)",
                  "Lifespan (years)",
                  "Adult lifespan (years)",
                  "Propagule size (mm)",
                  "Trophic level")

p_lfh_notebook=vector('list',7)
a=0
for (i in c(5,6,7,8,9,13,14)){
  a=a+1
  message(a)
  p_lfh_notebook[[a]] <- local({
    i <- i
    data_plot=data.frame(species=lfh$Species_plot,
                         y=lfh$div,
                         x=lfh[,i],
                         x2=lfh$Parental_Care)
    wikipedia=c("Coryphoblennius_galerita",
                "Coris_julis",
                "Dicentrarchus_labrax",
                "Diplodus_puntazzo",
                "Hippocampus_guttulatus",
                "Lophius_budegassa",
                "Lithognathus_mormyrus",
                "Merluccius_merluccius",
                "Mullus_surmuletus",
                "Pagellus_erythrinus",
                "Serranus_cabrilla",
                "Spondyliosoma_cantharus",
                "Symphodus_cinereus",
                "Sardina_pilchardus",
                "Sarda_sarda",
                "Syngnathus_typhle")
    
    data_plot$onclick <- sprintf("window.open(\"%s%s\")",
                              "http://en.wikipedia.org/wiki/", as.character(wikipedia) )
    m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
    m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
    m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')
    
    data_plot$mean=fitted(m1)
    data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
    data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]
    
    fit_m1_nopc_tmp=fitted(m1_nopc)
    sd_025_m1_no_pc_tmp=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
    sd_975_m1_no_pc_tmp=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]
    
    fit_m1_pc_tmp=fitted(m1_pc)
    sd_025_m1_pc_tmp=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
    sd_975_m1_pc_tmp=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]
    
    data_plot_tmp=data_plot
    data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
    data_plot_tmp$mean=fitted(m1_nopc)
    data_plot_tmp$sd_025=fitted(m1_nopc)*100-2*coefficients(summary(m1_nopc))$mean[2,2]
    data_plot_tmp$sd_975=fitted(m1_nopc)*100+2*coefficients(summary(m1_nopc))$mean[2,2]
    
    data_plot_tmp2=data_plot
    data_plot_tmp2=data_plot_tmp2[data_plot_tmp2$x2=="Yes",]
    data_plot_tmp2$mean=fitted(m1_pc)
    data_plot_tmp2$sd_025=fitted(m1_pc)*100-2*coefficients(summary(m1_pc))$mean[2,2]
    data_plot_tmp2$sd_975=fitted(m1_pc)*100+2*coefficients(summary(m1_pc))$mean[2,2]
    
    min_trait=min(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
    max_trait=max(data_plot$x[is.na(data_plot$y)==FALSE],na.rm=T)
    diff_max_min_gen=max(data_plot$y,na.rm=T)-min(data_plot$y,na.rm=T)
    diff_max_min_trait=max_trait-min_trait
    
    italic_species=c("italic('C. galerita')",
                     "italic('C. julis')",
                     "italic('D. labrax')",
                     "italic('D. puntazzo')",
                     "italic('H. guttulatus')",
                     "italic('L. budegassa')",
                     "italic('L. mormyrus')",
                     "italic('M. merluccius')",
                     "italic('M. surmuletus')",
                     "italic('P. erythrinus')",
                     "italic('S. cabrilla')",
                     "italic('S. cantharus')",
                     "italic('S. cinereus')",
                     "italic('S. pilchardus')",
                     "italic('S. sarda')",
                     "italic('S. typhle')")
    col=viridis(100)
    p<-ggplot(data_plot, 
              aes(x = x, 
                  y = y,
                  label=species,
                  group=x2)
              
    )+
      geom_line(data=data_plot_tmp,aes(x=x,
                                       y=mean*100,
                                       color=col[75]),
                size=2)+
      geom_line(aes(x=x,
                    y=data_plot$mean*100,
                    color=col[62]),
                size=2,
                linetype='dashed')+
      #geom_line(data=data_plot_tmp2,aes(x=x,
      #                                  y=mean*100,
      #                                  color=col[50]),
      #          size=2,
      #          linetype='dashed')+
      geom_ribbon(aes(ymin=data_plot$sd_025,
                      ymax=data_plot$sd_975),
                  #fill="grey",
                  fill=col[25],
                  alpha=0.1,
                  col="white") +
      #geom_ribbon(data=data_plot_tmp2,aes(ymin=sd_025,
      #                                    ymax=sd_975),
      #            fill="grey",
      #            #fill=col[75],
      #            alpha=0.1,
      #            col="white") +
      geom_ribbon(data=data_plot_tmp,aes(ymin=sd_025,
                                         ymax=sd_975),
                  #fill="grey",
                  fill=col[50],
                  alpha=0.1,
                  col="white") +
      #geom_point(size=2.5,
      #           aes(shape=x2)) +
      geom_point_interactive(
        aes(tooltip = paste("Heterozygosity:",round(y,3),
                            "\n Vernacular:",lfh$Vernacular,
                            "\n Body size (cm):",lfh$Body_Size,
                            "\n Trophic level:",lfh$Trophic_Level,
                            "\n Fecundity:",lfh$Fec,
                            "\n Propagule size (cm):",lfh$Propagule_Size,
                            "\n Age at maturity (years):",lfh$Age_Mat,
                            "\n Length at maturity (cm):",lfh$L_Mat,
                            "\n Lifespan (years):",lfh$Lifespan,
                            "\n Adult lifspan (years):",lfh$Adult_Lifespan,
                            "\n Hermaphroditism :",lfh$Hermaphrodism,
                            "\n Parental care:",x2),
            onclick=onclick,
            shape=x2),
        size = 2.5,
        #col=geno$color,
        position='jitter',
        alpha=0.75)+
      #geom_rangeframe()+
      theme_classic()+
      geom_text_repel(label=italic_species,col="black",parse=T)+
      #scale_x_continuous(breaks=seq(1,21),
      #                   labels=as.character(seq(1,21)))+
      scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                         labels=as.character(seq(0,1.5,by=0.25)))+
      xlab(legend_notebook[a])+
      ylab("Heterozygosity (%)")+
      scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
      scale_shape_manual(name = "Brooding behaviour",values=c(19, 18))+
      theme(legend.position = "bottom")
    #print(ex)
  })
  
}

## Check beta -- OPTIONAL ----
sapply(c("probit", "cloglog", "cauchit", "loglog"),
       function(x) lrtest(m1,update(m1, link = x))$`Pr(>Chisq)`[2])

sapply(c("log", "sqrt"),
       function(x) lrtest(m1,update(m1, link.phi = x))$`Pr(>Chisq)`[2])