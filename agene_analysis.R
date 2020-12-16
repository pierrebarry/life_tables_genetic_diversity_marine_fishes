#------------------------------------------------------------------------#
#                                                                        #
#                        AgeNe output analysis                           #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages ----
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(betareg)
library(readxl)
library(ggiraph)
library(ggrepel)
library(viridis)
time_code=0

# Load function -----
test_slope_different_1=function(fit=fit){
  # get from : https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r
  # Compute Summary with statistics      
  sfit<- summary(fit)
  # Compute t-student H0: intercept=5. The estimation of coefficients and their s.d. are in sfit$coefficients
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  return(pval)
}

test_intercept_different_0=function(fit=fit){
  # get from : https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r
  sfit<- summary(fit)
  # Compute t-student H0: intercept=5. The estimation of coefficients and their s.d. are in sfit$coefficients
  tstats <- (0-sfit$coefficients[1,1])/sfit$coefficients[1,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  return(pval)
}

# Load data ----
load(file="Data/div.Rdata")
load(file="Data/agene/agene_output.Rdata")
lfh<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="lfh"))
lfh$div=as.vector(div)
# Analysis ----
for (i in 1:5){
  agene_output[[i]]$body_size=lfh$Body_Size
}

agene_test=data.frame(Sim=seq(1,16),AgeMat=rep(c(rep(0,1),rep(1,1)),8) ,
                      Surv=rep(c(rep(0,2),rep(1,2)),4),
                      Fec=rep(c(rep(0,4),rep(1,4)),2),
                      Sex=rep(c(rep(0,8),rep(1,8))),
                      test_slope=rep(0,16),test_intercept=rep(0,16),
                      slope_mean=rep(0,16),slope_sd=rep(0,16),
                      inter_mean=rep(0,16),inter_sd=rep(0,16),
                      test_slope_nopc=rep(0,16),test_intercept_nopc=rep(0,16),
                      slope_mean_nopc=rep(0,16),slope_sd_nopc=rep(0,16),
                      inter_mean_nopc=rep(0,16),inter_sd_nopc=rep(0,16),
                      pvalue=rep(0,16),pvalue_nopc=rep(0,16)
                      )

sensibility_agene_slope=vector('list',16)
sensibility_agene_intercept=vector('list',16)
slope=vector('list',16)
intercept=vector('list',16)
pvalue=vector('list',16)


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

agene_output$onclick <- sprintf("window.open(\"%s%s\")",
                             "http://en.wikipedia.org/wiki/", as.character(wikipedia) )

p_agene=vector('list',16)
p_agene_pairwise=vector('list',16)
p_agene_pairwise_nopc=vector('list',16)

for (sim in 1:16){
  
  title=c("--",
          "AgeMat",
          "Surv",
          "AgeMat + Surv",
          "Fec",
          "AgeMat + Fec",
          "Surv + Fec",
          "AgeMat + Surv + Fec",
          "Sex",
          "AgeMat + Sex",
          "Surv + Sex",
          "AgeMat + Surv + Sex",
          "Fec + Sex",
          "AgeMat + Fec + Sex",
          "Surv + Fec + Sex",
          "AgeMat + Surv + Fec + Sex")
  
  p_agene[[sim]]<-local({
    data_plot=data.frame(species=lfh$Species_plot,
                         y=lfh$div,
                         x=agene_output[[4]][,sim+3],
                         x2=lfh$Parental_Care)
    data_plot$y=data_plot$y/(max(data_plot$y))
    data_plot$x=data_plot$x/(max(data_plot$x))
    
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
    
    div_agene<-ggplot(data_plot_whole, 
                      aes(x = x, 
                          y = mean*100,
                          label=species,
                      )) +
      geom_line(aes(colour=gg),
                size=1)+
      geom_point(data=data_plot,size=2.5,
                 aes(x=x,
                     y=y,
                     shape=x2))+
      geom_rangeframe()+
      theme_classic()+
      geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
      scale_x_continuous(breaks=seq(0,1,by=0.1),
                         labels=as.character(seq(0,1,by=0.1)),
                         limits = c(0,1))+
      scale_y_continuous(breaks=seq(0,1,by=0.1),
                         labels=as.character(seq(0,1,by=0.1)),
                         limits=c(0,1))+
      xlab(expression(alpha))+
      ylab("")+
      scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
      scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))+
      geom_abline(intercept=0,slope=1,col="red",lty=2)+
      annotate("text", x = 0.2025, y = 0.9, label = paste("Slope = ",round(coefficients(summary(m1_nopc))$mean[2,1],2)," (sd = ",round(coefficients(summary(m1_nopc))$mean[2,2],2),")",sep=""))+
      ggtitle(title[sim])
    
  })
  
  compare_het=data.frame(SP1=c(NA),
                         SP2=c(NA),
                         ratio_mean=c(NA),
                         ratio_truehet=c(NA),
                         cross=c(NA))
  
  Sp=c("Lbude","Hgutt","Dlabr","Scant","Dpunt","Lmorm","Cgale","Scine","Mmerl","Ssard","Styph","Peryt",
       "Msurm","Cjuli","Scabr","Spilc")
  
  for (j in 1:(nrow(agene_output[[1]])-1)){
    for (k in (j+1):nrow(agene_output[[1]])){
      vect=c(as.character(Sp[j]),
             as.character(Sp[k]),
             round(agene_output[[4]][which(agene_output[[4]]$Species_code==Sp[j]),sim+3]/agene_output[[4]][which(agene_output[[4]]$Species_code==Sp[k]),sim+3],3),
             round(div[which(names(div)==Sp[j])]/div[which(names(div)==Sp[k])],3)
              )
      compare_het=rbind(compare_het,vect)
    }
  }
  
  for (j in 1:(nrow(agene_output[[1]]))){
    for (k in 1:nrow(agene_output[[1]])){
      if (j!=k){
        vect=c(as.character(Sp[j]),
               as.character(Sp[k]),
               round(agene_output[[4]][which(agene_output[[4]]$Species_code==Sp[j]),sim+3]/agene_output[[4]][which(agene_output[[4]]$Species_code==Sp[k]),sim+3],3),
               round(div[which(names(div)==Sp[j])]/div[which(names(div)==Sp[k])],3)
        )
        compare_het=rbind(compare_het,vect)
      }

    }
  }
  
  compare_het=compare_het[-1,]
  compare_het$ratio_mean=as.numeric(compare_het$ratio_mean)
  compare_het$ratio_truehet=as.numeric(compare_het$ratio_truehet)
  
  col=c()
  for (i in 1:nrow(compare_het)){
    
    if (compare_het$SP1[i]=="Hgutt" |
        compare_het$SP2[i]=="Hgutt"){
      
      col[i]="red"
      
    } else {
      col[i]="blue"
    }
    
  }
  
  compare_het$col=col

  # Parental care
  compare_het_stat=compare_het
  fit<-lm(compare_het_stat$ratio_mean~compare_het_stat$ratio_truehet)
  
  p1<-ggplot(compare_het_stat,aes(ratio_truehet,ratio_mean))+
    theme_classic()+
    geom_rangeframe()+
    geom_point(size=2)+
    geom_smooth(method="lm",fullrange="T")+
    geom_abline(intercept=0,slope=1,col="red")+
    xlab("Ratio of estimated NeN")+
    ylab("Ratio of observed genetic diversity")+
    ylim(c(0,max(compare_het_stat$ratio_mean)))+
    xlim(c(0,1))+
    ggtitle("Whole dataset")
  
  p_agene_pairwise[[sim]] <- local({
    p1_notebook<-ggplot(compare_het_stat,aes(ratio_truehet,ratio_mean))+
      theme_classic()+
      geom_rangeframe()+
      geom_smooth(method="lm",fullrange="T")+
      geom_abline(intercept=0,slope=1,col="red")+
      geom_point_interactive(aes(tooltip = paste("Species 1:",compare_het_stat$SP1,
                                                 "\n Species 2:",compare_het_stat$SP2,
                                                 "\n Ratio observed:",compare_het_stat$ratio_truehet,
                                                 "\n Ratio estimated:",compare_het_stat$ratio_mean)),
                             size=2)+
      xlab("Ratio of estimated NeN")+
      ylab("Ratio of observed genetic diversity")+
      ylim(c(0,max(compare_het_stat$ratio_mean)))+
      xlim(c(0,1))+
      ggtitle("Whole dataset")
    
  })

  
  agene_test$test_slope[sim]=test_slope_different_1(fit=fit)
  agene_test$test_intercept[sim]=test_intercept_different_0(fit=fit)
  agene_test$slope_mean[sim]=coefficients(summary(fit))[2,1]
  agene_test$slope_sd[sim]=coefficients(summary(fit))[2,2]
  agene_test$inter_mean[sim]=coefficients(summary(fit))[1,1]
  agene_test$inter_sd[sim]=coefficients(summary(fit))[1,2] 
  agene_test$pvalue[sim]=coefficients(summary(fit))[2,4] 
  
  # No parental care
  compare_het_stat_nopc=compare_het
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Hgutt",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Hgutt",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Scine",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Scine",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Styph",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Styph",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Cgale",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Cgale",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Scant",]
  compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Scant",]
  
  fit<-lm(compare_het_stat_nopc$ratio_mean~compare_het_stat_nopc$ratio_truehet)
  
  p2<-ggplot(compare_het_stat_nopc,aes(ratio_truehet,ratio_mean))+
    theme_classic()+
    geom_rangeframe()+
    geom_point(size=2)+
    geom_smooth(data=compare_het_stat_nopc,method="lm",fullrange="T")+
    geom_abline(intercept=0,slope=1,col="red")+
    xlab("Ratio of estimated NeN")+
    ylab("Ratio of observed genetic diversity")+
    ylim(c(0,max(compare_het_stat_nopc$ratio_mean)))+
    xlim(c(0,1))+
    ggtitle("No parental care species")
  
  p_agene_pairwise_nopc[[sim]] <- local({
    p2_notebook<-ggplot(compare_het_stat_nopc,aes(ratio_truehet,ratio_mean))+
      theme_classic()+
      geom_rangeframe()+
      geom_smooth(method="lm",fullrange="T")+
      geom_abline(intercept=0,slope=1,col="red")+
      geom_point_interactive(aes(tooltip = paste("Species 1:",compare_het_stat_nopc$SP1,
                                                 "\n Species 2:",compare_het_stat_nopc$SP2,
                                                 "\n Ratio observed:",compare_het_stat_nopc$ratio_truehet,
                                                 "\n Ratio estimated:",compare_het_stat_nopc$ratio_mean)),
                             size=2)+
      xlab("Ratio of estimated NeN")+
      ylab("Ratio of observed genetic diversity")+
      ylim(c(0,max(compare_het_stat_nopc$ratio_mean)))+
      xlim(c(0,1))+
      ggtitle("Only non brooding species")
    
  })
  
  agene_test$test_slope_nopc[sim]=test_slope_different_1(fit=fit)
  agene_test$test_intercept_nopc[sim]=test_intercept_different_0(fit=fit)
  agene_test$slope_mean_nopc[sim]=coefficients(summary(fit))[2,1]
  agene_test$slope_sd_nopc[sim]=coefficients(summary(fit))[2,2]
  agene_test$inter_mean_nopc[sim]=coefficients(summary(fit))[1,1]
  agene_test$inter_sd_nopc[sim]=coefficients(summary(fit))[1,2] 
  agene_test$pvalue_nopc[sim]=coefficients(summary(fit))[2,4] 
  
  
}

# Plot suppmat figure ----
pdf(paste(wd,"/figures/agene_suppmat.pdf",sep=""),width=30,height=30)
print(ggarrange(p_agene[[1]],
                p_agene[[2]],
                p_agene[[3]],
                p_agene[[4]],
                p_agene[[5]],
                p_agene[[6]],
                p_agene[[7]],
                p_agene[[8]],
                p_agene[[9]],
                p_agene[[10]],
                p_agene[[11]],
                p_agene[[12]],
                p_agene[[13]],
                p_agene[[14]],
                p_agene[[15]],
                p_agene[[16]],
                nrow=4,
                ncol=4))
dev.off()

figure1<-annotate_figure(ggarrange(p_agene_pairwise[[1]], p_agene_pairwise_nopc[[1]], ncol = 2, nrow = 1),
                         top=title[1])
figure2<-annotate_figure(ggarrange(p_agene_pairwise[[2]], p_agene_pairwise_nopc[[2]], ncol = 2, nrow = 1),
                         top=title[2])
figure3<-annotate_figure(ggarrange(p_agene_pairwise[[3]], p_agene_pairwise_nopc[[3]], ncol = 2, nrow = 1),
                         top=title[3])
figure4<-annotate_figure(ggarrange(p_agene_pairwise[[4]], p_agene_pairwise_nopc[[4]], ncol = 2, nrow = 1),
                         top=title[4])
figure5<-annotate_figure(ggarrange(p_agene_pairwise[[5]], p_agene_pairwise_nopc[[5]], ncol = 2, nrow = 1),
                         top=title[5])
figure6<-annotate_figure(ggarrange(p_agene_pairwise[[6]], p_agene_pairwise_nopc[[6]], ncol = 2, nrow = 1),
                        top=title[6])
figure7<-annotate_figure(ggarrange(p_agene_pairwise[[7]], p_agene_pairwise_nopc[[7]], ncol = 2, nrow = 1),
                         top=title[7])
figure8<-annotate_figure(ggarrange(p_agene_pairwise[[8]], p_agene_pairwise_nopc[[8]], ncol = 2, nrow = 1),
                         top=title[8])
figure9<-annotate_figure(ggarrange(p_agene_pairwise[[9]], p_agene_pairwise_nopc[[9]], ncol = 2, nrow = 1),
                         top=title[9])
figure10<-annotate_figure(ggarrange(p_agene_pairwise[[10]], p_agene_pairwise_nopc[[10]], ncol = 2, nrow = 1),
                         top=title[10])
figure11<-annotate_figure(ggarrange(p_agene_pairwise[[11]], p_agene_pairwise_nopc[[11]], ncol = 2, nrow = 1),
                         top=title[11])
figure12<-annotate_figure(ggarrange(p_agene_pairwise[[12]], p_agene_pairwise_nopc[[12]], ncol = 2, nrow = 1),
                         top=title[12])
figure13<-annotate_figure(ggarrange(p_agene_pairwise[[13]], p_agene_pairwise_nopc[[13]], ncol = 2, nrow = 1),
                         top=title[13])
figure14<-annotate_figure(ggarrange(p_agene_pairwise[[14]], p_agene_pairwise_nopc[[14]], ncol = 2, nrow = 1),
                         top=title[14])
figure15<-annotate_figure(ggarrange(p_agene_pairwise[[15]], p_agene_pairwise_nopc[[15]], ncol = 2, nrow = 1),
                         top=title[15])
figure16<-annotate_figure(ggarrange(p_agene_pairwise[[16]], p_agene_pairwise_nopc[[16]], ncol = 2, nrow = 1),
                         top=title[16])
figure1<-annotate_figure(ggarrange(p_agene_pairwise[[1]], p_agene_pairwise_nopc[[1]], ncol = 2, nrow = 1),
                         top=title[1])
pdf(paste(wd,"/figures/agene_pairwise_suppmat.pdf",sep=""),width=30,height=30)
print(ggarrange(figure1,
      figure2,
      figure3,
      figure4,
      figure5,
      figure6,
      figure7,
      figure8,
      figure9,
      figure10,
      figure11,
      figure12,
      figure13,
      figure14,
      figure15,
      figure16,
      nrow=4,
      ncol=4))
dev.off()

# All possible combinations
combi=combn(seq(1,16),11)
slope=c()
pvalue=c()
rsquared=c()
if (time_code==1){
  for (i in 1:ncol(combi)){
    aaa<-betareg(I(lfh[combi[,i],]$div/100)~agene_output[[4]][combi[,i],]$Output16)
    slope=c(slope,coefficients(summary(aaa))$mean[2,1])
    pvalue=c(pvalue,coefficients(summary(aaa))$mean[2,4])
    rsquared=c(rsquared,aaa$pseudo.r.squared)
  }
  
  sim=16
  data_plot=data.frame(species=lfh$Species_plot,
                       y=lfh$div,
                       x=agene_output[[4]][,sim+3],
                       x2=lfh$Parental_Care)
  m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')

  dd=data.frame(x=slope)
  p<-ggplot(dd,aes(x=x)) +
    geom_density(fill="grey",alpha=0.15)+
    geom_vline(xintercept=as.numeric(coefficients(m1_nopc)[2]), color="red",
               linetype="dashed")+
    geom_vline(xintercept=quantile(slope,c(0.025)), color="blue")+
    geom_vline(xintercept=quantile(slope,c(0.975)), color="blue")+
    labs(x="Slope between Ne/N and genetic diversity", 
         y = "Frequency")+
    labs(tag="A")+
    theme_classic()
  print(p)
  
  dd=data.frame(x=rsquared)
  p1<-ggplot(dd,aes(x=x)) +
    geom_density(fill="grey",alpha=0.15)+
    geom_vline(xintercept=0.55, color="red",
               linetype="dashed")+
    geom_vline(xintercept=quantile(rsquared,c(0.025)), color="blue")+
    geom_vline(xintercept=quantile(rsquared,c(0.975)), color="blue")+
    labs(x="Pseudo R² between Ne/N and genetic diversity", 
         y = "Frequency")+  
    theme_classic()+
    labs(tag="B")
  print(p1)
  
  dd=data.frame(x=log(pvalue))
  p2<-ggplot(dd,aes(x=x)) +
    geom_density(fill="grey",alpha=0.15)+
    geom_vline(xintercept=log(0.000966), color="red",
               linetype="dashed")+
    geom_vline(xintercept=quantile(dd$x,c(0.025)), color="blue")+
    geom_vline(xintercept=quantile(dd$x,c(0.975)), color="blue")+
    labs(x="pvalue between adult lifespan and genetic diversity", 
         y = "Frequency")+  
    theme_classic()+
    labs(tag="B")
  print(p2)
  
  figure<-ggarrange(p,p1,ncol=1)
  pdf("C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT/figures/sensibility_agene_slope.pdf",width=7.5,height=5)
  print(figure)
  dev.off()
  
}

# Plot full model -----
data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=agene_output[[4]]$Output16,
                     x2=lfh$Parental_Care)
summary(lm(y~x,data=data_plot[data_plot$x2=="No",]))

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

div_agene<-ggplot(data_plot_whole, 
          aes(x = x, 
              y = mean*100,
              label=species,
              #group=x2)
          )) +
  geom_line(aes(colour=gg),
            size=1)+
  #geom_ribbon(aes(ymin=sd_025,
  #                ymax=sd_975,colour=gg),
  #            alpha=0.05,
  #            size=0.05) +
  geom_point(data=data_plot,size=2.5,
             aes(x=x,
                 y=y,
                 shape=x2))+
  geom_rangeframe()+
  theme_classic()+
  geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  scale_x_continuous(breaks=seq(0,0.7,by=0.1),
                     labels=as.character(seq(0,0.7,by=0.1)))+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                     labels=as.character(seq(0,1.5,by=0.25)))+
  xlab(expression(N[e]/N))+
  ylab("")+
  scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
  scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
  theme(legend.position = "bottom")+
  annotate("text", x = 0.2025, y = 1.325, label = "paste(italic(p), \" = 0.000966 \")",parse=T)+
  annotate("text", x = 0.2025, y = 1.265, label = "paste(italic(R) ^ 2, \" = 0.55 \")",parse=T)

x=7
pdf(paste(wd,"/figures/div_agene.pdf",sep=""),width=1*x,height=(2/3)*x)
print(div_agene)
dev.off()

# Pairwise-plot full model ----
sim=16
compare_het=data.frame(SP1=c(NA),
                       SP2=c(NA),
                       ratio_mean=c(NA),
                       ratio_truehet=c(NA),
                       cross=c(NA))

Sp=c("Lbude","Hgutt","Dlabr","Scant","Dpunt","Lmorm","Cgale","Scine","Mmerl","Ssard","Styph","Peryt",
     "Msurm","Cjuli","Scabr","Spilc")


for (j in 1:(nrow(agene_output[[1]])-1)){
  for (k in (j+1):nrow(agene_output[[1]])){
    vect=c(as.character(Sp[j]),
           as.character(Sp[k]),
           round(agene_output[[4]][which(agene_output[[4]]$Species_code==Sp[j]),sim+3]/agene_output[[4]][which(agene_output[[4]]$Species_code==Sp[k]),sim+3],3),
           round(div[which(names(div)==Sp[j])]/div[which(names(div)==Sp[k])],3)
    )
    compare_het=rbind(compare_het,vect)
  }
}
compare_het=compare_het[-1,]
compare_het$ratio_mean=as.numeric(compare_het$ratio_mean)
compare_het$ratio_truehet=as.numeric(compare_het$ratio_truehet)

# Parental care
compare_het_stat=compare_het
fit<-lm(compare_het_stat$ratio_mean~compare_het_stat$ratio_truehet)

pairwise_agene<-ggplot(compare_het_stat,aes(ratio_truehet,ratio_mean))+
  theme_classic()+
  geom_rangeframe()+
  geom_point(size=2)+
  geom_smooth(method="lm",fullrange="T")+
  geom_abline(intercept=0,slope=1,col="red")+
  xlab("")+
  ylab("")+
  ylim(c(0,max(compare_het_stat$ratio_mean)))+
  xlim(c(0,1))+
  ggtitle("Whole dataset")+
  annotate("text", x = 0.225, y = 4, label = "paste(italic(p), \" = 0.0205 \")",parse=T)+
  annotate("text", x = 0.225, y = 3.75, label = "Est. slope = 0.73")

# No parental care
compare_het_stat_nopc=compare_het
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Hgutt",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Hgutt",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Scine",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Scine",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Styph",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Styph",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Cgale",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Cgale",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP1!="Scant",]
compare_het_stat_nopc=compare_het_stat_nopc[compare_het_stat_nopc$SP2!="Scant",]

fit<-lm(compare_het_stat_nopc$ratio_mean~compare_het_stat_nopc$ratio_truehet)

pairwise_agene_nopc<-ggplot(compare_het_stat_nopc,aes(ratio_truehet,ratio_mean))+
  theme_classic()+
  geom_rangeframe()+
  geom_point(size=2)+
  geom_smooth(data=compare_het_stat_nopc,method="lm",fullrange="T")+
  geom_abline(intercept=0,slope=1,col="red")+
  xlab("")+
  ylab("")+
  ylim(c(0,max(compare_het_stat_nopc$ratio_mean)))+
  xlim(c(0,1))+
  ggtitle("Only non brooding species")+
  annotate("text", x = 0.23, y = 2.1, label = "paste(italic(p), \" = 0.000117 \")",parse=T)+
  annotate("text", x = 0.23, y = 1.969, label = "Est. slope = 0.94")

figure<-ggarrange(pairwise_agene,
                  pairwise_agene_nopc,ncol=2)

x=7
pdf(paste(wd,"/figures/div_agene_pairwise.pdf",sep=""),width=1*x,height=(2/3)*x)
print(figure)
dev.off()
