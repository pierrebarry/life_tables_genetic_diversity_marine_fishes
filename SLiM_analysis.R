#------------------------------------------------------------------------#
#                                                                        #
#                           SLiM analysis                                #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages ----
library(ggthemes)
library(ggpubr)
library(ggiraph)
library(ggrepel)
library(ggplot2)
library(viridis)
library(betareg)
library(readxl)
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
lfh<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet='lfh'))
load(file="Data/div.Rdata")
lfh$div=as.vector(div)
lfh$log_fec=log(lfh$Fec)
setwd("Data/forward_slim_old")
species=c("Cgale",
          "Cjuli",
          "Dlabr",
          "Dpunt",
          "Hgutt",
          "Lbude",
          "Lmorm",
          "Mmerl",
          "Msurm",
          "Peryt",
          "Scabr",
          "Scant",
          "Scine",
          "Spilc",
          "Ssard",
          "Styph")
# Check pop and genet file ----
output=paste("Output",seq(1,16),sep="")
K=2000
pop_species=data.frame(Species=species)
for (i in 1:16){
  pop_species=cbind(pop_species,rep(0,16))
}
colnames(pop_species)=c("Species",output)

est_species=data.frame(Species=species)
for (i in 1:16){
  est_species=cbind(est_species,rep(0,16))
}
colnames(est_species)=c("Species",output)
plot_true=1
    
for (sim in seq(1,8)){

    # Check pop file
    for (sp in seq(1,16)){
      
      print(paste(species[sp],"-",sim))
      pop_fit=c()
      pop_mean=c()
      test_pop=vector('list',1)
      
      for (i in seq(1,50)){
        if (file.exists(paste("Output",sim,"/poptotal_",species[sp],"_",K,"_",i-1,".txt",sep=""))){
          test_pop[[i]]<-read.table(paste("Output",sim,"/poptotal_",species[sp],"_",K,"_",i-1,".txt",sep=""),sep="\n")
        }
      }
      
      if (plot_true==1){
        
        plot(seq(100,24000,by=100),seq(0.0001,0.0075,length.out=length(seq(100,24000,by=100))),
             type="n",
             xlim=c(0,24000),
             ylim=c(0,max(unlist(test_pop),na.rm=T)),
             col="grey",
             xlab="Years",
             ylab="Neutral diversity",
             main=paste(species[sp]," - Sim:",sim,sep=""))
        
        
        x=seq(100,24000,by=100)
        
        for (i in seq(1,length(test_pop))){
          
          x<-seq(100,24000,by=100)
          y<-test_pop[[i]][1:240,1]
          lines(x,y,type="l",xlim=c(0,24000),col="grey",ylim=c(0,0.03),lwd=3)
          pop_fit=NA
          pop_mean=c(pop_mean,mean(y[1:240],na.rm=T))
          
        }
        
        mean_neutral=c()
        for (i in 1:nrow(test_pop[[1]])){
          mean_neutral=c(mean_neutral,mean(
            unlist(
              test_pop
            )[seq(i,length(test_pop)*nrow(test_pop[[1]])+1,by=nrow(test_pop[[1]]))],na.rm=T))
        }
        lines(x,as.vector(mean_neutral)[1:240],col="red",lwd=3)
        median_neutral=c()
        for (i in 1:nrow(test_pop[[1]])){
          median_neutral=c(median_neutral,median(
            unlist(
              test_pop
            )[seq(i,length(test_pop)*nrow(test_pop[[1]])+1,by=nrow(test_pop[[1]]))],na.rm=T))
        }
        lines(x,as.vector(median_neutral)[1:240],col="orange",lwd=3)
        
        text(17500,max(unlist(test_pop),na.rm=T)*0.8,
             paste("Mean :",round(mean(pop_mean),4),sep=""))
        
      } else {
        
        for (i in seq(1,length(test_pop))){
          y<-test_pop[[i]][1:240,1]
          
          pop_mean=c(pop_mean,mean(y[1:240],na.rm=T))
        }   
        
      }

      pop_species[sp,sim+1]=mean(pop_mean)
      
    }
    
    # Check genet file
    for (sp in seq(1,16)){
      
      genet_fit=c()
      genet_mean=c()
      test_genet=vector('list',1)
      for (i in seq(1,50)){
        if (file.exists(paste("Output",sim,"/div_",species[sp],"_",K,"_",i-1,".txt",sep=""))){
          test_genet[[i]]<-read.table(paste("Output",sim,"/div_",species[sp],"_",K,"_",i-1,".txt",sep=""),sep="\n")
        }
      }
      
      if (plot_true==1){
        
        plot(seq(100,24000,by=100),seq(0.0001,0.0075,length.out=length(seq(100,24000,by=100))),
             type="n",
             xlim=c(0,24000),
             ylim=c(0,max(unlist(test_genet),na.rm=T)),
             col="grey",
             xlab="Years",
             ylab="Neutral diversity",
             main=paste(species[sp]," - Sim:",sim,sep=""))
        
        x=seq(100,24000,by=100)
        
        for (i in seq(1,length(test_genet))){
          
          x<-seq(100,24000,by=100)
          y<-test_genet[[i]][1:240,1]
          lines(x,y,type="l",xlim=c(0,24000),col="grey",ylim=c(0,0.03),lwd=3)
          genet_fit=NA
          genet_mean=c(genet_mean,mean(y[150:200],na.rm=T))
          
        }
        
        mean_neutral=c()
        for (i in 1:nrow(test_genet[[1]])){
          mean_neutral=c(mean_neutral,mean(
            unlist(
              test_genet
            )[seq(i,length(test_genet)*nrow(test_genet[[1]])+1,by=nrow(test_genet[[1]]))],na.rm=T))
        }
        lines(x,as.vector(mean_neutral)[1:240],col="red",lwd=3)
        median_neutral=c()
        for (i in 1:nrow(test_genet[[1]])){
          median_neutral=c(median_neutral,median(
            unlist(
              test_genet
            )[seq(i,length(test_genet)*nrow(test_genet[[1]])+1,by=nrow(test_genet[[1]]))],na.rm=T))
        }
        lines(x,as.vector(median_neutral)[1:240],col="orange",lwd=3)
        
        text(5000,max(unlist(test_genet),na.rm=T)*0.9,
             paste("Mean :",round(mean(genet_mean),4),sep=""))
        
        
      } else {
        
        for (i in seq(1,length(test_genet))){
          
          x<-seq(100,24000,by=100)
          y<-test_genet[[i]][1:240,1]
          genet_mean=c(genet_mean,median(y[150:240],na.rm=T))
          
        }
        
      }

      est_species[sp,sim+1]=median(genet_mean)
      
    }

}

# Save new output (if wanted) ----
save_new_output=0
if (save_new_output==1){
  save(est_species,file="est_species.Rdata")
  save(pop_species,file="pop_species.Rdata")
} else {
  load(file="est_species.Rdata")
  load(file="pop_species.Rdata")
}

# Analysis of SLiM output -----
est_species=est_species[order(est_species$Species),]
pop_species=pop_species[order(pop_species$Species),]
forward_test=data.frame(Sim=seq(1,16),Surv=rep(c(rep(0,1),rep(1,1)),8) ,
                        Fec=rep(c(rep(0,2),rep(1,2)),4),
                        Sex=rep(c(rep(0,4),rep(1,4)),2),
                        AgeMat=rep(c(rep(0,8),rep(1,8))),
                        test_slope=rep(0,16),test_intercept=rep(0,16),
                        slope_mean=rep(0,16),slope_sd=rep(0,16),
                        inter_mean=rep(0,16),inter_sd=rep(0,16),
                        pvalue=rep(0,16),
                        test_slope_nopc=rep(0,16),test_intercept_nopc=rep(0,16),
                        slope_mean_nopc=rep(0,16),slope_sd_nopc=rep(0,16),
                        inter_mean_nopc=rep(0,16),inter_sd_nopc=rep(0,16),
                        pvalue_nopc=rep(0,16)
                        )

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

est_species$onclick <- sprintf("window.open(\"%s%s\")",
                                "http://en.wikipedia.org/wiki/", as.character(wikipedia) )

p_slim=vector('list',16)
p_slim_pairwise=vector('list',16)
p_slim_pairwise_nopc=vector('list',16)

for (sim in 1:8){
  
  title=c("--",
          "Surv",
          "Fec",
          "Surv + Fec",
          "Sex",
          "Surv + Sex",
          "Fec + Sex",
          "Surv + Fec + Sex",
          "AgeMat",
          "AgeMat + Surv",
          "AgeMat + Fec",
          "AgeMat + Surv + Fec",
          "AgeMat + Sex",
          "AgeMat + Surv + Sex",
          "AgeMat + Fec + Sex",
          "AgeMat + Surv + Fec + Sex")
  
  p_slim[[sim]]<-local({
    data_plot=data.frame(species=lfh$Species_plot,
                         y=lfh$div,
                         x=est_species[,sim+1],
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
  
  for (j in 1:(nrow(est_species)-1)){
    for (k in (j+1):nrow(est_species)){
      vect=c(as.character(Sp[j]),
             as.character(Sp[k]),
             round(est_species[which(est_species$Species==Sp[j]),sim+1]/est_species[which(est_species$Species==Sp[k]),sim+1],3),
             round(div[which(names(div)==Sp[j])]/div[which(names(div)==Sp[k])],3))
      compare_het=rbind(compare_het,vect)
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
    #geom_point_interactive(aes(x = ratio_truehet, y = ratio_mean,
    #                           tooltip =SP1, data_id = SP1))+
    xlab("Ratio of estimated NeN")+
    ylab("Ratio of observed genetic diversity")+
    ylim(c(0,max(compare_het_stat$ratio_mean)))+
    xlim(c(0,1))+
    ggtitle("Whole dataset")
  
  p_slim_pairwise[[sim]] <- local({
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
  
  
  forward_test$test_slope[sim]=test_slope_different_1(fit=fit)
  forward_test$test_intercept[sim]=test_intercept_different_0(fit=fit)
  forward_test$slope_mean[sim]=coefficients(summary(fit))[2,1]
  forward_test$slope_sd[sim]=coefficients(summary(fit))[2,2]
  forward_test$inter_mean[sim]=coefficients(summary(fit))[1,1]
  forward_test$inter_sd[sim]=coefficients(summary(fit))[1,2] 
  forward_test$pvalue[sim]=coefficients(summary(fit))[2,4] 
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
  
  p_slim_pairwise_nopc[[sim]] <- local({
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
  
  forward_test$test_slope_nopc[sim]=test_slope_different_1(fit=fit)
  forward_test$test_intercept_nopc[sim]=test_intercept_different_0(fit=fit)
  forward_test$slope_mean_nopc[sim]=coefficients(summary(fit))[2,1]
  forward_test$slope_sd_nopc[sim]=coefficients(summary(fit))[2,2]
  forward_test$inter_mean_nopc[sim]=coefficients(summary(fit))[1,1]
  forward_test$inter_sd_nopc[sim]=coefficients(summary(fit))[1,2] 
  forward_test$pvalue_nopc[sim]=coefficients(summary(fit))[2,4]
  
}


pdf(paste(wd,"/figures/slim_suppmat.pdf",sep=""),width=30,height=30)
print(ggarrange(p_slim[[1]],
                p_slim[[9]],
                p_slim[[2]],
                p_slim[[10]],
                p_slim[[3]],
                p_slim[[11]],
                p_slim[[4]],
                p_slim[[12]],
                p_slim[[5]],
                p_slim[[13]],
                p_slim[[6]],
                p_slim[[14]],
                p_slim[[7]],
                p_slim[[15]],
                p_slim[[8]],
                p_slim[[16]],
                nrow=4,
                ncol=4))
dev.off()

figure1<-annotate_figure(ggarrange(p_slim_pairwise[[1]], p_slim_pairwise_nopc[[1]], ncol = 2, nrow = 1),
                         top=title[1])
figure2<-annotate_figure(ggarrange(p_slim_pairwise[[2]], p_slim_pairwise_nopc[[2]], ncol = 2, nrow = 1),
                         top=title[2])
figure3<-annotate_figure(ggarrange(p_slim_pairwise[[3]], p_slim_pairwise_nopc[[3]], ncol = 2, nrow = 1),
                         top=title[3])
figure4<-annotate_figure(ggarrange(p_slim_pairwise[[4]], p_slim_pairwise_nopc[[4]], ncol = 2, nrow = 1),
                         top=title[4])
figure5<-annotate_figure(ggarrange(p_slim_pairwise[[5]], p_slim_pairwise_nopc[[5]], ncol = 2, nrow = 1),
                         top=title[5])
figure6<-annotate_figure(ggarrange(p_slim_pairwise[[6]], p_slim_pairwise_nopc[[6]], ncol = 2, nrow = 1),
                         top=title[6])
figure7<-annotate_figure(ggarrange(p_slim_pairwise[[7]], p_slim_pairwise_nopc[[7]], ncol = 2, nrow = 1),
                         top=title[7])
figure8<-annotate_figure(ggarrange(p_slim_pairwise[[8]], p_slim_pairwise_nopc[[8]], ncol = 2, nrow = 1),
                         top=title[8])
figure9<-annotate_figure(ggarrange(p_slim_pairwise[[9]], p_slim_pairwise_nopc[[9]], ncol = 2, nrow = 1),
                         top=title[9])
figure10<-annotate_figure(ggarrange(p_slim_pairwise[[10]], p_slim_pairwise_nopc[[10]], ncol = 2, nrow = 1),
                          top=title[10])
figure11<-annotate_figure(ggarrange(p_slim_pairwise[[11]], p_slim_pairwise_nopc[[11]], ncol = 2, nrow = 1),
                          top=title[11])
figure12<-annotate_figure(ggarrange(p_slim_pairwise[[12]], p_slim_pairwise_nopc[[12]], ncol = 2, nrow = 1),
                          top=title[12])
figure13<-annotate_figure(ggarrange(p_slim_pairwise[[13]], p_slim_pairwise_nopc[[13]], ncol = 2, nrow = 1),
                          top=title[13])
figure14<-annotate_figure(ggarrange(p_slim_pairwise[[14]], p_slim_pairwise_nopc[[14]], ncol = 2, nrow = 1),
                          top=title[14])
figure15<-annotate_figure(ggarrange(p_slim_pairwise[[15]], p_slim_pairwise_nopc[[15]], ncol = 2, nrow = 1),
                          top=title[15])
figure16<-annotate_figure(ggarrange(p_slim_pairwise[[16]], p_slim_pairwise_nopc[[16]], ncol = 2, nrow = 1),
                          top=title[16])
pdf(paste(wd,"/figures/slim_pairwise_suppmat.pdf",sep=""),width=30,height=30)
print(ggarrange(figure1,
                figure9,
                figure2,
                figure10,
                figure3,
                figure11,
                figure4,
                figure12,
                figure5,
                figure13,
                figure6,
                figure14,
                figure7,
                figure15,
                figure8,
                figure16,
                nrow=4,
                ncol=4))
dev.off()


# All possible contaminations----
combi=combn(seq(1,16),11)
slope=c()
pvalue=c()
rsquared=c()

if (time_code==1){
  for (i in 1:ncol(combi)){
    aaa<-betareg(I(lfh[combi[,i],]$div/100)~est_species[combi[,i],]$Output7)
    slope=c(slope,coefficients(summary(aaa))$mean[2,1])
    pvalue=c(pvalue,coefficients(summary(aaa))$mean[2,4])
    rsquared=c(rsquared,aaa$pseudo.r.squared)
  }
  
  dd=data.frame(x=slope)
  p<-ggplot(dd,aes(x=x)) +
    geom_density(fill="grey",alpha=0.15)+
    geom_vline(xintercept=as.numeric(forward_test$slope_mean_nopc[7]), color="red",
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
    geom_vline(xintercept=0.5158, color="red",
               linetype="dashed")+
    geom_vline(xintercept=quantile(rsquared,c(0.025)), color="blue")+
    geom_vline(xintercept=quantile(rsquared,c(0.975)), color="blue")+
    labs(x="Pseudo R² between adult lifespan and genetic diversity", 
         y = "Frequency")+  
    theme_classic()+
    labs(tag="A")
  print(p1)
  
  dd=data.frame(x=pvalue)
  p2<-ggplot(dd,aes(x=x)) +
    geom_density(fill="grey",alpha=0.15)+
    geom_vline(xintercept=log(0.00597), color="red",
               linetype="dashed")+
    geom_vline(xintercept=quantile(dd$x,c(0.025)), color="blue")+
    geom_vline(xintercept=quantile(dd$x,c(0.975)), color="blue")+
    labs(x="pvalue between adult lifespan and genetic diversity", 
         y = "Frequency")+  
    theme_classic()+
    labs(tag="B")
  print(p2)
  
  figure<-ggarrange(p1,p2,ncol=1)
  pdf("../../figures/sensibility_slim_slope.pdf",width=7.5,height=5)
  print(figure)
  dev.off()
  
  forward_test[c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16),c(1,5,2,3,4,6,7,12,13)]
  
}

# Full model plot ----
data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=est_species$Output8,
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

div_slim<-ggplot(data_plot_whole, 
                 aes(x = x, 
                     y = mean*100,
                     label=species
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
  scale_x_continuous(breaks=seq(0,0.08,by=0.01),
                     labels=as.character(seq(0,0.08,by=0.01))
                     )+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25),
                     labels=as.character(seq(0,1.5,by=0.25)))+
  xlab("Simulated heterozygosity (%)")+
  ylab("")+
  scale_color_manual(name = "Model",values=viridis(100)[c(47.5,72.5,97.5)],labels=c("Whole dataset","Only non brooding species","No parental care species"))+
  scale_shape_manual(name = "Brooding behaviour",values=c(19, 1))+
  theme(legend.position = "bottom")+
  annotate("text", x = 0.015, y = 1.325, label = "paste(italic(p), \" = 0.0115 \")",parse=T)+
  annotate("text", x = 0.015, y = 1.265, label = "paste(italic(R) ^ 2, \" = 0.4347 \")",parse=T)


x=7
pdf(paste(wd,"/figures/div_slim.pdf",sep=""),width=1*x,height=(2/3)*x)
print(div_slim)
dev.off()

# Full model Pairwise-plot ----
sim=8
compare_het=data.frame(SP1=c(NA),
                       SP2=c(NA),
                       ratio_mean=c(NA),
                       ratio_truehet=c(NA))

Sp=c("Lbude","Hgutt","Dlabr","Scant","Dpunt","Lmorm","Cgale","Scine","Mmerl","Ssard","Styph","Peryt",
     "Msurm","Cjuli","Scabr","Spilc")

for (j in 1:(nrow(est_species)-1)){
  for (k in (j+1):nrow(est_species)){
    vect=c(as.character(Sp[j]),
           as.character(Sp[k]),
           round(est_species[which(est_species$Species==Sp[j]),sim+1]/est_species[which(est_species$Species==Sp[k]),sim+1],3),
           round(div[which(names(div)==Sp[j])]/div[which(names(div)==Sp[k])],3))
    compare_het=rbind(compare_het,vect)
  }
}
compare_het=compare_het[-1,]
compare_het$ratio_mean=as.numeric(compare_het$ratio_mean)
compare_het$ratio_truehet=as.numeric(compare_het$ratio_truehet)
# Parental care
compare_het_stat=compare_het
fit<-lm(compare_het_stat$ratio_mean~compare_het_stat$ratio_truehet)

pairwise_slim<-ggplot(compare_het_stat[compare_het_stat$ratio_mean<10,],aes(ratio_truehet,ratio_mean))+
  theme_classic()+
  geom_rangeframe()+
  geom_point(size=2)+
  geom_smooth(method="lm",fullrange="T")+
  geom_abline(intercept=0,slope=1,col="red")+
  xlab("")+
  ylab("")+
  ylim(c(0,max(compare_het_stat[compare_het_stat$ratio_mean<10,]$ratio_mean)))+
  xlim(c(0,1))+
  ggtitle("Whole dataset")+
  annotate("text", x = 0.225, y = 5.85, label = "paste(italic(p), \" = 0.0476 \")",parse=T)+
  annotate("text", x = 0.25, y = 5.4, label = "Est. slope = 0.7329")

# No parental care
compare_het_stat=compare_het
compare_het_stat=compare_het_stat[compare_het_stat$SP1!="Hgutt",]
compare_het_stat=compare_het_stat[compare_het_stat$SP2!="Hgutt",]
compare_het_stat=compare_het_stat[compare_het_stat$SP1!="Scine",]
compare_het_stat=compare_het_stat[compare_het_stat$SP2!="Scine",]
compare_het_stat=compare_het_stat[compare_het_stat$SP1!="Styph",]
compare_het_stat=compare_het_stat[compare_het_stat$SP2!="Styph",]
compare_het_stat=compare_het_stat[compare_het_stat$SP1!="Cgale",]
compare_het_stat=compare_het_stat[compare_het_stat$SP2!="Cgale",]
compare_het_stat=compare_het_stat[compare_het_stat$SP1!="Scant",]
compare_het_stat=compare_het_stat[compare_het_stat$SP2!="Scant",]

fit<-lm(compare_het_stat$ratio_mean~compare_het_stat$ratio_truehet)

pairwise_slim_nopc<-ggplot(compare_het_stat,aes(ratio_truehet,ratio_mean))+
  theme_classic()+
  geom_rangeframe()+
  geom_point(size=2)+
  geom_smooth(method="lm",fullrange="T")+
  geom_abline(intercept=0,slope=1,col="red")+
  xlab("")+
  ylab("")+
  ylim(c(0,max(compare_het_stat$ratio_mean)))+
  xlim(c(0,1))+
  ggtitle("Only non brooding species")+
  annotate("text", x = 0.225, y = 2.7, label = "paste(italic(p), \" = 0.00357 \")",parse=T)+
  annotate("text", x = 0.25, y = 2.5, label = "Est. slope = 0.9493")


figure<-ggarrange(p1, p2, ncol = 2, nrow = 1)

x=7
pdf(paste(wd,"/figures/div_slim_pairwise.pdf",sep=""),width=1*x,height=(2/3)*x)
print(div_slim)
dev.off()

# Supp Mat figures ----
# Check pop file
par(mfrow=c(4,4),
    mar=c(2,2,2,2),
    mai=c(0.5,0.5,0.25,0.05))
plot_true=1
for (sim in seq(1,8)){
  for (sp in seq(1,16)){
    
    print(paste(species[sp],"-",sim))
    pop_fit=c()
    pop_mean=c()
    test_pop=vector('list',1)
    
    for (i in seq(1,50)){
      if (file.exists(paste("Output",sim,"/poptotal_",species[sp],"_",K,"_",i-1,".txt",sep=""))){
        test_pop[[i]]<-read.table(paste("Output",sim,"/poptotal_",species[sp],"_",K,"_",i-1,".txt",sep=""),sep="\n")
      }
    }
    
    if (plot_true==1){
      
      plot(seq(100,24000,by=100),seq(0.0001,0.0075,length.out=length(seq(100,24000,by=100))),
           type="n",
           xlim=c(0,24000),
           ylim=c(0,max(unlist(test_pop),na.rm=T)),
           col="grey",
           xlab="",
           yaxt="n",
           xaxt="n",
           ylab="",
           main=species[sp])
      axis(2, at=c(0,1000,2000),labels=c(0,1000,2000))
      axis(1, at=seq(0,25000,by=5000),labels=seq(0,25000,by=5000))
      
      
      x=seq(100,24000,by=100)
      
      for (i in seq(1,length(test_pop))){
        
        x<-seq(100,24000,by=100)
        y<-test_pop[[i]][1:240,1]
        lines(x,y,type="l",xlim=c(0,24000),col="grey",ylim=c(0,0.03),lwd=3)
        pop_fit=NA
        pop_mean=c(pop_mean,mean(y[1:240],na.rm=T))
        
      }
      
      mean_neutral=c()
      for (i in 1:nrow(test_pop[[1]])){
        mean_neutral=c(mean_neutral,mean(
          unlist(
            test_pop
          )[seq(i,length(test_pop)*nrow(test_pop[[1]])+1,by=nrow(test_pop[[1]]))],na.rm=T))
      }
      lines(x,as.vector(mean_neutral)[1:240],col="red",lwd=3)
      median_neutral=c()
      for (i in 1:nrow(test_pop[[1]])){
        median_neutral=c(median_neutral,median(
          unlist(
            test_pop
          )[seq(i,length(test_pop)*nrow(test_pop[[1]])+1,by=nrow(test_pop[[1]]))],na.rm=T))
      }
      lines(x,as.vector(median_neutral)[1:240],col="orange",lwd=3)
      
      text(17500,max(unlist(test_pop),na.rm=T)*0.7,
           paste("Mean :",round(mean(pop_mean),4),sep=""))
      
    } else {
      
      for (i in seq(1,length(test_pop))){
        y<-test_pop[[i]][1:240,1]
        
        pop_mean=c(pop_mean,mean(y[1:240],na.rm=T))
      }   
      
    }
    
  }
  
  mtext("Population size",side=2,line=-1.5,outer=TRUE,cex=1,las=0)
  mtext("Years",side=1,line=-1,outer=TRUE,cex=1,las=0)
  dev.copy2pdf(file=paste(wd,"/figures/pop_check_sim",sim,".pdf",sep=""), width = 15, height = 10)
  
}

# Check genet file
par(mfrow=c(4,4),
    mar=c(2,1,2,2),
    mai=c(0.5,0.5,0.25,0.05))
plot_true=1
for (sim in seq(1,8)){  
  for (sp in seq(1,16)){
    
    genet_fit=c()
    genet_mean=c()
    test_genet=vector('list',1)
    for (i in seq(1,50)){
      if (file.exists(paste("Output",sim,"/div_",species[sp],"_",K,"_",i-1,".txt",sep=""))){
        test_genet[[i]]<-read.table(paste("Output",sim,"/div_",species[sp],"_",K,"_",i-1,".txt",sep=""),sep="\n")
      }
    }
    
    if (plot_true==1){

      plot(seq(100,24000,by=100),seq(0.0001,0.0075,length.out=length(seq(100,24000,by=100))),
           type="n",
           xlim=c(0,24000),
           ylim=c(0,max(unlist(test_genet),na.rm=T)),
           col="grey",
           xlab="",
           xaxt="n",
           ylab="",
           main=species[sp])
      axis(1, at=seq(0,25000,by=5000),labels=seq(0,25000,by=5000))
      
      
      x=seq(100,24000,by=100)
      
      for (i in seq(1,length(test_genet))){
        
        x<-seq(100,24000,by=100)
        y<-test_genet[[i]][1:240,1]
        lines(x,y,type="l",xlim=c(0,24000),col="grey",ylim=c(0,0.03),lwd=3)
        genet_fit=NA
        genet_mean=c(genet_mean,mean(y[150:200],na.rm=T))
        
      }
      
      mean_neutral=c()
      for (i in 1:nrow(test_genet[[1]])){
        mean_neutral=c(mean_neutral,mean(
          unlist(
            test_genet
          )[seq(i,length(test_genet)*nrow(test_genet[[1]])+1,by=nrow(test_genet[[1]]))],na.rm=T))
      }
      lines(x,as.vector(mean_neutral)[1:240],col="red",lwd=3)
      median_neutral=c()
      for (i in 1:nrow(test_genet[[1]])){
        median_neutral=c(median_neutral,median(
          unlist(
            test_genet
          )[seq(i,length(test_genet)*nrow(test_genet[[1]])+1,by=nrow(test_genet[[1]]))],na.rm=T))
      }
      lines(x,as.vector(median_neutral)[1:240],col="orange",lwd=3)
      
      text(5000,max(unlist(test_genet),na.rm=T)*0.9,
           paste("Mean :",round(mean(genet_mean),4),sep=""))
      abline(h=0.08,col="chartreuse4",lty=2)
      
      
    } else {
      
      for (i in seq(1,length(test_genet))){
        
        x<-seq(100,24000,by=100)
        y<-test_genet[[i]][1:240,1]
        genet_mean=c(genet_mean,median(y[150:240],na.rm=T))
        
      }
      
    }
    
  }
  
  mtext("Genetic diversity (%)",side=2,line=-1.5,outer=TRUE,cex=1,las=0)
  mtext("Years",side=1,line=-1,outer=TRUE,cex=1,las=0)
  dev.copy2pdf(file=paste(wd,"/figures/simu_gen_sim",sim,".pdf",sep=""), width = 15, height = 10)
  
  
}
