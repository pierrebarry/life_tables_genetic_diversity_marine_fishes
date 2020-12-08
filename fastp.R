#------------------------------------------------------------------------#
#                                                                        #
#                               FASTP                                    #
#                                                                        #
#------------------------------------------------------------------------#

# Load packages -----
library(ggplot2)
library("RColorBrewer")
library(ggpubr)
library("gridExtra")
library(ggiraph)
library(ggExtra)
# Load data ----
load(file="Data/fastp_process.Rdata")
# Species plot ----
analogy_legend=c(NA,NA,NA,
                 "Total reads before filtering ",    
                 "Total bases before filtering (Gb)",
                 "GC content before filtering (%)",
                 "Q20 rate before filtering (%)",
                 "Q30 rate before filtering (%)",
                 "Total reads after filtering (10⁹)",
                 "Total bases after filtering (Gb)",
                 "GC content after filtering (%)",
                 "Q20 rate after filtering (%)",
                 "Q30 rate after filtering (%)",
                 "Passed filter reads (%)",
                 "Corrected reads (%)",
                 "Corrected bases (%)",
                 "Low quality reads (%)",
                 "Too many N reads (%)",
                 "Too short reads (%)",
                 "Low complexity reads (%)",
                 "Duplication rate (%)")

fastp_plot=data.frame(LOCATION=rep(fastp_process$LOCATION,18),
                      SPECIES=rep(fastp_process$SPECIES,18),
                      SAMPLE=rep(fastp_process$SAMPLE,18),
                      param=rep(NA,nrow(fastp_process)*18),
                      value=rep(0,nrow(fastp_process)*18),
                      color=rep(fastp_process$color,18))

a=1
for (i in 4:21){
  for (j in 1:nrow(fastp_process)){
    fastp_plot$param[a]=colnames(fastp_process)[i]
    fastp_plot$value[a]=fastp_process[j,i]
    a=a+1
  }
} 

fastp_plot=fastp_plot[is.na(fastp_plot$param)==FALSE,]


labels_species=c()

for (sp in 1:length(levels(fastp_process$SPECIES))){
  
  labels_species=c(labels_species,
                   paste(levels(fastp_process$SPECIES)[sp]," \n n = ",table(fastp_process$SPECIES)[sp],sep=" "))
  
}

labels_species=as.vector(labels_species)

fastp_process$LOCATION=factor(fastp_process$LOCATION,levels=c("Li","Mu","Fa","Ga"))

# Plot Supplementary Figure 1 ----
p1 <- ggplot(fastp_process, aes(x=SPECIES, y=TOTAL_READS_AFTER_FILTERING ),color=LOCATION) + 
  geom_point(aes(x=SPECIES, y=TOTAL_READS_AFTER_FILTERING,color=color),size = 1.5,position="jitter",
             alpha=0.75,show.legend = F)+
  geom_boxplot(outlier.shape = NA,alpha=0,show.legend = F)+
  geom_violin(color="black",alpha=0,show.legend=F) +
  geom_jitter(position=position_jitter(0.2),color=fastp_process$color,size=0.00001,alpha=0,show.legend = F)+
  scale_x_discrete(name="Species",
                   labels=labels_species)+
  theme_bw()+
  ylab(expression("Nb. of reads (x10"^9*"bp)"))+
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Murcia","Faro","Gulf of Gascogne"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p1<-ggMarginal(p1, type="density",margins='y',alpha=0.75,colour = 'black',fill = 'chartreuse4',cex=0.5)

p2 <- ggplot(fastp_process, aes(x=SPECIES, y=Q30_RATE_AFTER_FILTERING),color=fastp_process$color) + 
  geom_point(aes(x=SPECIES, y=Q30_RATE_AFTER_FILTERING,color=color),size = 1.5,position="jitter",
             alpha=0.75,show.legend = F)+
  geom_boxplot(outlier.shape = NA,alpha=0,show.legend = F)+
  geom_violin(color="black",alpha=0,show.legend=F) +
  geom_jitter(position=position_jitter(0.2),color=fastp_process$color,size=0.00001,alpha=0,show.legend = F)+
  scale_x_discrete(name="Species",
                   labels=labels_species)+
  theme_bw()+
  ylab(expression("Q30 (%)"))+
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Murcia","Faro","Gulf of Gascogne"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2<-ggMarginal(p2, type="density",margins='y',alpha=0.75,colour = 'black',fill = 'chartreuse4',cex=0.5)

p3 <- ggplot(fastp_process, aes(x=SPECIES, y=DUPLICATION_RATE),color=fastp_process$color) + 
  geom_point(aes(x=SPECIES, y=DUPLICATION_RATE,color=color),size = 1.5,position="jitter",
             alpha=0.75,show.legend = F)+
  geom_boxplot(outlier.shape = NA,alpha=0,show.legend = F)+
  geom_violin(color="black",alpha=0,show.legend=F) +
  geom_jitter(position=position_jitter(0.2),color=fastp_process$color,size=0.00001,alpha=0,show.legend = F)+
  scale_x_discrete(name="Species",
                   labels=labels_species)+
  theme_bw()+
  ylab(expression("Duplication rate (%)"))+
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Murcia","Faro","Gulf of Gascogne"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3<-ggMarginal(p3, type="density",margins='y',alpha=0.75,colour = 'black',fill = 'chartreuse4',cex=0.5)

p4 <- ggplot(fastp_process, aes(x=SPECIES, y=GC_CONTENT_AFTER_FILTERING),color=LOCATION) + 
  geom_point(aes(x=SPECIES, y=GC_CONTENT_AFTER_FILTERING,color=LOCATION),size = 1.5,position="jitter",
             alpha=0.75,show.legend = T)+
  geom_boxplot(outlier.shape = NA,alpha=0,show.legend = F)+
  geom_violin(color="black",alpha=0,show.legend=F) +
  geom_jitter(position=position_jitter(0.2),color=fastp_process$color,size=0.00001,alpha=0,show.legend = F)+
  scale_x_discrete(name="Species",
                   labels=labels_species)+
  theme_bw()+
  ylab("GC content (%)")+
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(legend.position="bottom")
p4<-ggMarginal(p4, type="density",margins='y',alpha=0.75,colour = 'black',fill = 'chartreuse4',cex=0.5)

pdf(paste("figures/Fastp.pdf",sep=""),width=10,height=8.5)
grid.arrange(p1,
             p2,
             p3,
             p4,
             ncol = 1,
             heights = c(1,1,1,1.5))
dev.off()


# Fastp stats ----
#get quantity of reads after filtering
mean(fastp_process$TOTAL_READS_AFTER_FILTERING)
sd(fastp_process$TOTAL_READS_AFTER_FILTERING)
anova(lm(TOTAL_READS_AFTER_FILTERING~SPECIES,data=fastp_process))

anova(lm(TOTAL_READS_AFTER_FILTERING~LOCATION,data=fastp_process))
scale_centered_species=c()
for(i in 1:nrow(fastp_process)){
  scale_centered_species[i]=(fastp_process$TOTAL_READS_AFTER_FILTERING[i]-mean(fastp_process[fastp_process$SPECIES==fastp_process$SPECIES[i],]$TOTAL_READS_AFTER_FILTERING))/(sd(fastp_process[fastp_process$SPECIES==fastp_process$SPECIES[i],]$TOTAL_READS_AFTER_FILTERING))
}
anova(lm(scale_centered_species~fastp_process$LOCATION))

#get Q30 after filtering
mean(fastp_process$Q30_RATE_AFTER_FILTERING)
sd(fastp_process$Q30_RATE_AFTER_FILTERING)
anova(lm(Q30_RATE_AFTER_FILTERING~SPECIES,data=fastp_process))

anova(lm(Q30_RATE_AFTER_FILTERING~LOCATION,data=fastp_process))
scale_centered_species=c()
for(i in 1:nrow(fastp_process)){
  scale_centered_species[i]=(fastp_process$Q30_RATE_AFTER_FILTERING[i]-mean(fastp_process[fastp_process$SPECIES==fastp_process$SPECIES[i],]$Q30_RATE_AFTER_FILTERING))/(sd(fastp_process[fastp_process$SPECIES==fastp_process$SPECIES[i],]$Q30_RATE_AFTER_FILTERING))
}
summary(lm(scale_centered_species~fastp_process$LOCATION))
anova(lm(scale_centered_species~fastp_process$LOCATION))

pairwise.t.test(fastp_process$Q30_RATE_AFTER_FILTERING,fastp_process$SPECIES)
summary(lm(Q30_RATE_AFTER_FILTERING~SPECIES*LOCATION,data=fastp_process))

#get mean sd of duplication rate
mean(fastp_process$DUPLICATION_RATE)
sd(fastp_process$DUPLICATION_RATE)
anova(lm(DUPLICATION_RATE~SPECIES,data=fastp_process))
scale_centered_species=c()
for(i in 1:nrow(fastp_process)){
  scale_centered_species[i]=(fastp_process$DUPLICATION_RATE[i]-mean(fastp_process[fastp_process$LOCATION==fastp_process$LOCATION[i],]$DUPLICATION_RATE))/(sd(fastp_process[fastp_process$LOCATION==fastp_process$LOCATION[i],]$DUPLICATION_RATE))
}
anova(lm(scale_centered_species~fastp_process$SPECIES))

anova(lm(DUPLICATION_RATE~LOCATION,data=fastp_process))
scale_centered_species=c()
for(i in 1:nrow(fastp_process)){
  scale_centered_species[i]=(fastp_process$DUPLICATION_RATE[i]-mean(fastp_process[fastp_process$SPECIES==fastp_process$SPECIES[i],]$DUPLICATION_RATE))/(sd(fastp_process[fastp_process$SPECIES==fastp_process$SPECIES[i],]$DUPLICATION_RATE))
}
anova(lm(scale_centered_species~fastp_process$LOCATION))

pairwise.t.test(fastp_process$DUPLICATION_RATE,fastp_process$SPECIES)
summary(lm(DUPLICATION_RATE~SPECIES*LOCATION,data=fastp_process))

# Plot notebook ----

stats=c("Total reads",
        "Total bases",
        "GC content (%)",
        "Q20",
        "Q30",
        "Passed filter reads",
        "Corrected reads",
        "Corrected bases",
        "Low quality reads",
        "Too many N reads",
        "Too short reads",
        "Low complexity reads",
        "Duplication rate (%)")

p<-vector('list',length(stats))

for (j in seq(9,21)){
  
  data_tmp=fastp_process[,c(1,2,3,j)]
  colnames(data_tmp)[4]="val"
  
  p[[j-8]] <- local({
    ppp <- ggplot(data_tmp, aes(x=SPECIES, y=val,
                  text=paste("Species:",SPECIES,sep="")
                  )) + 
      geom_boxplot(outlier.shape = NA,
                   outlier.size = 0.000001,
                   alpha=0,
                   show.legend = F)+
      geom_violin(color="black",alpha=0,show.legend=F) +
      scale_x_discrete(name="Species",
                       labels=labels_species)+
      geom_point_interactive(
        aes(tooltip = paste("Species:",SPECIES,
                            "\n Sample:",SAMPLE,
                            "\n Location:",LOCATION,
                            "\n Val:",val),data_id = SAMPLE), 
        size = 1,
        col=fastp_process$color,
        position='jitter',
        alpha=0.75)+ 
      theme_classic()+
      ylab(stats[j-8])+
      scale_colour_manual(name = "",
                          labels = c("Gulf of Gascogne","Faro","Murcia","Gulf of Lion"),
                          values = brewer.pal(n = 4, name = "RdBu"))+
      theme(legend.position="bottom")
  })
  
  
  
  
}

names(p)=stats
p_fastp=p
