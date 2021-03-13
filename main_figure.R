#------------------------------------------------------------------------#
#                                                                        #
#                           FIGURES                                      #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages ----
library(heatmap3)
library(ggplotify)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(png)
library(grid)
library(ggpubr)
library(ggiraph)
library(igraph)
library(magick)
library(emdbook)
library(readxl)
library(betareg)
library(viridis)
library(ggthemes)
library(ggrepel)
# Load functions ----
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

# Figure 1 -----
load(file="Data/Summary_GenomeScope.Rdata")
load(file="Data/div.Rdata")
ss=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]
scale_centered_species=c()
for(i in 1:nrow(ss)){
  scale_centered_species[i]=(ss$Heterozygosity[i]-mean(ss[ss$Species==ss$Species[i],]$Heterozygosity))/(sd(ss[ss$Species==ss$Species[i],]$Heterozygosity))
}
div_loc=tapply(scale_centered_species,
               list(ss$Species,ss$Location),
               mean)

row.names(div_loc)=c("C. galerita",
                     "C. julis",
                     "D. labrax",
                     "D. puntazzo",
                     "H. guttulatus",
                     "L. budegassa",
                     "L. mormyrus",
                     "M. merluccius",
                     "M. surmuletus",
                     "P. erythrinus",
                     "S. cabrilla",
                     "S. cantharus",
                     "S. cinereus",
                     "S. pilchardus",
                     "S. sarda",
                     "S. typhle")

colnames(div_loc)=c("Algarve",
                    "Bay of \n Biscay ",
                    "Gulf of \n Lion   ",
                    "Costa \n Calida")
heat<-as.ggplot(pheat<-pheatmap(div_loc,
         col=RColorBrewer::brewer.pal(10, "Spectral"),
         labels_row = as.expression(lapply(rownames(div_loc), function(a) bquote(italic(.(a))))),
         margins=c(10,10),
         fontsize_col=9,
         cexRow=1.25,
         legend=F,
         treeheight_row = 80,
         treeheight_col=10,
         cellheight = 22.2))
load(file="Data/Summary_GenomeScope.Rdata")
ss=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]
Summary_GenomeScope$Species=factor(Summary_GenomeScope$Species,levels=levels(Summary_GenomeScope$Species)[rev(pheat$tree_row$order)])
color_med_atl=data.frame(Location=c("Li","Mu","Fa","Ga"),
                         Col=brewer.pal(n = 4, name = "RdBu"))
Summary_GenomeScope$Species=factor(Summary_GenomeScope$Species)
color=c()
for (i in 1:nrow(Summary_GenomeScope)){
  color[i]=as.character(color_med_atl$Col[which(as.character(Summary_GenomeScope$Location[i])==color_med_atl$Location)])
}
Summary_GenomeScope$col=color

Summary_GenomeScope$Location=factor(Summary_GenomeScope$Location,levels=c("Li","Mu","Fa","Ga"))

for (i in 1:nrow(Summary_GenomeScope)){
  Summary_GenomeScope$div[i]=as.vector(div)[which(as.character(Summary_GenomeScope$Species[i])==names(div))]
}

Summary_GenomeScope$div=as.numeric(Summary_GenomeScope$div)
p_div_main<-ggplot(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], 
                   aes(x=Species, 
                       y=Heterozygosity,
                       text=paste("Sample:",Sample," \n Heterozygosity:",Heterozygosity," %"))) + 
  geom_point(col=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]$col,
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T) + 
  ylim(c(0,2))+
  geom_segment(aes(x=Species, 
                   xend=Species, 
                   y=0, 
                   yend=div,
                   text1=paste("Species:",Species," \n Heterozygosity:",div," %")))  +
  coord_flip()+
  geom_point(aes(x=Species,y=div),
             col="orange",
             size=3)+
  ylab("Genetic diversity (%)")+
  xlab("") +
  theme_bw()+
  theme(
    axis.text.y=element_blank(),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_line(colour="white"),
    panel.grid.minor = element_line(colour="white"),
    legend.position="bottom") +
  scale_color_identity("",
                       guide="legend",
                       labels = c("Gulf of Lion",
                                  "Murcia",
                                  "Faro",
                                  "Bay of Biscay"),
                       breaks = c("#CA0020",
                                  "#F4A582",
                                  "#92C5DE",
                                  "#0571B0"))

df <- data.frame()
fake_plot<-ggplot(df) + geom_point() + xlim(0,1) + ylim(1,17)
for (j in 1:length(levels(Summary_GenomeScope$Species))){
  
  img<-readPNG(paste("Data/",levels(Summary_GenomeScope$Species)[j],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  fake_plot<-fake_plot+  
    annotation_custom(g, 
                      xmin=0, 
                      xmax=1, 
                      ymin=j,
                      ymax=j+1)
  
}

fake_plot<-fake_plot+
  theme_void()+
  ylim(0.25,16.4)
main_sp<-ggarrange(fake_plot,p_div_main,ncol=2,widths = c(0.1,0.9),heights = c(0.9,1))

gatk_div=read.table("Data/diversity_Dlabr.het",header=T)
hh=(gatk_div$N_SITES-gatk_div$O.HOM.)/gatk_div$N_SITES
indiv_het=hh*(11.5/675)*100
gatk_div=cbind(gatk_div,indiv_het)
geno=Summary_GenomeScope[Summary_GenomeScope$Species=="Dlabr",]
indiv_geno=rep(0,20)

for (i in 1:20){
  
  for (j in 1:20){
    
    if (gatk_div$INDV[j]==geno$Sample[i]){
      
      indiv_geno[i]=gatk_div$indiv_het[j]
      
    }
    
  }
  
}

geno=cbind(geno,indiv_geno)
color=c()
for (i in 1:nrow(geno)){
  color[i]=as.character(color_med_atl$Col[which(as.character(geno$Location[i])==color_med_atl$Location)])
}
geno$col=color
geno$indiv_geno=as.numeric(geno$indiv_geno)
geno$Heterozygosity=as.numeric(geno$Heterozygosity)
my.formula <- y~x
ll<-lm(Heterozygosity~indiv_geno,data=geno)
pvalue=summary(ll)$coefficients[2,4]
summary(ll)

dlabr <- readPNG(paste("Data/Dlabr.png",sep=""))
g <- rasterGrob(dlabr, interpolate=TRUE)

p1<-ggplot(geno,aes(x=indiv_geno,y=Heterozygosity))+
  theme_classic()+
  stat_smooth(method="lm", se=T,col="black",linetype="dashed")+
  geom_point_interactive(
    aes(tooltip = paste("Sample:",Sample,
                        "\n Location:",Location,
                        "\n GenomeScope:",round(Heterozygosity,3),
                        "\n GATK:",round(indiv_geno,3)),data_id = Location,
        color=Location),
    size = 3,
    position='jitter',
    alpha=0.75)+
  annotate("text", x = 0.37, y = 0.34, label = "paste(italic(p), \" = 4.45e-10 \")",parse=T)+
  annotate("text", x = 0.37, y = 0.32, label = "paste(italic(R) ^ 2, \" = 0.89 \")",parse=T)+
  ylab("")+
  xlab("")+
  ggtitle(expression(paste("Sea bass (", italic("D. labrax"), ")")))+
  annotation_custom(g, ymin=0.40, ymax=0.44, xmin=0.29, xmax=0.32) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(breaks = c("Li", "Mu", "Fa","Ga"),
                     labels=c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                     values=as.vector(color_med_atl$Col))

gatk_div=read.table("Data/het_Spilc.het",header=T)
hh=(gatk_div$N_SITES-gatk_div$O.HOM.)/gatk_div$N_SITES
indiv_het=hh*10
gatk_div=cbind(gatk_div,indiv_het)
geno=Summary_GenomeScope[Summary_GenomeScope$Species=="Spilc",]
indiv_geno=rep(0,20)

for (i in 1:20){
  for (j in 1:20){
    if (gatk_div$INDV[j]==geno$Sample[i]){
      indiv_geno[i]=indiv_het[j]
    }
  }
}

geno=cbind(geno,indiv_geno)
color=c()
for (i in 1:nrow(geno)){
  color[i]=as.character(color_med_atl$Col[which(as.character(geno$Location[i])==color_med_atl$Location)])
}
geno$col=color
geno$indiv_geno=as.numeric(geno$indiv_geno)
geno$Heterozygosity=as.numeric(geno$Heterozygosity)
my.formula <- y~x
ll<-lm(Heterozygosity~indiv_geno,data=geno[geno$Heterozygosity<=1.8,])
pvalue=summary(ll)$coefficients[2,4]
summary(ll)

spilc <- readPNG(paste("Data/Spilc.png",sep=""))
g <- rasterGrob(spilc, interpolate=TRUE)
geno$shape=rep(0,20)
for (i in 1:20){
  if (geno$Heterozygosity[i]>=1.8){
    geno$shape[i]=2
  } else {
    geno$shape[i]=1
  }
}
geno$shape=factor(geno$shape)
p2<-ggplot(geno,aes(x=indiv_geno,y=Heterozygosity,color=Location))+
  theme_classic()+
  stat_smooth(data=geno[geno$Heterozygosity<=1.8,],method="lm", se=T,col="black",linetype="dashed")+
  geom_point_interactive(
    aes(tooltip = paste("Sample:",Sample,
                        "\n Location:",Location,
                        "\n GenomeScope:",round(Heterozygosity,3),
                        "\n GATK:",round(indiv_geno,3)),data_id = Location,
        color=Location,
        shape=shape),
    size = 3,
    position='jitter',
    alpha=0.75)+ 
  annotate("text", x = 1.4, y = 1.3, label = "paste(italic(p), \" = 0.0363 \")",parse=T)+
  annotate("text", x = 1.4, y = 1.25, label = "paste(italic(R) ^ 2, \" = 0.246 \")",parse=T)+
  ylab("")+
  xlab("")+
  ggtitle(expression(paste("European pilchard (", italic("S. pilchardus"), ")")))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(breaks = c("Li", "Mu", "Fa","Ga"),
                     labels=c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                     values=as.vector(color_med_atl$Col))+
  annotation_custom(g, ymin=1.8, ymax=2.0, xmin=1.275, xmax=1.3) 

figure1<-ggarrange(p1, p2, nrow = 2, ncol = 1,common.legend = TRUE,legend="top")
figure1<-annotate_figure(figure1,
                         left = text_grob("Genetic diversity estimated with GenomeScope (%)",rot = 90,vjust=2),
                         bottom = text_grob("Genetic diversity estimated after mapping and variant calling (%)",vjust=-1.5)
)


pdf(paste("figures/main.pdf",sep=""),width=11,height=11)
print(ggarrange(map,
                figure1,
                heat,
                main_sp,
                labels=c("A","B","C","D"),
                nrow=2,
                ncol=2,
                vjust=c(1.5,1.5,0,0)
                )
    )
dev.off()


# Figure 3 ----
nodes <- c('0','1','2','3','4','5','6','7','8','9')
x <- seq(0,19,by=2)
y <- rep(90,10)
from <- c('0','1','2','3','4','5','6','7','8',
          '2','3','4','5','6','7','8','9')
to <- c('1','2','3','4','5','6','7','8','9',
        '0','0','0','0','0','0','0','0')
NodeList <- data.frame(nodes, x ,y)
EdgeList <- data.frame(from, to)
a<- graph_from_data_frame(vertices = NodeList, 
                          d= EdgeList, 
                          directed = TRUE)
a <-simplify(a, remove.multiple = FALSE, remove.loops = TRUE)
V(a)$size <- 8
V(a)$frame.color <- "black"
V(a)$color <- "orange"

label=c(expression('S'['0->1']),
        expression('S'['1->2']),
        expression('S'['2->3']),
        expression('S'['3->4']),
        expression('S'['4->5']),
        expression('S'['5->6']),
        expression('S'['6->7']),
        expression('S'['7->8']),
        expression('S'['8->9']),
        expression('f'['2']),
        expression('f'['3']),
        expression('f'['4']),
        expression('f'['5']),
        expression('f'['6']),
        expression('f'['7']),
        expression('f'['8']),
        expression('f'['9']))


L_t=c()
Linf=50
k=0.485
t0=0.44
L=10
for (i in 0:9){
  L_t[i+1]=Linf*(1-exp(-(k*(i-t0))))
}
data_Lt=data.frame(Age=seq(0,9),
                   Lt=L_t)

x=seq(0,L)
b=-L/(-4.60517^(1/1))
Sx=c(1,rep(NA,length(x)-1))
for (i in 2:length(x)){
  Sx[i]=exp(-((x[i]/b)^(1)))
}

age=seq(1,10)
fec=rep(0,10)
Maturity=2
for (j in Maturity:10){
  fec[j]=1*exp((j-2+1)*0.5)
}
fec[2:10]=fec[2:10]/max(fec[2:10])

fec=c(0,fec)
age=c(0,age)

data_Sx=data.frame(Age=rep(x,2),
                   value=c(Sx,fec),
                   lifetime=c(rep("Survival",11),
                              rep("Fecundity",11))
)

raster <- image_read_svg('Data/1146060.svg', width = 350)

run_image=1
if (run_image==1){
  pdf(paste(wd,"/figures/leslie.pdf",sep=""),width=8.25,height=9)
  
  plot(a, edge.arrow.size=.4,
       vertex.color = rgb(0.8,0.4,0.3,0.8),
       vertex.label.family="Times",
       vertex.label.font=2, 
       vertex.label.color="white",
       edge.curved=10,
       frame=FALSE,
       asp=0.075,
       edge.label=label,
       edge.color=c(rep("#56B4E9",9),
                    rep("#E69F00",8)),
       edge.lty=2,
       edge.label.x=c(seq(-0.875,1,by=0.22),
                      seq(-0.89,1.015,by=0.22)[2:9]),
       edge.label.y=c(rep(0,9),
                      #rep(-2.25,8)
                      seq(-2.2,-2.5,length.out=8)),
       edge.label.color=c(rep("#56B4E9",9),
                          rep("#E69F00",8)),
       edge.color = "black",
       edge.label.family="Times",
       edge.label.font=2, 
       edge.label.color="black",
       ylim=c(-15,10))
  
  
  
  y_max=seq(0.5,1.5,length.out=10)+1
  x_decal=lseq(0.025,0.115,length.out=10)
  count=0
  for (j in seq(-1,1,length.out=10)){
    count=count+1
    rasterImage(raster,
                xleft=j-x_decal[count],
                xright=j+x_decal[count], 
                ybottom=1,
                ytop=((x_decal[count]*4)/0.25)+1)
  }
  
  qp <- ggplot(data=data_Lt,aes(x=Age,y=Lt))+
    geom_line(size=1.15)+
    theme_classic()+
    xlab("")+
    ylab("Length")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  sx <- ggplot(data=data_Sx,aes(x=Age,y=value))+
    geom_line(aes(color=lifetime),size=1.25)+
    theme_classic()+
    xlab("")+
    scale_y_continuous(name = "Probability of \n survival \n  until age x", 
                       sec.axis = sec_axis(~./5,breaks=seq(0,1,by=0.1), name = "Relative fecundity", 
                                           labels = seq(0,5,by=0.5))) + 
    scale_x_continuous(breaks=seq(0,10,by=1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank())+
    theme(
      legend.position = c(.6, 1),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6),
      legend.title=element_blank()
    )+
    scale_color_manual(values=c("#E69F00", "#56B4E9"))
  
  print(sx, vp=viewport(0.51, 0.75, .99, .2))
  
  dev.off()
  
  
}
#All plots
load(file="Data/agene/agene_output.Rdata")
lfh<-as.data.frame(read_excel("Data/GENETIC_DIVERSITY_DATA.xlsx",sheet="lfh"))
lfh$div=as.vector(div)
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
load(file="Data/forward_slim/est_species.Rdata")
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
load(file="Data/forward_slim_old/est_species.Rdata")
data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=est_species$Output16,
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


corr_plot<-ggarrange(div_agene,div_slim,labels=c("A","B"),common.legend = T, legend='top')                                    
corr_plot<-annotate_figure(corr_plot,
                           left = text_grob("Observed heterozygosity (%)",rot = 90,vjust=2.5),
)
#pairwise_plot_agene<-ggarrange(pairwise_agene,
#                               pairwise_agene_nopc,
#                               labels=c("C","D"),
#                               nrow=1)
#pairwise_plot_agene<-annotate_figure(pairwise_plot_agene,
#                                     left = text_grob(expression(paste("Ratio of observed ", N[e]/N)),rot = 90,vjust=2),
#)
#pairwise_plot_slim<-ggarrange(pairwise_slim,
#                              pairwise_slim_nopc,
#                              labels=c("E","F"),
#                              nrow=1)
#pairwise_plot_slim<-annotate_figure(pairwise_plot_slim,
#                                    left = text_grob("Ratio of simulated genetic diversity",rot = 90,vjust=2.5),
#)
#pairwise_plot<-ggarrange(pairwise_plot_agene,
#                         pairwise_plot_slim,
#                         nrow=1)
#pairwise_plot<-annotate_figure(pairwise_plot,
#                               bottom=text_grob("Ratio of observed genetic diversity",vjust=-1),
#)

data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=agene_output[[4]]$Output1,
                     x2=lfh$Parental_Care)
data_plot=data_plot[data_plot$x2=="No",]
data_plot$y=data_plot$y/(max(data_plot$y))
data_plot$x=data_plot$x/(max(data_plot$x))

summary(lm(y~x,data=data_plot[data_plot$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

data_plot_tmp=data.frame(species=lfh$Species_plot,
                         y=lfh$div,
                         x=agene_output[[4]]$Output16,
                         x2=lfh$Parental_Care)
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$y=data_plot_tmp$y/(max(data_plot_tmp$y))
data_plot_tmp$x=data_plot_tmp$x/(max(data_plot_tmp$x))

summary(lm(y~x,data=data_plot_tmp[data_plot_tmp$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot_tmp,link='logit')

data_plot_tmp$mean=fitted(m1)
data_plot_tmp$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

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
data_plot_whole$gg=c(rep(1,11),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

#data_plot_whole=data_plot_whole[data_plot_whole$x2=="No" & data_plot_whole$gg=="1",]

#tt_tmp=data.frame(species=lfh$Species_plot,
#                     y=lfh$div,
#☻                     x=agene_output[[4]]$Output16,
#                     x2=lfh$Parental_Care)
#tt_tmp$y=tt_tmp$y/(max(tt_tmp$y))
#tt_tmp$x=tt_tmp$x/(max(tt_tmp$x))
#tt_tmp=tt_tmp[tt_tmp$x2=="No",]
#tt_tmp=cbind(tt_tmp,data_plot_whole[,c(5:8)])
#tt_tmp$gg=rep(2,nrow(tt_tmp))

#data_plot_whole=rbind(data_plot_whole,tt_tmp)

agene_p1<-ggplot(data_plot_whole, 
                  aes(x = x, 
                      y = mean*100,
                      label=species,
                      group=gg
                      #group=x2)
                  )) +
  geom_line(aes(alpha=gg,
    colour=gg),
            size=2)+
  #geom_ribbon(aes(ymin=sd_025,
  #                ymax=sd_975,colour=gg),
  #            alpha=0.05,
  #            size=0.05) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 alpha=gg,
                 color=gg),
             shape=19)+
  geom_rangeframe()+
  theme_classic()+
  #geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  #scale_x_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  #scale_y_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  xlab("")+
  ylab("")+
  ggtitle("--")+
  xlim(c(0,1))+
  ylim(c(0,1))+
  scale_color_manual(name = "Simulation",values=c(viridis(100)[c(47.5,72.5,97.5)][2],"grey"),labels=c("Focal model","Complete model"))+
  scale_alpha_manual(name="Simulation",values=c(1,0.5),labels=c("Focal model","Complete model"))+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept=0,slope=1,col="red",lty=2)
  #annotate("text", x = 0.2025, y = 1.325, label = "paste(italic(p), \" = 0.000966 \")",parse=T)+
  #annotate("text", x = 0.2025, y = 1.265, label = "paste(italic(R) ^ 2, \" = 0.55 \")",parse=T)

data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=agene_output[[4]]$Output2,
                     x2=lfh$Parental_Care)
data_plot=data_plot[data_plot$x2=="No",]
data_plot$y=data_plot$y/(max(data_plot$y))
data_plot$x=data_plot$x/(max(data_plot$x))

summary(lm(y~x,data=data_plot[data_plot$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

data_plot_tmp=data.frame(species=lfh$Species_plot,
                         y=lfh$div,
                         x=agene_output[[4]]$Output16,
                         x2=lfh$Parental_Care)
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$y=data_plot_tmp$y/(max(data_plot_tmp$y))
data_plot_tmp$x=data_plot_tmp$x/(max(data_plot_tmp$x))

summary(lm(y~x,data=data_plot_tmp[data_plot_tmp$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot_tmp,link='logit')

data_plot_tmp$mean=fitted(m1)
data_plot_tmp$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

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
data_plot_whole$gg=c(rep(1,11),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

#data_plot_whole=data_plot_whole[data_plot_whole$x2=="No" & data_plot_whole$gg=="1",]

#tt_tmp=data.frame(species=lfh$Species_plot,
#                     y=lfh$div,
#☻                     x=agene_output[[4]]$Output16,
#                     x2=lfh$Parental_Care)
#tt_tmp$y=tt_tmp$y/(max(tt_tmp$y))
#tt_tmp$x=tt_tmp$x/(max(tt_tmp$x))
#tt_tmp=tt_tmp[tt_tmp$x2=="No",]
#tt_tmp=cbind(tt_tmp,data_plot_whole[,c(5:8)])
#tt_tmp$gg=rep(2,nrow(tt_tmp))

#data_plot_whole=rbind(data_plot_whole,tt_tmp)

agene_p2<-ggplot(data_plot_whole, 
                 aes(x = x, 
                     y = mean*100,
                     label=species,
                     group=gg
                     #group=x2)
                 )) +
  geom_line(aes(alpha=gg,
                colour=gg),
            size=2)+
  #geom_ribbon(aes(ymin=sd_025,
  #                ymax=sd_975,colour=gg),
  #            alpha=0.05,
  #            size=0.05) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 alpha=gg,
                 color=gg),
             shape=19)+
  geom_rangeframe()+
  theme_classic()+
  #geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  #scale_x_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  #scale_y_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  xlab("")+
  ylab("")+
  ggtitle("Age at maturity")+
  xlim(c(0,1))+
  ylim(c(0,1))+
  scale_color_manual(name = "Simulation",values=c(viridis(100)[c(47.5,72.5,97.5)][2],"grey"),labels=c("Focal model","Complete model"))+
  scale_alpha_manual(name="Simulation",values=c(1,0.5),labels=c("Focal model","Complete model"))+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept=0,slope=1,col="red",lty=2)
#annotate("text", x = 0.2025, y = 1.325, label = "paste(italic(p), \" = 0.000966 \")",parse=T)+
#annotate("text", x = 0.2025, y = 1.265, label = "paste(italic(R) ^ 2, \" = 0.55 \")",parse=T)

data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=agene_output[[4]]$Output3,
                     x2=lfh$Parental_Care)
data_plot=data_plot[data_plot$x2=="No",]
data_plot$y=data_plot$y/(max(data_plot$y))
data_plot$x=data_plot$x/(max(data_plot$x))

summary(lm(y~x,data=data_plot[data_plot$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

data_plot_tmp=data.frame(species=lfh$Species_plot,
                         y=lfh$div,
                         x=agene_output[[4]]$Output16,
                         x2=lfh$Parental_Care)
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$y=data_plot_tmp$y/(max(data_plot_tmp$y))
data_plot_tmp$x=data_plot_tmp$x/(max(data_plot_tmp$x))

summary(lm(y~x,data=data_plot_tmp[data_plot_tmp$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot_tmp,link='logit')

data_plot_tmp$mean=fitted(m1)
data_plot_tmp$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

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
data_plot_whole$gg=c(rep(1,11),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

#data_plot_whole=data_plot_whole[data_plot_whole$x2=="No" & data_plot_whole$gg=="1",]

#tt_tmp=data.frame(species=lfh$Species_plot,
#                     y=lfh$div,
#☻                     x=agene_output[[4]]$Output16,
#                     x2=lfh$Parental_Care)
#tt_tmp$y=tt_tmp$y/(max(tt_tmp$y))
#tt_tmp$x=tt_tmp$x/(max(tt_tmp$x))
#tt_tmp=tt_tmp[tt_tmp$x2=="No",]
#tt_tmp=cbind(tt_tmp,data_plot_whole[,c(5:8)])
#tt_tmp$gg=rep(2,nrow(tt_tmp))

#data_plot_whole=rbind(data_plot_whole,tt_tmp)

agene_p3<-ggplot(data_plot_whole, 
                 aes(x = x, 
                     y = mean*100,
                     label=species,
                     group=gg
                     #group=x2)
                 )) +
  geom_line(aes(alpha=gg,
                colour=gg),
            size=2)+
  #geom_ribbon(aes(ymin=sd_025,
  #                ymax=sd_975,colour=gg),
  #            alpha=0.05,
  #            size=0.05) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 alpha=gg,
                 color=gg),
             shape=19)+
  geom_rangeframe()+
  theme_classic()+
  #geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  #scale_x_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  #scale_y_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  xlab("")+
  ylab("")+
  ggtitle("Survival")+
  xlim(c(0,1))+
  ylim(c(0,1))+
  scale_color_manual(name = "Simulation",values=c(viridis(100)[c(47.5,72.5,97.5)][2],"grey"),labels=c("Focal model","Complete model"))+
  scale_alpha_manual(name="Simulation",values=c(1,0.5),labels=c("Focal model","Complete model"))+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept=0,slope=1,col="red",lty=2)
#annotate("text", x = 0.2025, y = 1.325, label = "paste(italic(p), \" = 0.000966 \")",parse=T)+
#annotate("text", x = 0.2025, y = 1.265, label = "paste(italic(R) ^ 2, \" = 0.55 \")",parse=T)

data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=agene_output[[4]]$Output4,
                     x2=lfh$Parental_Care)
data_plot=data_plot[data_plot$x2=="No",]
data_plot$y=data_plot$y/(max(data_plot$y))
data_plot$x=data_plot$x/(max(data_plot$x))

summary(lm(y~x,data=data_plot[data_plot$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

data_plot_tmp=data.frame(species=lfh$Species_plot,
                         y=lfh$div,
                         x=agene_output[[4]]$Output16,
                         x2=lfh$Parental_Care)
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$y=data_plot_tmp$y/(max(data_plot_tmp$y))
data_plot_tmp$x=data_plot_tmp$x/(max(data_plot_tmp$x))

summary(lm(y~x,data=data_plot_tmp[data_plot_tmp$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot_tmp,link='logit')

data_plot_tmp$mean=fitted(m1)
data_plot_tmp$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

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
data_plot_whole$gg=c(rep(1,11),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

#data_plot_whole=data_plot_whole[data_plot_whole$x2=="No" & data_plot_whole$gg=="1",]

#tt_tmp=data.frame(species=lfh$Species_plot,
#                     y=lfh$div,
#☻                     x=agene_output[[4]]$Output16,
#                     x2=lfh$Parental_Care)
#tt_tmp$y=tt_tmp$y/(max(tt_tmp$y))
#tt_tmp$x=tt_tmp$x/(max(tt_tmp$x))
#tt_tmp=tt_tmp[tt_tmp$x2=="No",]
#tt_tmp=cbind(tt_tmp,data_plot_whole[,c(5:8)])
#tt_tmp$gg=rep(2,nrow(tt_tmp))

#data_plot_whole=rbind(data_plot_whole,tt_tmp)

agene_p4<-ggplot(data_plot_whole, 
                 aes(x = x, 
                     y = mean*100,
                     label=species,
                     group=gg
                     #group=x2)
                 )) +
  geom_line(aes(alpha=gg,
                colour=gg),
            size=2)+
  #geom_ribbon(aes(ymin=sd_025,
  #                ymax=sd_975,colour=gg),
  #            alpha=0.05,
  #            size=0.05) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 alpha=gg,
                 color=gg),
             shape=19)+
  geom_rangeframe()+
  theme_classic()+
  #geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  #scale_x_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  #scale_y_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  xlab("")+
  ylab("")+
  ggtitle("Age at Maturity + Survival")+
  xlim(c(0,1))+
  ylim(c(0,1))+
  scale_color_manual(name = "Simulation",values=c(viridis(100)[c(47.5,72.5,97.5)][2],"grey"),labels=c("Focal model","Complete model"))+
  scale_alpha_manual(name="Simulation",values=c(1,0.5),labels=c("Focal model","Complete model"))+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept=0,slope=1,col="red",lty=2)
#annotate("text", x = 0.2025, y = 1.325, label = "paste(italic(p), \" = 0.000966 \")",parse=T)+
#annotate("text", x = 0.2025, y = 1.265, label = "paste(italic(R) ^ 2, \" = 0.55 \")",parse=T)


data_plot=data.frame(species=lfh$Species_plot,
                     y=lfh$div,
                     x=agene_output[[4]]$Output5,
                     x2=lfh$Parental_Care)
data_plot=data_plot[data_plot$x2=="No",]
data_plot$y=data_plot$y/(max(data_plot$y))
data_plot$x=data_plot$x/(max(data_plot$x))

summary(lm(y~x,data=data_plot[data_plot$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot,link='logit')
m1_nopc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="No",],link='logit')
m1_pc <- betareg(I(y/100)~x,data=data_plot[data_plot$x2=="Yes",],link='logit')

data_plot$mean=fitted(m1)
data_plot$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

data_plot_tmp=data.frame(species=lfh$Species_plot,
                         y=lfh$div,
                         x=agene_output[[4]]$Output16,
                         x2=lfh$Parental_Care)
data_plot_tmp=data_plot_tmp[data_plot_tmp$x2=="No",]
data_plot_tmp$y=data_plot_tmp$y/(max(data_plot_tmp$y))
data_plot_tmp$x=data_plot_tmp$x/(max(data_plot_tmp$x))

summary(lm(y~x,data=data_plot_tmp[data_plot_tmp$x2=="No",]))

m1<-betareg(I(y/100)~x,data=data_plot_tmp,link='logit')

data_plot_tmp$mean=fitted(m1)
data_plot_tmp$sd_025=fitted(m1)*100-2*coefficients(summary(m1))$mean[2,2]
data_plot_tmp$sd_975=fitted(m1)*100+2*coefficients(summary(m1))$mean[2,2]

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
data_plot_whole$gg=c(rep(1,11),rep(2,11))
data_plot_whole$gg=factor(data_plot_whole$gg)

#data_plot_whole=data_plot_whole[data_plot_whole$x2=="No" & data_plot_whole$gg=="1",]

#tt_tmp=data.frame(species=lfh$Species_plot,
#                     y=lfh$div,
#☻                     x=agene_output[[4]]$Output16,
#                     x2=lfh$Parental_Care)
#tt_tmp$y=tt_tmp$y/(max(tt_tmp$y))
#tt_tmp$x=tt_tmp$x/(max(tt_tmp$x))
#tt_tmp=tt_tmp[tt_tmp$x2=="No",]
#tt_tmp=cbind(tt_tmp,data_plot_whole[,c(5:8)])
#tt_tmp$gg=rep(2,nrow(tt_tmp))

#data_plot_whole=rbind(data_plot_whole,tt_tmp)

agene_p5<-ggplot(data_plot_whole, 
                 aes(x = x, 
                     y = mean*100,
                     label=species,
                     group=gg
                     #group=x2)
                 )) +
  geom_line(aes(alpha=gg,
                colour=gg),
            size=2)+
  #geom_ribbon(aes(ymin=sd_025,
  #                ymax=sd_975,colour=gg),
  #            alpha=0.05,
  #            size=0.05) +
  geom_point(data=data_plot_whole,size=2.5,
             aes(x=x,
                 y=y,
                 alpha=gg,
                 color=gg),
             shape=19)+
  geom_rangeframe()+
  theme_classic()+
  #geom_text_repel(data=data_plot,aes(x=x,y=y),label=italic_species,col="black",parse=T)+
  #scale_x_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  #scale_y_continuous(breaks=seq(0,1,by=0.1),
  #                   labels=as.character(seq(0,1,by=0.1)))+
  xlab("")+
  ylab("")+
  ggtitle("Fecundity")+
  xlim(c(0,1))+
  ylim(c(0,1))+
  scale_color_manual(name = "Simulation",values=c(viridis(100)[c(47.5,72.5,97.5)][2],"grey"),labels=c("Focal model","Complete model"))+
  scale_alpha_manual(name="Simulation",values=c(1,0.5),labels=c("Focal model","Complete model"))+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))+
  geom_abline(intercept=0,slope=1,col="red",lty=2)
#annotate("text", x = 0.2025, y = 1.325, label = "paste(italic(p), \" = 0.000966 \")",parse=T)+
#annotate("text", x = 0.2025, y = 1.265, label = "paste(italic(R) ^ 2, \" = 0.55 \")",parse=T)

agene_unique_plot<-ggarrange(agene_p1,
                               agene_p2,
                               agene_p3,
                               agene_p4,
                               agene_p5,
labels=c("C","D","E","F","G"),
nrow=1,
common.legend = T)

agene_unique_plot<-annotate_figure(agene_unique_plot,
left = text_grob("Scaled genetic diversity",rot = 90,vjust=1),
bottom = text_grob(paste("Scaled",expression(Ne/N),sep=" "),vjust=-0.5))
                 

all_plot<-ggarrange(
  corr_plot,
agene_unique_plot,
  nrow=2)

#width=12,height=7.5
pdf(paste(wd,"/figures/vk_plot.pdf",sep=""),width=15,height=9.375)
print(all_plot)
dev.off()

#surv_fec<-ggarrange(
#  surv_plot,
#  fecun_plot,
#  nrow=2)
#pdf(paste(wd,"/figures/sur_fec.pdf",sep=""),width=12,height=7.5)
#print(surv_fec)
#dev.off()
