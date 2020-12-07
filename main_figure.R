#------------------------------------------------------------------------#
#                                                                        #
#                           FIGURES                                      #
#                                                                        #
#------------------------------------------------------------------------#

# Figure 1 -----
library(heatmap3)
library(ggplotify)
library(pheatmap)

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

heat<-as.ggplot(pheat<-pheatmap(div_loc,
         col=RColorBrewer::brewer.pal(10, "Spectral"),
         #col=plasma(500),
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
# Genetic diversity plot ----

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
#♠  ylim(0.5,16.25)
main_sp<-ggarrange(fake_plot,p_div_main,ncol=2,widths = c(0.1,0.9),heights = c(0.9,1))



wd="C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT"
setwd(wd)
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
library(grid)
library(png)
library(emdbook)
library(TeachingDemos)
library(igraph)
library(magick)
library(ggplot2)

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
#fec[Maturity[i]:Lifespan[i]]=fec[Maturity[i]:Lifespan[i]]+abs(min(fec[Maturity[i]:Lifespan[i]]))
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
corr_plot<-ggarrange(div_agene,div_slim,labels=c("A","B"),common.legend = T, legend='top')                                    
corr_plot<-annotate_figure(corr_plot,
                           left = text_grob("Observed heterozygosity (%)",rot = 90,vjust=2.5),
)
pairwise_plot_agene<-ggarrange(pairwise_agene,
                               pairwise_agene_nopc,
                               labels=c("C","D"),
                               nrow=1)
pairwise_plot_agene<-annotate_figure(pairwise_plot_agene,
                                     left = text_grob(expression(paste("Ratio of observed ", N[e]/N)),rot = 90,vjust=2),
)
pairwise_plot_slim<-ggarrange(pairwise_slim,
                              pairwise_slim_nopc,
                              labels=c("E","F"),
                              nrow=1)
pairwise_plot_slim<-annotate_figure(pairwise_plot_slim,
                                    left = text_grob("Ratio of simulated genetic diversity",rot = 90,vjust=2.5),
)
pairwise_plot<-ggarrange(pairwise_plot_agene,
                         pairwise_plot_slim,
                         nrow=1)
pairwise_plot<-annotate_figure(pairwise_plot,
                               bottom=text_grob("Ratio of observed genetic diversity",vjust=-1),
)
all_plot<-ggarrange(
  corr_plot,
  pairwise_plot,
  nrow=2)

pdf(paste(wd,"/figures/vk_plot.pdf",sep=""),width=12,height=7.5)
print(all_plot)
dev.off()

surv_fec<-ggarrange(
  surv_plot,
  fecun_plot,
  nrow=2)
pdf(paste(wd,"/figures/sur_fec.pdf",sep=""),width=12,height=7.5)
print(surv_fec)
dev.off()
