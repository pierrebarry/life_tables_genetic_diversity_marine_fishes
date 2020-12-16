#------------------------------------------------------------------------#
#                                                                        #
#                         GENOMESCOPE                                    #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages ----
library(ggpubr)
library("RColorBrewer")
library(ggplot2)
library(png)
library(grid)
library(plotly)
library(ggiraph)
# Load data ----
load(file="Data/Summary_GenomeScope.Rdata")
load(file="Data/div.Rdata")
# GenomeScope ----
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

p_div_genomescope<-ggplotly(p_div_main,tooltip=c("text","text1"),height = 650)  %>% partial_bundle() 
img1<-readPNG(paste("Data/","Cgale",".png",sep=""))
img2<-readPNG(paste("Data/","Cjuli",".png",sep=""))
img3<-readPNG(paste("Data/","Dlabr",".png",sep=""))
img4<-readPNG(paste("Data/","Dpunt",".png",sep=""))
img5<-readPNG(paste("Data/","Hgutt",".png",sep=""))
img6<-readPNG(paste("Data/","Lbude",".png",sep=""))
img7<-readPNG(paste("Data/","Lmorm",".png",sep=""))
img8<-readPNG(paste("Data/","Mmerl",".png",sep=""))
img9<-readPNG(paste("Data/","Msurm",".png",sep=""))
img10<-readPNG(paste("Data/","Peryt",".png",sep=""))
img11<-readPNG(paste("Data/","Scabr",".png",sep=""))
img12<-readPNG(paste("Data/","Scant",".png",sep=""))
img13<-readPNG(paste("Data/","Scine",".png",sep=""))
img14<-readPNG(paste("Data/","Spilc",".png",sep=""))
img15<-readPNG(paste("Data/","Ssard",".png",sep=""))
img16<-readPNG(paste("Data/","Styph",".png",sep=""))
j=1
p_div_genomescope<-  p_div_genomescope %>%
  layout(
    images = list(
      list(
        source = raster2uri(as.raster(img1)),
        x = -0.001, y = j/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img2)),
        x = -0.001, y = (j+1)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img3)),
        x = -0.001, y = (j+2)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img4)),
        x = -0.001, y = (j+3)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img5)),
        x = 0, y = (j+4)/16, 
        sizex = 0.5, sizey = 0.05
      ),
      list(
        source = raster2uri(as.raster(img6)),
        x = -0.001, y = (j+5)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img7)),
        x = -0.001, y = (j+6)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img8)),
        x = -0.001, y = (j+7)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img9)),
        x = -0.001, y = (j+8)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img10)),
        x = -0.001, y = (j+9)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img11)),
        x = -0.001, y = (j+10)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img12)),
        x = -0.001, y = (j+11)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img13)),
        x = -0.001, y = (j+12)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img14)),
        x = -0.001, y = (j+13)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img15)),
        x = -0.001, y = (j+14)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img16)),
        x = -0.001, y = (j+15-0.5)/16, 
        sizex = 0.1, sizey = 0.1
      )
    )
  ) 

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
p_div_main<-p_div_main+ylim(c(0,2))
p_div_main<-ggarrange(fake_plot,p_div_main,ncol=2,widths = c(0.1,0.9),heights = c(0.9,1))
pdf(paste("figures/Diversity_lollipop.pdf",sep=""),width=7.5,height=5)
p_div_main
dev.off()

p<-ggplot(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], 
          aes(x=Species, 
              y=Heterozygosity)) + 
  geom_point(aes(col=Location),
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T)+
  geom_boxplot(outlier.shape = NA,
               alpha=0)+
  geom_violin(color="black",
              alpha=0) +
  ylab("Heterozygosity (%)")+
  ylim(c(min(Summary_GenomeScope$Heterozygosity),
         max(Summary_GenomeScope$Heterozygosity)+0.1*max(Summary_GenomeScope$Heterozygosity)))+
  xlab("Species") +
  scale_x_discrete(labels=c("Coryphoblennius \n galerita",
                            "Coris \n julis",
                            "Dicentrarchus \n labrax",
                            "Diplodus \n puntazzo",
                            "Hippocampus \n guttulatus",
                            "Lophius \n budegassa",
                            "Lithognathus \n mormyrus",
                            "Merluccius \n merluccius",
                            "Mullus \n surmuletus",
                            "Pagellus \n erythrinus",
                            "Serranus \n cabrilla",
                            "Spondyliosoma \n cantharus",
                            "Symphodus \n cinereus",
                            "Sardina \n pilchardus",
                            "Sarda \n sarda",
                            "Syngnathus \n typhle"
  ))+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_line(colour="white")) +
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(legend.position="bottom")

for (j in 1:length(levels(Summary_GenomeScope$Species))){
  
  img<-readPNG(paste("Data/",levels(Summary_GenomeScope$Species)[j],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom(g, xmin=j-0.4, 
                      xmax=j+0.4, 
                      ymin=max(Summary_GenomeScope$Heterozygosity)+0.01*max(Summary_GenomeScope$Heterozygosity),
                      ymax=max(Summary_GenomeScope$Heterozygosity)+0.1*max(Summary_GenomeScope$Heterozygosity)) 
  
}


pdf(paste("figures/Diversity.pdf",sep=""),width=7.5,height=5)
p
dev.off()

# Genome length plot ----
genome_length=with(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], tapply(Genome_Length, Species, mean))
Summary_GenomeScope$genome_length=rep(0,nrow(Summary_GenomeScope))
for (i in 1:nrow(Summary_GenomeScope)){
  Summary_GenomeScope$genome_length[i]=as.numeric(genome_length[which(names(genome_length)==Summary_GenomeScope$Species[i])])
}
p_div_main<-ggplot(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], 
                   aes(x=Species, 
                       y=Genome_Length,
                       text=paste("Sample:",Sample," \n Genome length:",Genome_Length/1000000," Mb"))) + 
  geom_point(col=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]$col,
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T) + 
  ylim(c(-5e7,8e+8))+
  geom_segment(aes(x=Species, 
                   xend=Species, 
                   y=0, 
                   yend=genome_length,
                   text1=paste("Species:",Species," \n Genome length:",Genome_Length/1000000," Mp")))  +
  coord_flip()+
  geom_point(aes(x=Species,y=genome_length),
             col="orange",
             size=3)+
  ylab("Genome length (bp)")+
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
p_genome_length_genomescope<-ggplotly(p_div_main,tooltip=c("text","text1"),height = 650)  %>% partial_bundle() 

img1<-readPNG(paste("Data/","Cgale",".png",sep=""))
img2<-readPNG(paste("Data/","Cjuli",".png",sep=""))
img3<-readPNG(paste("Data/","Dlabr",".png",sep=""))
img4<-readPNG(paste("Data/","Dpunt",".png",sep=""))
img5<-readPNG(paste("Data/","Hgutt",".png",sep=""))
img6<-readPNG(paste("Data/","Lbude",".png",sep=""))
img7<-readPNG(paste("Data/","Lmorm",".png",sep=""))
img8<-readPNG(paste("Data/","Mmerl",".png",sep=""))
img9<-readPNG(paste("Data/","Msurm",".png",sep=""))
img10<-readPNG(paste("Data/","Peryt",".png",sep=""))
img11<-readPNG(paste("Data/","Scabr",".png",sep=""))
img12<-readPNG(paste("Data/","Scant",".png",sep=""))
img13<-readPNG(paste("Data/","Scine",".png",sep=""))
img14<-readPNG(paste("Data/","Spilc",".png",sep=""))
img15<-readPNG(paste("Data/","Ssard",".png",sep=""))
img16<-readPNG(paste("Data/","Styph",".png",sep=""))

j=1
p_genome_length_genomescope<-  p_genome_length_genomescope %>%
  layout(
    images = list(
      list(
        source = raster2uri(as.raster(img1)),
        x = -0.001, y = j/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img2)),
        x = -0.001, y = (j+1)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img3)),
        x = -0.001, y = (j+2)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img4)),
        x = -0.001, y = (j+3)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img5)),
        x = 0, y = (j+4)/16, 
        sizex = 0.5, sizey = 0.05
      ),
      list(
        source = raster2uri(as.raster(img6)),
        x = -0.001, y = (j+5)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img7)),
        x = -0.001, y = (j+6)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img8)),
        x = -0.001, y = (j+7)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img9)),
        x = -0.001, y = (j+8)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img10)),
        x = -0.001, y = (j+9)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img11)),
        x = -0.001, y = (j+10)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img12)),
        x = -0.001, y = (j+11)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img13)),
        x = -0.001, y = (j+12)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img14)),
        x = -0.001, y = (j+13)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img15)),
        x = -0.001, y = (j+14)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img16)),
        x = -0.001, y = (j+15-0.5)/16, 
        sizex = 0.1, sizey = 0.1
      )
    )
  ) 



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
p_genome_length_main<-p_div_main+ylim(c(0,8e+8))
p_genome_length_main<-ggarrange(fake_plot,p_genome_length_main,ncol=2,widths = c(0.1,0.9),heights = c(0.9,1))
pdf(paste("figures/Genome_length_lollipop.pdf",sep=""),width=7.5,height=5)
p_genome_length_main
dev.off()

p<-ggplot(Summary_GenomeScope, 
          aes(x=Species, 
              y=Genome_Length)) + 
  geom_point(aes(col=Location),
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T)+
  geom_boxplot(outlier.shape = NA,
               alpha=0)+
  geom_violin(color="black",
              alpha=0) +
  ylab("Genome Length (bp)")+
  ylim(c(min(Summary_GenomeScope$Genome_Length),
         max(Summary_GenomeScope$Genome_Length)+0.1*max(Summary_GenomeScope$Genome_Length)))+
  xlab("Species") +
  scale_x_discrete(labels=c("Coryphoblennius \n galerita",
                            "Coris \n julis",
                            "Dicentrarchus \n labrax",
                            "Diplodus \n puntazzo",
                            "Hippocampus \n guttulatus",
                            "Lophius \n budegassa",
                            "Lithognathus \n mormyrus",
                            "Merluccius \n merluccius",
                            "Mullus \n surmuletus",
                            "Pagellus \n erythrinus",
                            "Serranus \n cabrilla",
                            "Spondyliosoma \n cantharus",
                            "Symphodus \n cinereus",
                            "Sardina \n pilchardus",
                            "Sarda \n sarda",
                            "Syngnathus \n typhle"
  ))+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_line(colour="white")) +
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(legend.position="bottom")

for (j in 1:length(levels(Summary_GenomeScope$Species))){
  
  img<-readPNG(paste("Data/",levels(Summary_GenomeScope$Species)[j],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom(g, xmin=j-0.4, 
                      xmax=j+0.4, 
                      ymin=max(Summary_GenomeScope$Genome_Length)+0.01*max(Summary_GenomeScope$Genome_Length),
                      ymax=max(Summary_GenomeScope$Genome_Length)+0.1*max(Summary_GenomeScope$Genome_Length)) 
  
}


pdf(paste("figures/Genome_length.pdf",sep=""),width=7.5,height=5)
p
dev.off()




# Genome repeat plot ----
genome_repeat=with(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], tapply(Genome_Repeat, Species, mean))
Summary_GenomeScope$genome_repeat=rep(0,nrow(Summary_GenomeScope))
for (i in 1:nrow(Summary_GenomeScope)){
  Summary_GenomeScope$genome_repeat[i]=as.numeric(genome_repeat[which(names(genome_repeat)==Summary_GenomeScope$Species[i])])
}
p_div_main<-ggplot(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], 
                   aes(x=Species, 
                       y=Genome_Repeat,
                       text=paste("Sample:",Sample," \n Genome repeat length:",Genome_Repeat/1000000," Mb"))) + 
  geom_point(col=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]$col,
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T) + 
  ylim(c(-2e7,2.5e+8))+
  geom_segment(aes(x=Species, 
                   xend=Species, 
                   y=0, 
                   yend=genome_repeat,
                   text1=paste("Species:",Species," \n Genome repeat length:",Genome_Repeat/1000000," Mp")))  +
  coord_flip()+
  geom_point(aes(x=Species,y=genome_repeat),
             col="orange",
             size=3)+
  ylab("Genome repeat length (bp)")+
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
p_genome_repeat_genomescope<-ggplotly(p_div_main,tooltip=c("text","text1"),height = 650)  %>% partial_bundle() 

img1<-readPNG(paste("Data/","Cgale",".png",sep=""))
img2<-readPNG(paste("Data/","Cjuli",".png",sep=""))
img3<-readPNG(paste("Data/","Dlabr",".png",sep=""))
img4<-readPNG(paste("Data/","Dpunt",".png",sep=""))
img5<-readPNG(paste("Data/","Hgutt",".png",sep=""))
img6<-readPNG(paste("Data/","Lbude",".png",sep=""))
img7<-readPNG(paste("Data/","Lmorm",".png",sep=""))
img8<-readPNG(paste("Data/","Mmerl",".png",sep=""))
img9<-readPNG(paste("Data/","Msurm",".png",sep=""))
img10<-readPNG(paste("Data/","Peryt",".png",sep=""))
img11<-readPNG(paste("Data/","Scabr",".png",sep=""))
img12<-readPNG(paste("Data/","Scant",".png",sep=""))
img13<-readPNG(paste("Data/","Scine",".png",sep=""))
img14<-readPNG(paste("Data/","Spilc",".png",sep=""))
img15<-readPNG(paste("Data/","Ssard",".png",sep=""))
img16<-readPNG(paste("Data/","Styph",".png",sep=""))

j=1
p_genome_repeat_genomescope<-  p_genome_repeat_genomescope %>%
  layout(
    images = list(
      list(
        source = raster2uri(as.raster(img1)),
        x = -0.001, y = j/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img2)),
        x = -0.001, y = (j+1)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img3)),
        x = -0.001, y = (j+2)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img4)),
        x = -0.001, y = (j+3)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img5)),
        x = 0, y = (j+4)/16, 
        sizex = 0.5, sizey = 0.05
      ),
      list(
        source = raster2uri(as.raster(img6)),
        x = -0.001, y = (j+5)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img7)),
        x = -0.001, y = (j+6)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img8)),
        x = -0.001, y = (j+7)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img9)),
        x = -0.001, y = (j+8)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img10)),
        x = -0.001, y = (j+9)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img11)),
        x = -0.001, y = (j+10)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img12)),
        x = -0.001, y = (j+11)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img13)),
        x = -0.001, y = (j+12)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img14)),
        x = -0.001, y = (j+13)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img15)),
        x = -0.001, y = (j+14)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img16)),
        x = -0.001, y = (j+15-0.5)/16, 
        sizex = 0.1, sizey = 0.1
      )
    )
  ) 



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
p_genome_repeat<-p_div_main+ylim(c(0,2.5e+8))
p_genome_repeat<-ggarrange(fake_plot,p_genome_repeat,ncol=2,widths = c(0.1,0.9),heights = c(0.9,1))
pdf(paste("figures/Genome_repeat_lollipop.pdf",sep=""),width=7.5,height=5)
p_genome_repeat
dev.off()

p<-ggplot(Summary_GenomeScope, 
          aes(x=Species, 
              y=Genome_Repeat)) + 
  geom_point(aes(col=Location),
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T)+
  geom_boxplot(outlier.shape = NA,
               alpha=0)+
  geom_violin(color="black",
              alpha=0) +
  ylab("Genome repeat length (bp)")+
  ylim(c(min(Summary_GenomeScope$Genome_Repeat),
         max(Summary_GenomeScope$Genome_Repeat)+0.1*max(Summary_GenomeScope$Genome_Repeat)))+
  xlab("Species") +
  scale_x_discrete(labels=c("Coryphoblennius \n galerita",
                            "Coris \n julis",
                            "Dicentrarchus \n labrax",
                            "Diplodus \n puntazzo",
                            "Hippocampus \n guttulatus",
                            "Lophius \n budegassa",
                            "Lithognathus \n mormyrus",
                            "Merluccius \n merluccius",
                            "Mullus \n surmuletus",
                            "Pagellus \n erythrinus",
                            "Serranus \n cabrilla",
                            "Spondyliosoma \n cantharus",
                            "Symphodus \n cinereus",
                            "Sardina \n pilchardus",
                            "Sarda \n sarda",
                            "Syngnathus \n typhle"
  ))+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_line(colour="white")) +
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(legend.position="bottom")

for (j in 1:length(levels(Summary_GenomeScope$Species))){
  
  img<-readPNG(paste("Data/",levels(Summary_GenomeScope$Species)[j],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom(g, xmin=j-0.4, 
                      xmax=j+0.4, 
                      ymin=max(Summary_GenomeScope$Genome_Repeat)+0.01*max(Summary_GenomeScope$Genome_Repeat),
                      ymax=max(Summary_GenomeScope$Genome_Repeat)+0.1*max(Summary_GenomeScope$Genome_Repeat)) 
  
}


pdf(paste("figures/Genome_repeat.pdf",sep=""),width=7.5,height=5)
p
dev.off()




# Genome unique plot ----
genome_unique=with(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], tapply(Genome_Unique, Species, mean))
Summary_GenomeScope$genome_unique=rep(0,nrow(Summary_GenomeScope))
for (i in 1:nrow(Summary_GenomeScope)){
  Summary_GenomeScope$genome_unique[i]=as.numeric(genome_unique[which(names(genome_unique)==Summary_GenomeScope$Species[i])])
}
p_div_main<-ggplot(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], 
                   aes(x=Species, 
                       y=Genome_Unique,
                       text=paste("Sample:",Sample," \n Genome repeat length:",Genome_Unique/1000000," Mb"))) + 
  geom_point(col=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]$col,
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T) + 
  ylim(c(-5e7,6.5e+8))+
  geom_segment(aes(x=Species, 
                   xend=Species, 
                   y=0, 
                   yend=genome_unique,
                   text1=paste("Species:",Species," \n Genome repeat length:",Genome_Unique/1000000," Mp")))  +
  coord_flip()+
  geom_point(aes(x=Species,y=genome_unique),
             col="orange",
             size=3)+
  ylab("Genome unique length (bp)")+
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
p_genome_unique_genomescope<-ggplotly(p_div_main,tooltip=c("text","text1"),height = 650)  %>% partial_bundle() 

img1<-readPNG(paste("Data/","Cgale",".png",sep=""))
img2<-readPNG(paste("Data/","Cjuli",".png",sep=""))
img3<-readPNG(paste("Data/","Dlabr",".png",sep=""))
img4<-readPNG(paste("Data/","Dpunt",".png",sep=""))
img5<-readPNG(paste("Data/","Hgutt",".png",sep=""))
img6<-readPNG(paste("Data/","Lbude",".png",sep=""))
img7<-readPNG(paste("Data/","Lmorm",".png",sep=""))
img8<-readPNG(paste("Data/","Mmerl",".png",sep=""))
img9<-readPNG(paste("Data/","Msurm",".png",sep=""))
img10<-readPNG(paste("Data/","Peryt",".png",sep=""))
img11<-readPNG(paste("Data/","Scabr",".png",sep=""))
img12<-readPNG(paste("Data/","Scant",".png",sep=""))
img13<-readPNG(paste("Data/","Scine",".png",sep=""))
img14<-readPNG(paste("Data/","Spilc",".png",sep=""))
img15<-readPNG(paste("Data/","Ssard",".png",sep=""))
img16<-readPNG(paste("Data/","Styph",".png",sep=""))

j=1
p_genome_unique_genomescope<-  p_genome_unique_genomescope %>%
  layout(
    images = list(
      list(
        source = raster2uri(as.raster(img1)),
        x = -0.001, y = j/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img2)),
        x = -0.001, y = (j+1)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img3)),
        x = -0.001, y = (j+2)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img4)),
        x = -0.001, y = (j+3)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img5)),
        x = 0, y = (j+4)/16, 
        sizex = 0.5, sizey = 0.05
      ),
      list(
        source = raster2uri(as.raster(img6)),
        x = -0.001, y = (j+5)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img7)),
        x = -0.001, y = (j+6)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img8)),
        x = -0.001, y = (j+7)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img9)),
        x = -0.001, y = (j+8)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img10)),
        x = -0.001, y = (j+9)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img11)),
        x = -0.001, y = (j+10)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img12)),
        x = -0.001, y = (j+11)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img13)),
        x = -0.001, y = (j+12)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img14)),
        x = -0.001, y = (j+13)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img15)),
        x = -0.001, y = (j+14)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img16)),
        x = -0.001, y = (j+15-0.5)/16, 
        sizex = 0.1, sizey = 0.1
      )
    )
  ) 



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
p_genome_unique<-p_div_main+ylim(c(0,6.5e+8))
p_genome_unique<-ggarrange(fake_plot,p_genome_unique,ncol=2,widths = c(0.1,0.9),heights = c(0.9,1))
pdf(paste("figures/Genome_unique_lollipop.pdf",sep=""),width=7.5,height=5)
p_genome_unique
dev.off()

p<-ggplot(Summary_GenomeScope, 
          aes(x=Species, 
              y=Genome_Unique)) + 
  geom_point(aes(col=Location),
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T)+
  geom_boxplot(outlier.shape = NA,
               alpha=0)+
  geom_violin(color="black",
              alpha=0) +
  ylab("Genome unique length (bp)")+
  ylim(c(min(Summary_GenomeScope$Genome_Unique),
         max(Summary_GenomeScope$Genome_Unique)+0.1*max(Summary_GenomeScope$Genome_Unique)))+
  xlab("Species") +
  scale_x_discrete(labels=c("Coryphoblennius \n galerita",
                            "Coris \n julis",
                            "Dicentrarchus \n labrax",
                            "Diplodus \n puntazzo",
                            "Hippocampus \n guttulatus",
                            "Lophius \n budegassa",
                            "Lithognathus \n mormyrus",
                            "Merluccius \n merluccius",
                            "Mullus \n surmuletus",
                            "Pagellus \n erythrinus",
                            "Serranus \n cabrilla",
                            "Spondyliosoma \n cantharus",
                            "Symphodus \n cinereus",
                            "Sardina \n pilchardus",
                            "Sarda \n sarda",
                            "Syngnathus \n typhle"
  ))+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_line(colour="white")) +
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(legend.position="bottom")

for (j in 1:length(levels(Summary_GenomeScope$Species))){
  
  img<-readPNG(paste("Data/",levels(Summary_GenomeScope$Species)[j],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom(g, xmin=j-0.4, 
                      xmax=j+0.4, 
                      ymin=max(Summary_GenomeScope$Genome_Unique)+0.01*max(Summary_GenomeScope$Genome_Unique),
                      ymax=max(Summary_GenomeScope$Genome_Unique)+0.1*max(Summary_GenomeScope$Genome_Unique)) 
  
}


pdf(paste("figures/Genome_unique.pdf",sep=""),width=7.5,height=5)
p
dev.off()

# Model fit plot ----
model_fit=with(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], tapply(Model_Fit, Species, mean))
Summary_GenomeScope$model_fit=rep(0,nrow(Summary_GenomeScope))
for (i in 1:nrow(Summary_GenomeScope)){
  Summary_GenomeScope$model_fit[i]=as.numeric(model_fit[which(names(model_fit)==Summary_GenomeScope$Species[i])])
}
p_div_main<-ggplot(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], 
                   aes(x=Species, 
                       y=Model_Fit,
                       text=paste("Sample:",Sample," \n Model fit:",Model_Fit," %"))) + 
  geom_point(col=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]$col,
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T) + 
  ylim(c(-5,100))+
  geom_segment(aes(x=Species, 
                   xend=Species, 
                   y=0, 
                   yend=model_fit,
                   text1=paste("Species:",Species," \n Model fit:",Model_Fit," %")))  +
  coord_flip()+
  geom_point(aes(x=Species,y=model_fit),
             col="orange",
             size=3)+
  ylab("Model fit (%)")+
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
p_model_fit_genomescope<-ggplotly(p_div_main,tooltip=c("text","text1"),height = 650)  %>% partial_bundle() 

img1<-readPNG(paste("Data/","Cgale",".png",sep=""))
img2<-readPNG(paste("Data/","Cjuli",".png",sep=""))
img3<-readPNG(paste("Data/","Dlabr",".png",sep=""))
img4<-readPNG(paste("Data/","Dpunt",".png",sep=""))
img5<-readPNG(paste("Data/","Hgutt",".png",sep=""))
img6<-readPNG(paste("Data/","Lbude",".png",sep=""))
img7<-readPNG(paste("Data/","Lmorm",".png",sep=""))
img8<-readPNG(paste("Data/","Mmerl",".png",sep=""))
img9<-readPNG(paste("Data/","Msurm",".png",sep=""))
img10<-readPNG(paste("Data/","Peryt",".png",sep=""))
img11<-readPNG(paste("Data/","Scabr",".png",sep=""))
img12<-readPNG(paste("Data/","Scant",".png",sep=""))
img13<-readPNG(paste("Data/","Scine",".png",sep=""))
img14<-readPNG(paste("Data/","Spilc",".png",sep=""))
img15<-readPNG(paste("Data/","Ssard",".png",sep=""))
img16<-readPNG(paste("Data/","Styph",".png",sep=""))

j=1
p_model_fit_genomescope<-  p_model_fit_genomescope %>%
  layout(
    images = list(
      list(
        source = raster2uri(as.raster(img1)),
        x = -0.001, y = j/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img2)),
        x = -0.001, y = (j+1)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img3)),
        x = -0.001, y = (j+2)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img4)),
        x = -0.001, y = (j+3)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img5)),
        x = 0, y = (j+4)/16, 
        sizex = 0.5, sizey = 0.05
      ),
      list(
        source = raster2uri(as.raster(img6)),
        x = -0.001, y = (j+5)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img7)),
        x = -0.001, y = (j+6)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img8)),
        x = -0.001, y = (j+7)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img9)),
        x = -0.001, y = (j+8)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img10)),
        x = -0.001, y = (j+9)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img11)),
        x = -0.001, y = (j+10)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img12)),
        x = -0.001, y = (j+11)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img13)),
        x = -0.001, y = (j+12)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img14)),
        x = -0.001, y = (j+13)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img15)),
        x = -0.001, y = (j+14)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img16)),
        x = -0.001, y = (j+15-0.5)/16, 
        sizex = 0.1, sizey = 0.1
      )
    )
  ) 



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
p_model_fit<-p_div_main+ylim(c(0,100))
p_model_fit<-ggarrange(fake_plot,p_model_fit,ncol=2,widths = c(0.1,0.9),heights = c(0.9,1))
pdf(paste("figures/Model_fit_lollipop.pdf",sep=""),width=7.5,height=5)
p_model_fit
dev.off()

p<-ggplot(Summary_GenomeScope, 
          aes(x=Species, 
              y=Model_Fit)) + 
  geom_point(aes(col=Location),
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T)+
  geom_boxplot(outlier.shape = NA,
               alpha=0)+
  geom_violin(color="black",
              alpha=0) +
  ylab("Model fit (%)")+
  ylim(c(min(Summary_GenomeScope$Model_Fit),
         max(Summary_GenomeScope$Model_Fit)+0.1*max(Summary_GenomeScope$Model_Fit)))+
  xlab("Species") +
  scale_x_discrete(labels=c("Coryphoblennius \n galerita",
                            "Coris \n julis",
                            "Dicentrarchus \n labrax",
                            "Diplodus \n puntazzo",
                            "Hippocampus \n guttulatus",
                            "Lophius \n budegassa",
                            "Lithognathus \n mormyrus",
                            "Merluccius \n merluccius",
                            "Mullus \n surmuletus",
                            "Pagellus \n erythrinus",
                            "Serranus \n cabrilla",
                            "Spondyliosoma \n cantharus",
                            "Symphodus \n cinereus",
                            "Sardina \n pilchardus",
                            "Sarda \n sarda",
                            "Syngnathus \n typhle"
  ))+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_line(colour="white")) +
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(legend.position="bottom")

for (j in 1:length(levels(Summary_GenomeScope$Species))){
  
  img<-readPNG(paste("Data/",levels(Summary_GenomeScope$Species)[j],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom(g, xmin=j-0.4, 
                      xmax=j+0.4, 
                      ymin=max(Summary_GenomeScope$Model_Fit)+0.01*max(Summary_GenomeScope$Model_Fit),
                      ymax=max(Summary_GenomeScope$Model_Fit)+0.1*max(Summary_GenomeScope$Model_Fit)) 
  
}


pdf(paste("figures/Model_fit.pdf",sep=""),width=7.5,height=5)
p
dev.off()



# Read error plot ----
read_error=with(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], tapply(Read_error_rate, Species, mean))
Summary_GenomeScope$read_error=rep(0,nrow(Summary_GenomeScope))
for (i in 1:nrow(Summary_GenomeScope)){
  Summary_GenomeScope$read_error[i]=as.numeric(read_error[which(names(read_error)==Summary_GenomeScope$Species[i])])
}
p_div_main<-ggplot(Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",], 
                   aes(x=Species, 
                       y=Read_error_rate,
                       text=paste("Sample:",Sample," \n Read error rate:",Read_error_rate*100," %"))) + 
  geom_point(col=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]$col,
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T) + 
  ylim(c(-0.05,1))+
  geom_segment(aes(x=Species, 
                   xend=Species, 
                   y=0, 
                   yend=read_error,
                   text1=paste("Species:",Species," \n Read error rate:",Read_error_rate*100," %")))  +
  coord_flip()+
  geom_point(aes(x=Species,y=read_error),
             col="orange",
             size=3)+
  ylab("Read error rate (%)")+
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
p_read_error_genomescope<-ggplotly(p_div_main,tooltip=c("text","text1"),height = 650)  %>% partial_bundle() 

img1<-readPNG(paste("Data/","Cgale",".png",sep=""))
img2<-readPNG(paste("Data/","Cjuli",".png",sep=""))
img3<-readPNG(paste("Data/","Dlabr",".png",sep=""))
img4<-readPNG(paste("Data/","Dpunt",".png",sep=""))
img5<-readPNG(paste("Data/","Hgutt",".png",sep=""))
img6<-readPNG(paste("Data/","Lbude",".png",sep=""))
img7<-readPNG(paste("Data/","Lmorm",".png",sep=""))
img8<-readPNG(paste("Data/","Mmerl",".png",sep=""))
img9<-readPNG(paste("Data/","Msurm",".png",sep=""))
img10<-readPNG(paste("Data/","Peryt",".png",sep=""))
img11<-readPNG(paste("Data/","Scabr",".png",sep=""))
img12<-readPNG(paste("Data/","Scant",".png",sep=""))
img13<-readPNG(paste("Data/","Scine",".png",sep=""))
img14<-readPNG(paste("Data/","Spilc",".png",sep=""))
img15<-readPNG(paste("Data/","Ssard",".png",sep=""))
img16<-readPNG(paste("Data/","Styph",".png",sep=""))

j=1
p_read_error_genomescope<-  p_read_error_genomescope %>%
  layout(
    images = list(
      list(
        source = raster2uri(as.raster(img1)),
        x = -0.001, y = j/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img2)),
        x = -0.001, y = (j+1)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img3)),
        x = -0.001, y = (j+2)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img4)),
        x = -0.001, y = (j+3)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img5)),
        x = 0, y = (j+4)/16, 
        sizex = 0.5, sizey = 0.05
      ),
      list(
        source = raster2uri(as.raster(img6)),
        x = -0.001, y = (j+5)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img7)),
        x = -0.001, y = (j+6)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img8)),
        x = -0.001, y = (j+7)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img9)),
        x = -0.001, y = (j+8)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img10)),
        x = -0.001, y = (j+9)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img11)),
        x = -0.001, y = (j+10)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img12)),
        x = -0.001, y = (j+11)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img13)),
        x = -0.001, y = (j+12)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img14)),
        x = -0.001, y = (j+13)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img15)),
        x = -0.001, y = (j+14)/16, 
        sizex = 0.1, sizey = 0.1
      ),
      list(
        source = raster2uri(as.raster(img16)),
        x = -0.001, y = (j+15-0.5)/16, 
        sizex = 0.1, sizey = 0.1
      )
    )
  ) 



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
p_read_error_rate<-p_div_main+ylim(c(0,1))
p_read_error_rate<-ggarrange(fake_plot,p_read_error_rate,ncol=2,widths = c(0.1,0.9),heights = c(0.9,1))
pdf(paste("figures/Read_error_lollipop.pdf",sep=""),width=7.5,height=5)
p_read_error_rate
dev.off()

p<-ggplot(Summary_GenomeScope, 
          aes(x=Species, 
              y=Read_error_rate)) + 
  geom_point(aes(col=Location),
             size = 1.5,
             position="jitter",
             alpha=0.75,
             show.legend = T)+
  geom_boxplot(outlier.shape = NA,
               alpha=0)+
  geom_violin(color="black",
              alpha=0) +
  ylab("Read error rate (%)")+
  ylim(c(min(Summary_GenomeScope$Read_error_rate),
         max(Summary_GenomeScope$Read_error_rate)+0.1*max(Summary_GenomeScope$Read_error_rate)))+
  xlab("Species") +
  scale_x_discrete(labels=c("Coryphoblennius \n galerita",
                            "Coris \n julis",
                            "Dicentrarchus \n labrax",
                            "Diplodus \n puntazzo",
                            "Hippocampus \n guttulatus",
                            "Lophius \n budegassa",
                            "Lithognathus \n mormyrus",
                            "Merluccius \n merluccius",
                            "Mullus \n surmuletus",
                            "Pagellus \n erythrinus",
                            "Serranus \n cabrilla",
                            "Spondyliosoma \n cantharus",
                            "Symphodus \n cinereus",
                            "Sardina \n pilchardus",
                            "Sarda \n sarda",
                            "Syngnathus \n typhle"
  ))+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_line(colour="white")) +
  scale_colour_manual(name = "",
                      labels = c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                      values = brewer.pal(n = 4, name = "RdBu"))+
  theme(legend.position="bottom")

for (j in 1:length(levels(Summary_GenomeScope$Species))){
  
  img<-readPNG(paste("Data/",levels(Summary_GenomeScope$Species)[j],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom(g, xmin=j-0.4, 
                      xmax=j+0.4, 
                      ymin=max(Summary_GenomeScope$Read_error_rate)+0.01*max(Summary_GenomeScope$Read_error_rate),
                      ymax=max(Summary_GenomeScope$Read_error_rate)+0.1*max(Summary_GenomeScope$Read_error_rate)) 
  
}


pdf(paste("figures/Read_error.pdf",sep=""),width=7.5,height=5)
p
dev.off()

# Supp figure ----
pdf(paste("figures/genomescope_suppmat.pdf",sep=""),width=12.5,height=12.5)
ggarrange(p_genome_length_main,
          p_genome_unique,
          p_genome_repeat,
          p_model_fit,
          p_read_error_rate,
          labels=c("A","B","C","D","E"),
          label.y=1.015,
          nrow=3,ncol=2)
dev.off()

# Verification of GenomeScope ------
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

# Plotly
p1<-ggplot(geno,aes(x=indiv_geno,y=Heterozygosity))+
  theme_classic()+
  stat_smooth(method="lm", se=T,col="black",linetype="dashed")+
  geom_point(aes(color=Location,
                 text=paste("Sample:",Sample," \n Heterozygosity - GenomeScope:",Heterozygosity," % \n GATK :",indiv_geno)),
             size=3,
             position="jitter",
             alpha=0.75)+
  annotate("text", x = 0.37, y = 0.34, label = "paste(italic(p), \" = 4.45e-10 \")",parse=T)+
  annotate("text", x = 0.37, y = 0.32, label = "paste(italic(R) ^ 2, \" = 0.89 \")",parse=T)+
  ylab("")+
  xlab("")+
  annotation_custom(g, ymin=0.40, ymax=0.44, xmin=0.29, xmax=0.32) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(breaks = c("Li", "Mu", "Fa","Ga"),
                     labels=c("Gulf of Lion","Mar Menor","Algarve","Bay of Biscay"),
                     values=as.vector(color_med_atl$Col))
p_GATK_Dlabr<-ggplotly(p1,tooltip=c("text"),height = 250)  %>% partial_bundle() 

p_GATK_Dlabr<-  p_GATK_Dlabr %>%
  layout(
    images = list(
      list(
        source = raster2uri(as.raster(img3)),
        y = 0.3, x = 0.425, 
        sizex = 0.5, sizey = 0.5
      )))

# Plot
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

# Spilc
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

# Plotly
p2<-ggplot(geno,aes(x=indiv_geno,y=Heterozygosity))+
  theme_classic()+
  stat_smooth(data=geno[geno$Heterozygosity<=1.8,],method="lm", se=T,col="black",linetype="dashed")+
  geom_point(aes(color=Location,
                 text=paste("Sample:",Sample," \n Heterozygosity - GenomeScope:",Heterozygosity," % \n GATK :",indiv_geno)),
             size=3,
             position="jitter",
             alpha=0.75)+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(breaks = c("Li", "Mu", "Fa","Ga"),
                     labels=c("Gulf of Lion","Mar Menor","Algarve","Bay of Biscay"),
                     values=as.vector(color_med_atl$Col))
p_GATK_Spilc<-ggplotly(p2,tooltip=c("text"),height = 250)  %>% partial_bundle() 

p_GATK_Spilc<-  p_GATK_Spilc %>%
  layout(
    images = list(
      list(
        source = raster2uri(as.raster(img14)),
        y = 1.8, x = 1.275, 
        sizex = 0.1, sizey = 0.1
      )))
# Plot
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


gatk_fig=figure1
pdf(paste("figures/GenomeScope_GATK.pdf",sep=""),width=10,height=7.5)
print(gatk_fig)
dev.off()
# GenomeScope sensibility ----
load(file="Data/Sensibility_test.Rdata")

p<- ggplot(Sensibility_test, aes(x=Kmer, y=Hetero)) + 
  geom_line() +
  theme_classic()+
  xlab("Length of k-mer")+
  ylab("Heterozygosity (%)")+
  geom_pointrange(aes(ymin=Hetero_min,ymax=Hetero_max),size=0.5)+
  geom_vline(xintercept=21,linetype="dashed", color = "red")

fig_sensib=p
pdf(paste("figures/Sensibility_GenomeScope.pdf",sep=""),width=10,height=5)
print(fig_sensib)
dev.off()

# Export genetic diversity values ----
ss=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]
round(tapply(ss$Heterozygosity,ss$Species,FUN=median),3)
round(tapply(ss$Heterozygosity,ss$Species,FUN=mean),3)
round(tapply(ss$Heterozygosity,ss$Species,FUN=sd),3)
div<-tapply(ss$Heterozygosity,ss$Species,FUN=median)
genome_length<-tapply(ss$Genome_Length,ss$Species,FUN=median)

save(div,file="Data/div.Rdata")
save(genome_length,file="Data/genome_length.Rdata")

