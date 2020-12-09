#------------------------------------------------------------------------#
#                                                                        #
#                     Mapping statistics with samtools                   #
#                                                                        #
#------------------------------------------------------------------------#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load packages -----
library(tidyverse)
library(plotly)
library(RColorBrewer)
library(png)
library(grid)
library(cowplot)
annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){ layer(data = data, stat = StatIdentity, position = PositionIdentity, 
                                                                                 geom = ggplot2:::GeomCustomAnn,
                                                                                 inherit.aes = TRUE, params = list(grob = grob, 
                                                                                                                   xmin = xmin, xmax = xmax, 
                                                                                                                   ymin = ymin, ymax = ymax))}
# Collect data ----
SAMPLES=c("DlabrLi1","DlabrLi2","DlabrLi3","DlabrLi4","DlabrLi5",
          "DlabrMu1","DlabrMu2","DlabrMu3","DlabrMu4","DlabrMu6",
          "DlabrFa1","DlabrFa3","DlabrFa4","DlabrFa5","DlabrFa6",
          "DlabrGa2","DlabrGa3","DlabrGa4","DlabrGa5","DlabrGa6",
          "SpilcLi2","SpilcLi3","SpilcLi4","SpilcLi5","SpilcLi6",
          "SpilcMu1","SpilcMu2","SpilcMu3","SpilcMu4","SpilcMu6",
          "SpilcFa1","SpilcFa3","SpilcFa4","SpilcFa5","SpilcFa6",
          "SpilcGa1","SpilcGa3","SpilcGa4","SpilcGa5","SpilcGa6")
SPECIES=c(rep("Dlabr",20),
          rep("Spilc",20))



color_med_atl=data.frame(Location=c("Li","Mu","Fa","Ga"),
                         Col=brewer.pal(n = 4, name = "RdBu"))

if (length(SAMPLES)==0 | length(SPECIES)==0){
  print("No samples or species provided")
} else {
  
  # Loop ----
  globalstats=data.frame(SAMPLES=c(NA),
                         SPECIES=c(NA),
                         LOCATION=c(NA),
                         stats=c(NA),
                         value=c(NA),
                         Color=c(NA))
  
  readlength_data=data.frame(SAMPLES=c(NA),
                             SPECIES=c(NA),
                             LOCATION=c(NA),
                             read_length=c(NA),
                             Count=c(NA),
                             Color=c(NA))
  
  insertsize_data=data.frame(SAMPLES=c(NA),
                             SPECIES=c(NA),
                             LOCATION=c(NA),
                             insert_size=c(NA),
                             allpairs=c(NA),
                             Color=c(NA))
  
  coverage_data=data.frame(SAMPLES=c(NA),
                           SPECIES=c(NA),
                           LOCATION=c(NA),
                           Coverage=c(NA),
                           Count=c(NA),
                           Color=c(NA))
  
  gc_coverage_data=data.frame(SAMPLES=c(NA),
                              SPECIES=c(NA),
                              LOCATION=c(NA),
                              GC=c(NA),
                              Count=c(NA),
                              Color=c(NA))
  
  insertion_data=data.frame(SAMPLES=c(NA),
                            SPECIES=c(NA),
                            LOCATION=c(NA),
                            insertion=c(NA),
                            count=c(NA),
                            Color=c(NA))
  
  deletion_data=data.frame(SAMPLES=c(NA),
                           SPECIES=c(NA),
                           LOCATION=c(NA),
                           deletion=c(NA),
                           count=c(NA),
                           Color=c(NA))
  
  for (i in 1:length(SAMPLES)){
    
    flagstat<-readLines(paste("Data/samtools/",SAMPLES[i],"_flagstat.txt",sep=""))
    stats<- readLines(paste("Data/samtools/",SAMPLES[i],"_samtools_stats.stats",sep=""))
    
    # GLOBAL STATS
    sn <- grep("^SN",stats, value=TRUE)
    sn <- separate(data.frame(sn),col=1, into=c("ID", "Name","Value"), sep="\t")[,-1]
    for (j in 1:nrow(sn)){
      globalstats=rbind(globalstats,
                        c(SAMPLES[i],
                          substr(SAMPLES[i],0,5),
                          substr(SAMPLES[i],6,7),
                          sn$Name[j],
                          sn$Value[j],
                          as.character(color_med_atl[which(substr(SAMPLES[i],6,7)==color_med_atl$Location),2])))
    }
    
    # READ LENGTH
    rl <- grep("^RL",stats, value=TRUE)
    rl <- separate(data.frame(rl),col=1, into=c("ID", "read_length", "count"), sep="\t")[,-1]
    for (j in 1:length(rl$count)){
      readlength_data=rbind(readlength_data,
                            c(SAMPLES[i],
                              substr(SAMPLES[i],0,5),
                              substr(SAMPLES[i],6,7),
                              rl$read_length[j],
                              rl$count[j],
                              as.character(color_med_atl[which(substr(SAMPLES[i],6,7)==color_med_atl$Location),2])))
    }
    
    
    # INSERT SIZE
    is <- grep("^IS",stats, value=TRUE)
    is <- separate(data.frame(is),col=1, into=c("ID", "insert size","all pairs", "inward", "outward", "other"), sep="\t")[,-1]
    for (j in 1:length(is$`insert size`)){
      insertsize_data=rbind(insertsize_data,
                            c(SAMPLES[i],
                              substr(SAMPLES[i],0,5),
                              substr(SAMPLES[i],6,7),
                              is$`insert size`[j],
                              is$`all pairs`[j],
                              as.character(color_med_atl[which(substr(SAMPLES[i],6,7)==color_med_atl$Location),2])))
    }
    
    
    # COVERAGE PLOT
    cov <- grep("^COV",stats, value=TRUE)
    cov <- separate(data.frame(cov),col=1, into=c("ID","value","Coverage", "count"), sep="\t")[,-c(1,2)]
    for (j in 1:length(cov$count)){
      coverage_data=rbind(coverage_data,
                          c(SAMPLES[i],
                            substr(SAMPLES[i],0,5),
                            substr(SAMPLES[i],6,7),
                            cov$Coverage[j],
                            cov$count[j],
                            as.character(color_med_atl[which(substr(SAMPLES[i],6,7)==color_med_atl$Location),2])))
    }
    
    # GC Coverage
    gcd <- grep("^GCD",stats, value=TRUE)
    gcd <- separate(data.frame(gcd),col=1, into=c("ID","GC","Unique_Sequence","10th","25th","50th","75th","90th"), sep="\t")[,-1]
    for (j in 1:length(gcd$GC)){
      gc_coverage_data=rbind(gc_coverage_data,
                             c(SAMPLES[i],
                               substr(SAMPLES[i],0,5),
                               substr(SAMPLES[i],6,7),
                               gcd$GC[j],
                               gcd$Unique_Sequence[j],
                               as.character(color_med_atl[which(substr(SAMPLES[i],6,7)==color_med_atl$Location),2])))
    }
    
    #INSERTION DATA
    # Indel distribution
    id <- grep("^ID",stats, value=TRUE)
    id <- separate(data.frame(id),col=1, into=c("ID", "length", "insertion_count", "deletion_count"), sep="\t")[,-1]
    for (j in 1:length(id$length)){
      insertion_data=rbind(insertion_data,
                           c(SAMPLES[i],
                             substr(SAMPLES[i],0,5),
                             substr(SAMPLES[i],6,7),
                             id$length[j],
                             id$insertion_count[j],
                             as.character(color_med_atl[which(substr(SAMPLES[i],6,7)==color_med_atl$Location),2])))
    }
    # Deletion distribution
    for (j in 1:length(id$length)){
      deletion_data=rbind(deletion_data,
                          c(SAMPLES[i],
                            substr(SAMPLES[i],0,5),
                            substr(SAMPLES[i],6,7),
                            id$length[j],
                            id$deletion_count[j],
                            as.character(color_med_atl[which(substr(SAMPLES[i],6,7)==color_med_atl$Location),2])))
    }
    
    
  }
  
  # Plot statistics ----
  #GLOBAL STATS
  globalstats=globalstats[-1,]
  globalstats$value=as.numeric(globalstats$value)
  
  for (i in 1:length(SAMPLES)){
    inter=globalstats[(globalstats$stats=="reads mapped:" | 
                   globalstats$stats=="reads unmapped:") & globalstats$SAMPLES==SAMPLES[i],]
  }
  
  # Reads mapped and unmapped
  p<-ggplot(data=globalstats[globalstats$stats=="reads mapped:" | 
                               globalstats$stats=="reads unmapped:",], aes(x=SAMPLES,y=value,fill=stats)) 
  
  x=levels(as.factor(unique(globalstats$SPECIES)))[1]
  img<-readPNG(paste("Data/",x,".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom2(g, xmin=18.5, xmax=20.5, ymin=1.3e+8, ymax=1.55e+8, data=globalstats[globalstats$stats=="reads mapped:" | 
                                                                                         globalstats$stats=="reads unmapped:",][1,])
  x=levels(as.factor(unique(globalstats$SPECIES)))[2]
  img<-readPNG(paste("Data/",x,".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)  
  
  p<-p+  
    annotation_custom2(g, xmin=18.5, xmax=20.5, ymin=1.3e+8, ymax=1.55e+8, data=globalstats[globalstats$stats=="reads mapped:" | 
                                                                                         globalstats$stats=="reads unmapped:",][41,])
  
  globalstats$SPECIES_plot=factor(globalstats$SPECIES,labels=c("italic('Dicentrarchus labrax')",
                                                                   "italic('Sardina pilchardus')"))

  p<-p+
    geom_bar(stat="identity",alpha=0.25)+
    ggtitle("Reads mapped - unmapped")+
    xlab("Samples")+
    ylab("Number of reads")+
    facet_grid(SPECIES_plot~.,scales='free',labeller=label_parsed)+
    coord_flip()+
    geom_vline(xintercept=seq(5.5,15.5,by=5),linetype="dashed",col="chartreuse4")+
    scale_fill_manual(values = c("blue","red"),name="")+
    theme_classic()+
    ylim(c(0,1.5e+8))+
    theme(strip.text.x = element_text(size = 8)) 

  pdf(paste(wd,"/figures/samtools.pdf",sep=""),width=7.5,height=7.5)
  print(p)
  dev.off()
  
  p <- ggplotly(p)
  
  # Reads mapped
  p<-ggplot(data=globalstats[globalstats$stats=="reads mapped:",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="reads mapped:",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Reads mapped")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value)*0.01,
           (max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)))+
    ggtitle("Reads mapped")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
     
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))

    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
  
  }
  
  
  pdf(paste(wd,"/figures/reads_mapped.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p
  
  # Reads mapped and pairds
  p<-ggplot(data=globalstats[globalstats$stats=="reads mapped and paired:",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="reads mapped and paired:",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Reads mapped and paired")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="reads mapped and paired:",]$value)-min(globalstats[globalstats$stats=="reads mapped and paired:",]$value)*0.01,
           (max(globalstats[globalstats$stats=="reads mapped and paired:",]$value)-min(globalstats[globalstats$stats=="reads mapped and paired:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)))+
    ggtitle("Reads mapped and paired")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))
    
    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
    
  }
  
  
  pdf(paste(wd,"/figures/reads_mapped_paired.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p
  
  # Reads unmapped
  p<-ggplot(data=globalstats[globalstats$stats=="reads unmapped:",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="reads unmapped:",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Reads unmapped")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="reads unmapped:",]$value)-min(globalstats[globalstats$stats=="reads unmapped:",]$value)*0.01,
           (max(globalstats[globalstats$stats=="reads unmapped:",]$value)-min(globalstats[globalstats$stats=="reads unmapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads unmapped:",]$value)))+
    ggtitle("Reads unmapped")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))
    
    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
    
  }
  
  
  pdf(paste(wd,"/figures/reads_unmapped.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p 
  
  # Reads properly paired
  p<-ggplot(data=globalstats[globalstats$stats=="reads properly paired:",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="reads properly paired:",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Reads properly paired")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="reads properly paired:",]$value)-min(globalstats[globalstats$stats=="reads properly paired:",]$value)*0.01,
           (max(globalstats[globalstats$stats=="reads properly paired:",]$value)-min(globalstats[globalstats$stats=="reads properly paired:",]$value))*0.1+max(globalstats[globalstats$stats=="reads properly paired:",]$value)))+
    ggtitle("Reads properly paired")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))
    
    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
    
  }
  
  
  pdf(paste(wd,"/figures/reads_properly_paired.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p 
  
  # Average length
  p<-ggplot(data=globalstats[globalstats$stats=="average length:",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="average length:",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Average length")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="average length:",]$value)-min(globalstats[globalstats$stats=="average length:",]$value)*0.01,
           (max(globalstats[globalstats$stats=="average length:",]$value)-min(globalstats[globalstats$stats=="average length:",]$value))*0.1+max(globalstats[globalstats$stats=="average length:",]$value)))+
    ggtitle("Average length")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))
    
    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
    
  }
  
  
  pdf(paste(wd,"/figures/average_length.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p 
  
  # Average quality
  p<-ggplot(data=globalstats[globalstats$stats=="average quality:",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="average quality:",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Average quality")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="average quality:",]$value)-min(globalstats[globalstats$stats=="average quality:",]$value)*0.01,
           (max(globalstats[globalstats$stats=="average quality:",]$value)-min(globalstats[globalstats$stats=="average quality:",]$value))*0.1+max(globalstats[globalstats$stats=="average quality:",]$value)))+
    ggtitle("Average quality")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))
    
    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
    
  }
  
  
  pdf(paste(wd,"/figures/average_quality.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p 
  
  # Insert size average
  p<-ggplot(data=globalstats[globalstats$stats=="insert size average:",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="insert size average:",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Insert size average")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="insert size average:",]$value)-min(globalstats[globalstats$stats=="insert size average:",]$value)*0.01,
           (max(globalstats[globalstats$stats=="insert size average:",]$value)-min(globalstats[globalstats$stats=="insert size average:",]$value))*0.1+max(globalstats[globalstats$stats=="insert size average:",]$value)))+
    ggtitle("Insert size average")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))
    
    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
    
  }
  
  
  pdf(paste(wd,"/figures/insert_size_average.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p 
  
  # Pairs on different chromosomes
  p<-ggplot(data=globalstats[globalstats$stats=="pairs on different chromosomes:",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="pairs on different chromosomes:",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Pair on different chromosomes")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value)-min(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value)*0.01,
           (max(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value)-min(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value))*0.1+max(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value)))+
    ggtitle("Pair on different chromosomes")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))
    
    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
    
  }
  
  
  pdf(paste(wd,"/figures/pair_diff_chromo.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p 
  
  # Percentage of properly paired reads
  p<-ggplot(data=globalstats[globalstats$stats=="percentage of properly paired reads (%):",], aes(x=SPECIES, y=value)) +
    geom_boxplot(outlier.shape = NA,alpha=0.75)+
    geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$Color,aes(text=SAMPLES))+
    theme_bw()+
    ylab("Percentage of properly paired reads")+
    xlab("Species")+
    ylim(c(min(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value)-min(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value)*0.01,
           (max(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value)-min(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value))*0.1+max(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value)))+
    ggtitle("Percentage of properly paired reads")+
    scale_x_discrete(labels=c(
      "Dicentrarchus \n labrax",
      "Sardina \n pilchardus"
    ))+
    theme(axis.text.x = element_text(face = "italic"))
  
  pimage <- axis_canvas(p, axis = 'x')
  
  for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    
    x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    img<-readPNG(paste("Data/",x,".png",sep=""))
    
    pimage <- pimage +
      draw_image(img, x = c(0.5,1.5)[j], scale = 0.5)
    
  }
  
  
  pdf(paste(wd,"/figures/percentage_properly_paired.pdf",sep=""),width=7.5,height=7.5)
  print(ggdraw(insert_xaxis_grob(p, pimage, position = "bottom")))
  dev.off()
  
  p <- ggplotly(p)
  p 
  
  # READ LENGTH
  readlength_data=readlength_data[-1,]
  readlength_data$read_length=as.numeric(readlength_data$read_length)
  readlength_data$Count=as.numeric(readlength_data$Count)
  p <- ggplot(readlength_data, aes(x=read_length, y=Count,group=SAMPLES)) +
    geom_line(color=readlength_data$Color)+
    theme_classic()+
    facet_grid(SPECIES ~ LOCATION)
  pdf(paste(wd,"/figures/read_length.pdf",sep=""),width=15,height=7.5)
  print(p)
  dev.off()
  p <- ggplotly(p)
  p
  
  # INSERT SIZE
  insertsize_data=insertsize_data[-1,]
  insertsize_data$insert_size=as.numeric(insertsize_data$insert_size)
  insertsize_data$allpairs=as.numeric(insertsize_data$allpairs)
  p <- ggplot(insertsize_data, aes(x=insert_size, y=allpairs,group=SAMPLES,color=Color)) +
    geom_line(color=insertsize_data$Color)+
    theme_classic()+
    facet_grid(SPECIES ~ LOCATION)
  pdf(paste(wd,"/figures/insert_size.pdf",sep=""),width=15,height=7.5)
  print(p)
  dev.off()
  p <- ggplotly(p)
  p
  
  # COVERAGE PLOT
  coverage_data=coverage_data[-1,]
  coverage_data$Coverage=as.numeric(coverage_data$Coverage)
  coverage_data$Count=as.numeric(coverage_data$Count)
  p <- ggplot(coverage_data, aes(x=Coverage, y=Count,group=SAMPLES,color=Color))
  
  x=levels(as.factor(unique(coverage_data$SPECIES)))[1]
  img<-readPNG(paste("Data/",x,".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom2(g, xmin=75, xmax=100, ymin=4e+7, ymax=5e+7, data=coverage_data[1,])
  x=levels(as.factor(unique(coverage_data$SPECIES)))[2]
  img<-readPNG(paste("Data/",x,".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)  
  
  p<-p+  
    annotation_custom2(g, xmin=75, xmax=100, ymin=4e+7, ymax=5e+7, data=coverage_data[(40040/2)+1,])
  
  coverage_data$SPECIES_plot=factor(coverage_data$SPECIES,labels=c("italic('Dicentrarchus labrax')",
                                                               "italic('Sardina pilchardus')"))
  
  
  p <- p +
    geom_line(color=coverage_data$Color)+
    theme_classic()+
    facet_grid(SPECIES_plot ~ .,labeller =label_parsed )+
    scale_x_continuous(
      limits=c(0,100),
      breaks = seq(0,100,by=10),
      label = seq(0,100,by=10)
    )
  pdf(paste(wd,"/figures/coverage.pdf",sep=""),width=15,height=7.5)
  print(p)
  dev.off()
  p <- ggplotly(p)
  p
  
  # GC coverage data
  gc_coverage_data=gc_coverage_data[-1,]
  gc_coverage_data$GC=as.numeric(gc_coverage_data$GC)
  gc_coverage_data$Count=as.numeric(gc_coverage_data$Count)
  gc_coverage_data$SPECIES=factor(gc_coverage_data$SPECIES)
  p <- ggplot(gc_coverage_data, aes(x=GC, y=Count,group=SAMPLES,color=Color))
  
  x=levels(as.factor(unique(gc_coverage_data$SPECIES)))[1]
  img<-readPNG(paste("Data/",x,".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+  
    annotation_custom2(g, xmin=75, xmax=100, ymin=4e+7, ymax=5e+7, data=gc_coverage_data[1,])
  x=levels(as.factor(unique(gc_coverage_data$SPECIES)))[2]
  img<-readPNG(paste("Data/",x,".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)  
  
  p<-p+  
    annotation_custom2(g, xmin=75, xmax=100, ymin=4e+7, ymax=5e+7, data=gc_coverage_data[(40040/2)+1,])
  
  gc_coverage_data$SPECIES_plot=factor(gc_coverage_data$SPECIES,labels=c("italic('Dicentrarchus labrax')",
                                                                   "italic('Sardina pilchardus')"))
  
  
  p <- p +
    geom_line(color=gc_coverage_data$Color)+
    theme_classic()+
    facet_grid(SPECIES_plot ~ .,labeller =label_parsed )+
  pdf(paste(wd,"/figures/gc_coverage.pdf",sep=""),width=15,height=7.5)
  print(p)
  dev.off()
  p <- ggplotly(p)
  p
  
  # Insertion data
  insertion_data=insertion_data[-1,]
  insertion_data$insertion=as.numeric(insertion_data$insertion)
  insertion_data$count=as.numeric(insertion_data$count)
  p <- ggplot(insertion_data, aes(x=insertion, y=count,group=SAMPLES,color=Color)) +
    geom_line(color=insertion_data$Color)+
    scale_x_continuous()+
    theme_bw()+
    facet_grid(SPECIES ~ .)
  
  p <- ggplotly(p)
  p
  
  # Deletion data
  deletion_data=deletion_data[-1,]
  deletion_data$deletion=as.numeric(deletion_data$deletion)
  deletion_data$count=as.numeric(deletion_data$count)
  p <- ggplot(deletion_data, aes(x=deletion, y=count,group=SAMPLES,color=Color)) +
    geom_line(color=deletion_data$Color)+
    scale_x_continuous()+
    theme_bw()+
    facet_grid(SPECIES ~ .)
  
  p <- ggplotly(p)
  p

 
}

