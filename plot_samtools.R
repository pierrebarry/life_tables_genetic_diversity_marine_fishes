
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

  library(knitr)
  library(tidyverse)
  library(reshape2)
  library(gridExtra)
  library(plotly)
  library(ggthemes)
  library(RColorBrewer)
  library(ggiraph)
  library(png)
  library(grid)
  library(ggExtra)
  library(magick)
  library(rsvg)
  
  color_med_atl=data.frame(Location=c("Li","Mu","Fa","Ga"),
                           Col=brewer.pal(n = 4, name = "RdBu"))
  
  if (length(SAMPLES)==0 | length(SPECIES)==0){
    print("No samples or species provided")
  } else {
    
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
      
      #print(SAMPLES[i])
      
      #if (file.exists(paste(wd))){
      
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
      for (j in 1:length(id$length)){
        deletion_data=rbind(deletion_data,
                            c(SAMPLES[i],
                              substr(SAMPLES[i],0,5),
                              substr(SAMPLES[i],6,7),
                              id$length[j],
                              id$deletion_count[j],
                              as.character(color_med_atl[which(substr(SAMPLES[i],6,7)==color_med_atl$Location),2])))
      }
      
      
      #}
      
    }
    
    #GLOBAL STATS
    globalstats=globalstats[-1,]
    globalstats$value=as.numeric(globalstats$value)
    globalstats[(globalstats$stats=="reads mapped:" | 
                  globalstats$stats=="reads unmapped:") & globalstats$SAMPLES==SAMPLES[i],]
    
    for (i in 1:length(SAMPLES)){
      inter=globalstats[(globalstats$stats=="reads mapped:" | 
                     globalstats$stats=="reads unmapped:") & globalstats$SAMPLES==SAMPLES[i],]
      print(c(SAMPLES[i],
              inter$value[1],
              inter$value[1]/(inter$value[1]+inter$value[2])))
    }
    
    # Reads mapped and unmapped
    p<-ggplot(data=globalstats[globalstats$stats=="reads mapped:" | 
                                 globalstats$stats=="reads unmapped:",], aes(x=SAMPLES,y=value,fill=stats)) +
      geom_bar(stat="identity",alpha=0.8)+
      ggtitle("Reads mapped - unmapped")+
      xlab("Samples")+
      ylab("Reads")+
      facet_grid(SPECIES~.,scales='free')+
      coord_flip()+
      geom_vline(xintercept=seq(5.5,15.5,by=5),linetype="dashed",col="chartreuse4")
    pdf(paste(wd,"/figures/samtools.pdf",sep=""),width=7.5,height=7.5)
    print(p)
    dev.off()
    p <- ggplotly(p)
    p
    

    
    # Reads mapped
    p<-ggplot(data=globalstats[globalstats$stats=="reads mapped:",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="reads mapped:",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Reads mapped")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value)*0.01,
             (max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)))+
      ggtitle("Reads mapped")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # Reads mapped and pairds
    p<-ggplot(data=globalstats[globalstats$stats=="reads mapped and paired:",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="reads mapped and paired:",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Reads mapped and pairs")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="reads mapped and paired:",]$value)-min(globalstats[globalstats$stats=="reads mapped and paired:",]$value)*0.01,
             (max(globalstats[globalstats$stats=="reads mapped and paired:",]$value)-min(globalstats[globalstats$stats=="reads mapped and paired:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped and paired:",]$value)))+
      ggtitle("Reads mapped and pairs")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # Reads unmapped
    p<-ggplot(data=globalstats[globalstats$stats=="reads unmapped:",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="reads unmapped:",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Reads unmapped")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="reads unmapped:",]$value)-min(globalstats[globalstats$stats=="reads unmapped:",]$value)*0.01,
             (max(globalstats[globalstats$stats=="reads unmapped:",]$value)-min(globalstats[globalstats$stats=="reads unmapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads unmapped:",]$value)))+
      ggtitle("Reads unmapped")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # Reads properly paired
    p<-ggplot(data=globalstats[globalstats$stats=="reads properly paired:",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="reads properly paired:",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Reads properly paired")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="reads properly paired:",]$value)-min(globalstats[globalstats$stats=="reads properly paired:",]$value)*0.01,
             (max(globalstats[globalstats$stats=="reads properly paired:",]$value)-min(globalstats[globalstats$stats=="reads properly paired:",]$value))*0.1+max(globalstats[globalstats$stats=="reads properly paired:",]$value)))+
      ggtitle("Reads properly paired")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # Average length
    p<-ggplot(data=globalstats[globalstats$stats=="average length:",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="average length:",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Average length")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="average length:",]$value)-min(globalstats[globalstats$stats=="average length:",]$value)*0.01,
             (max(globalstats[globalstats$stats=="average length:",]$value)-min(globalstats[globalstats$stats=="average length:",]$value))*0.1+max(globalstats[globalstats$stats=="average length:",]$value)))+
      ggtitle("Average length")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # Average quality
    p<-ggplot(data=globalstats[globalstats$stats=="average quality:",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="average quality:",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Average length")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="average quality:",]$value)-min(globalstats[globalstats$stats=="average quality:",]$value)*0.01,
             (max(globalstats[globalstats$stats=="average quality:",]$value)-min(globalstats[globalstats$stats=="average quality:",]$value))*0.1+max(globalstats[globalstats$stats=="average quality:",]$value)))+
      ggtitle("Average length")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # Insert size average
    p<-ggplot(data=globalstats[globalstats$stats=="insert size average:",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="insert size average:",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Insert size average")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="insert size average:",]$value)-min(globalstats[globalstats$stats=="insert size average:",]$value)*0.01,
             (max(globalstats[globalstats$stats=="insert size average:",]$value)-min(globalstats[globalstats$stats=="insert size average:",]$value))*0.1+max(globalstats[globalstats$stats=="insert size average:",]$value)))+
      ggtitle("Insert size average")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # Pairs on different chromosomes
    p<-ggplot(data=globalstats[globalstats$stats=="pairs on different chromosomes:",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="pairs on different chromosomes:",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Pairs on different chromosome")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value)-min(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value)*0.01,
             (max(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value)-min(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value))*0.1+max(globalstats[globalstats$stats=="pairs on different chromosomes:",]$value)))+
      ggtitle("Pairs on different chromosome")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # Percentage of properly paired reads
    p<-ggplot(data=globalstats[globalstats$stats=="percentage of properly paired reads (%):",], aes(x=SPECIES, y=value)) +
      geom_violin(color="black",alpha=0.5) +
      geom_boxplot(outlier.shape = NA,alpha=0.75)+
      geom_jitter(position=position_jitter(0.2),alpha=0.75,size=1.5,color=globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$Color,aes(text=SAMPLES))+
      theme_bw()+
      ylab("Percentage of properly paired reads")+
      xlab("Species")+
      ylim(c(min(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value)-min(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value)*0.01,
             (max(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value)-min(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value))*0.1+max(globalstats[globalstats$stats=="percentage of properly paired reads (%):",]$value)))+
      ggtitle("Percentage of properly paired reads")+
      scale_x_discrete(labels=c(#"Coryphoblennius \n galerita",
        #"Coris \n julis",
        "Dicentrarchus \n labrax",
        #"Diplodus \n puntazzo",
        #"Engraulis \n encrasicolus",
        #"Hippocampus \n guttulatus",
        #"Lophius \n budegassa",
        #"Lithognathus \n mormyrus",
        #"Merluccius \n merluccius",
        #"Mullus \n surmuletus",
        #"Pagellus \n erythrinus",
        #"Serranus \n cabrilla",
        #"Spondyliosoma \n cantharus",
        #"Symphodus \n cinereus",
        "Sardina \n pilchardus"
        #"Sarda \n sarda",
        #"Syngnathus \n typhle"
      ))
    
    #for (j in 1:length(levels(as.factor(unique(globalstats$SPECIES))))){
    #   
    #  x=levels(as.factor(unique(globalstats$SPECIES)))[j]
    #   img<-readPNG(paste("C:/Users/pierr/Documents/PROJETS/COGEDIV/OUTPUT/PHOTOS/",x,".png",sep=""))
    #   g <- rasterGrob(img, interpolate=TRUE)
    
    #   p<-p+  
    #    annotation_custom(g, xmin=j-0.4, xmax=j+0.4, 
    #                     ymin=max(globalstats[globalstats$stats=="reads mapped:",]$value)+0.01*max(globalstats[globalstats$stats=="reads mapped:",]$value),
    #                     ymax=(max(globalstats[globalstats$stats=="reads mapped:",]$value)-min(globalstats[globalstats$stats=="reads mapped:",]$value))*0.1+max(globalstats[globalstats$stats=="reads mapped:",]$value)) 
    
    # }
    
    p <- ggplotly(p)
    p
    
    # READ LENGTH
    readlength_data=readlength_data[-1,]
    readlength_data$read_length=as.numeric(readlength_data$read_length)
    readlength_data$Count=as.numeric(readlength_data$Count)
    p <- ggplot(readlength_data, aes(x=read_length, y=Count,group=SAMPLES)) +
      geom_line(color=readlength_data$Color)+
      theme_bw()+
      facet_grid(SPECIES ~ LOCATION)
    p <- ggplotly(p)
    p
    
    
    # INSERT SIZE
    insertsize_data=insertsize_data[-1,]
    insertsize_data$insert_size=as.numeric(insertsize_data$insert_size)
    insertsize_data$allpairs=as.numeric(insertsize_data$allpairs)
    p <- ggplot(insertsize_data, aes(x=insert_size, y=allpairs,group=SAMPLES,color=Color)) +
      geom_line(color=insertsize_data$Color)+
      theme_bw()+
      facet_grid(SPECIES ~ LOCATION)
    
    p <- ggplotly(p)
    p
    
    # COVERAGE PLOT
    coverage_data=coverage_data[-1,]
    coverage_data$Coverage=as.numeric(coverage_data$Coverage)
    coverage_data$Count=as.numeric(coverage_data$Count)
    p <- ggplot(coverage_data, aes(x=Coverage, y=Count,group=SAMPLES,color=Color)) +
      geom_line(color=coverage_data$Color)+
      scale_x_continuous()+
      theme_bw()+
      facet_grid(SPECIES ~ .)+
      xlim(c(0,150))
    
    p <- ggplotly(p)
    p
    
    # GC coverage data
    gc_coverage_data=gc_coverage_data[-1,]
    gc_coverage_data$GC=as.numeric(gc_coverage_data$GC)
    gc_coverage_data$Count=as.numeric(gc_coverage_data$Count)
    p <- ggplot(gc_coverage_data, aes(x=GC, y=Count,group=SAMPLES,color=Color)) +
      geom_line(color=gc_coverage_data$Color)+
      scale_x_continuous()+
      theme_bw()+
      facet_grid(SPECIES ~ .)
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
    
    #ACGT
    for (sp in 1:length(unique(SPECIES))){
      par(mfrow=c(4,5))
      for (sam in 1:length(SAMPLES[SPECIES==unique(SPECIES)[sp]])){
        stats<- readLines(paste(wd,"/",unique(SPECIES)[sp],"/Mapping/",SAMPLES[SPECIES==unique(SPECIES)[sp]][sam],"/",SAMPLES[SPECIES==unique(SPECIES)[sp]][sam],"_samtools_stats.stats",sep=""))
        actg <- grep("^GCC",stats, value=TRUE)
        actg <- separate(data.frame(actg),col=1, into=c("ID", "cycle", "A", "C", "G", "T", "N", "O"), sep="\t")[,-1]
        plot(actg$cycle,actg$A,type="l",col="chartreuse3",ylim=c(0,100),xlab='Cycle',ylab="GC(%)",main=SAMPLES[SPECIES==unique(SPECIES)[sp]][sam])
        lines(actg$cycle,actg$T,type="l",col="dodgerblue2")
        lines(actg$cycle,actg$C,type="l",col="black")
        lines(actg$cycle,actg$G,type="l",col="red")
      }
    }
    
    #GC content
    for (sp in 1:length(unique(SPECIES))){
      par(mfrow=c(4,5))
      for (sam in 1:length(SAMPLES[SPECIES==unique(SPECIES)[sp]])){
        stats<- readLines(paste(wd,"/",unique(SPECIES)[sp],"/Mapping/",SAMPLES[SPECIES==unique(SPECIES)[sp]][sam],"/",SAMPLES[SPECIES==unique(SPECIES)[sp]][sam],"_samtools_stats.stats",sep=""))
        gc <- grep("^GCF|^GCL",stats, value=TRUE)
        gc <- separate(data.frame(gc),col=1, into=c("Pair", "GC", "Count"), sep="\t")
        plot(gc[gc$Pair=="GCF",]$GC,gc[gc$Pair=="GCF",]$Count,type="l",col="blue",xlab="GC(%)",ylab="Count",main=SAMPLES[SPECIES==unique(SPECIES)[sp]][sam])
        lines(gc[gc$Pair=="GCL",]$GC,gc[gc$Pair=="GCL",]$Count,type="l",col="red")
      }
    }
    
    
    
    
  }






sn <- grep("^SN",stats, value=TRUE)
sn <- separate(data.frame(sn),col=1, into=c("ID", "Name","Value"), sep="\t")[,-1]
kable(sn, caption="Summary numbers")

# Quality
fq <- grep("^FFQ|^LFQ",data, value=TRUE)
fq <- separate(data.frame(fq),col=1, into=c("Pair", "Cycle", seq(43)), sep="\t")
plot(fq$Cycle,fq$'38')


