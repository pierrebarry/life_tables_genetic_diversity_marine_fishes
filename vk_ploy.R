### All plots ----
corr_plot<-ggarrange(div_agene,div_slim,labels=c("A","B"),common.legend = T, legend='top')                                    
corr_plot<-annotate_figure(corr_plot,
                           left = text_grob("Observed heterozygosity (%)",rot = 90,vjust=2.5),
)
pairwise_plot_agene<-ggarrange(pairwise_agene,
                               pairwise_agene_nopc,
                               labels=c("C","D"),
                               nrow=1)
pairwise_plot_agene<-annotate_figure(pairwise_plot_agene,
                                     left = text_grob("Ratio of observed Ne/N",rot = 90,vjust=2),
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
