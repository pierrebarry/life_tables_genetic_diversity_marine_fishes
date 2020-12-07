## Historical contingencies -----
library(gplots)
wd="C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT"
setwd(wd)
load(file="Data/Summary_GenomeScope.Rdata")
ss=Summary_GenomeScope[Summary_GenomeScope$Sample!="DpuntMu5",]
scale_centered_species=c()
for(i in 1:nrow(ss)){
  scale_centered_species[i]=(ss$Heterozygosity[i]-mean(ss[ss$Species==ss$Species[i],]$Heterozygosity))/(sd(ss[ss$Species==ss$Species[i],]$Heterozygosity))
}

#pp_sp

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
pdf(paste(wd,"/figures/historical.pdf",sep=""),width=5,height=5)
heatmap(div_loc,
        col=RColorBrewer::brewer.pal(5000, "Spectral"),
        #col=plasma(500),
        labRow = as.expression(lapply(rownames(div_loc), function(a) bquote(italic(.(a))))),
        margins=c(5,7),
        cexCol=1.5
)
dev.off()


