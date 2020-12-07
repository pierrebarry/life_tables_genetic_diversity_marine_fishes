#------------------------------------------------------------------------#
#                                                                        #
#     This R script generates figures and data presented                 #
#                         in the  article :                              #
#                                                                        #
#     The impact of life-history traits on genetic diversity             #
#                       in marine fishes                                 #
#                                                                        #
#     Pierre Barry, Thomas Broquet, Pierre-Alexandre Gagnaire            #
#                                                                        #
#  Each lines executes a specific R script present in .                  #
#                                                                        #
#------------------------------------------------------------------------#
Salut !
# Choose your working directory for LaTeX images and R scripts -----
wd="C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT"
setwd(wd)
code_file="CoGeDiv_GeneticDiversity/"
library(rmarkdown)
## Map
source(file=paste(code_file,"sampling_map.R",sep=""))
# Fastp output ----
source(file=paste(code_file,"fastp.R",sep=""))
# GenomeScope : estimation, sensibility test, comparison with variant calling-----
source(file=paste(code_file,"GenomeScope.R",sep=""))
# Life-history traits and genetic diversity ----
#source(file=paste(code_file,"lfh_diversity.R",sep=""))
# AgeNe
#source(file=paste(code_file,"agene_analysis.R",sep=""))
# Forward simulation
#source(file=paste(code_file,"SLiM_analysis.R",sep=""))
# Simulation lifetimes
#source(file=paste(code_file,"sim_lifetime_analysis.R",sep=""))
# Historical contingencies
#source(file=paste(code_file,"historical_contigencies.R",sep=""))

## samtools
#source(file=paste(code_file,"plot_samtools.R",sep=""))

render("C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT/Notebook/notebook.Rmd","all")

