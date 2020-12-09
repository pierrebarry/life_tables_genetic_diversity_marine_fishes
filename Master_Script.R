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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Choose your working directory for LaTeX images and R scripts -----
wd="C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT"
setwd(wd)
code_file="CoGeDiv_GeneticDiversity/"
library(rmarkdown)
# Map ----
source(file=paste(code_file,"sampling_map.R",sep=""))
# Fastp output ----
source(file=paste(code_file,"fastp.R",sep=""))
# GenomeScope : estimation, sensibility test, comparison with variant calling-----
source(file=paste(code_file,"genomeScope.R",sep=""))
# Life-history traits and genetic diversity ----
source(file=paste(code_file,"lfh_diversity.R",sep=""))
# Mapping statistics ----
source(file=paste(code_file,"plot_samtools.R",sep=""))
# Historical contingencies ----
source(file=paste(code_file,"historical_contigencies.R",sep=""))
# Run AgeNe ----
source(file=paste(code_file,"agene.R",sep=""))
# Get SLiM input ----
source(file=paste(code_file,"forward_simulation_slim.R",sep=""))
# Analyse AgeNe output ----
source(file=paste(code_file,"agene_analysis.R",sep=""))
# Analyses of forward simulation output ----
source(file=paste(code_file,"SLiM_analysis.R",sep=""))
# Run simulated life tables with AgeNe ----
source(file=paste(code_file,"sim_lifetime.R",sep=""))
# Analyses of simulated life tables output ----
source(file=paste(code_file,"sim_lifetime_analysis.R",sep=""))
# Plot figure of the main text of the paper ----
source(file=paste(code_file,"main_figure.R",sep=""))
# RMARKDOWN ----
render("C:/Users/ordinateur/ownCloud/COGEDIV/ARTICLE/Genetic_diversity_LHT/Notebook/notebook.Rmd","all")

