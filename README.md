# :fish: Life tables drive genetic diversity in marine fishes :dna:

Scripts and files to generate results and output of the article 

> Barry P., Broquet T., Gagnaire P.-A., 2021. Life tables drive genetic diversity in marine fishes. Evolution Letters

## Sampling information

:file_folder: Files:
- `sampling.Rdata`: sampling information of sequenced individuals readable in R
- `Data/GENETIC_DIVERSITY_DATA.xlsx` (sheet : SAMPLING) : sampling information of sequenced invidiausl readable in Excel format.
- Various `Data/*.png` files representing the 16 species of the study from Iglesias, S. (2013). Actinopterygians from the North-Eastern Atlantic and the Mediterranean (A 750
Natural Classification Based on Collection Specimens, with DNA Barcodes and Standardized 751
Photographs), Volume I (Plates), Provisional Version 09. (reproduced with permission).

:bar_chart: Scripts:
- `Data/sampling_map.R` : display sampling locations of all individuals

## Pre-processing fasta files

:file_folder: Files :
- `Data/fastp_process.Rdata`: sequence data after quality correction

:bar_chart: Scripts :
- `Snakefile_fastp` : snakefile running fastp command below :
```ruby
fastp -i {INPUT.R1} -I {INPUT.R2} -o {OUTPUT.FASTP_R1} -O {OUTPUT.FASTP_R2} --trim_poly_g --correction --low_complexity_filter --html {OUTPUT.REPORT_HTML} --json {OUTPUT.REPORTJSON} --report_title {SAMPLE} --thread 8 --dont_overwrite
```
- `Extract_fastp_info.py` : extract individual quality statistics from json output by fastp.
- `fastp.R` : display quality statistics of individual paired-end raw sequences

## Estimate individual genome-wide genetic diversity

:file_folder: Files :
- `Data/Summary_GenomeScope.txt` : output of individual genome-wide statistics estimates with GenomeScope.
- `Data/Sensibility_test.Rdata` : test of choice of k in jellyfish on the estimation of individual genetic diversity readable in R.
- `Data/diversity_Dlabr.het` : estimate of individual genome-wide genetic diversity with vcftools after mapping to reference genome and variant calling with GATK.
- `Data/het_Spilc.het` : estimate of individual genome-wide genetic diversity with vcftools after mapping to reference genome and variant calling with GATK.
- `samtools/*_flagstat.txt` and `samtools/*_samtools_stats.stats` : individual statistics of mapping statistics.

:bar_chart: Scripts :
- `Snakefile_genomescope` : snakefile running estimation of individual genetic diversity with command below :
```ruby
jellyfish count -C -m 21 -s 1000000 -t 8 {INPUT.FASTQ} -o {OUTPUT.JELLYFISH}
jellyfish histo -t 10 {INPUT.JELLYFISH} > {OUTPUT.HISTO}
R --vanilla --slave --args {OUTPUT.HISTO} 21 150 OutputGenomeScope 1000000 Summary {SAMPLE} < ~/GenomeScope/genomescope_cluster.R
```
- `genomescope_cluster.R` : `GenomeScope` algorithm to estimate individual genetic diversity (adapted from original scripts available [here](https://github.com/schatzlab/genomescope) and published in Vurture, G. W., Sedlazeck, F. J., Nattestad, M., Underwood, C. J., Fang, H., Gurtowski, J., and Schatz, M. C. (2017). GenomeScope : Fast reference-free genome profiling from shortreads. Bioinformatics, 33(14) :2202–2204) 
- `genomescope.R` : display individual estimated and standard deviation of genome-wide genetic diversity, genome length, genome unique length, genome repeat length and model fit estimated by  `GenomeScope`
- `historical_contingencies.R` : heatmap clustering of intraspecific variance in genetic diversity.
- `Snakefile_samtools_GATK` : snakefile running mapping, marking duplicates and variant calling.
- `plot_samtools` : display mapping statistics after individual mapping to reference genome.
- `vcftools_het` : estimate individual genetic diversity with vcftools from VCF generated after mapping and variant calling for D.labrax and S.pilchardus.

## Simple determinants of genetic diversity

:file_folder: Files :
- `Data/GENETIC_DIVERSITY_DATA.xlsx` (sheet : lfh) : life-history traits values for the  species (see supplementary files for references)
- `Data/div.Rdata` : species's median genetic diversity. 
- `Data/genome_length.Rdata` : species's median genome length.

:bar_chart: Scripts :
- `lfh_diversity.R` : correlates nine life history traits (body size, trophic level, fecundity, propagule size, age at first maturity, lifespan, adult lifespan, hermaphroditism, brooding strategy) with median species genetic diversity.

## Estimation Ne/N and correlation with genetic diversity 

:file_folder: Files :
- `GENETIC_DIVERSITY_DATA.xlsx` (sheet : SLiM) : life tables characteristics of the  species

:bar_chart: Scripts :
- `plot_age_surv_fec.R` : show  sex-specific age-specific survival and fecundity, cumulative survival for the 16 species retrieving from the litterature.

### AgeNe

:file_folder: Files : 
- `agene/AgeNe.exe` : AgeNe algorithm from Waples, R. S., Do, C., and Chopelet, J. (2011). Calculating Ne and Ne/N in age-structured populations : A hybrid Felsenstein-Hill approach. Ecology, 92(7) :1513–1522 (download [here](https://figshare.com/articles/dataset/Supplement_1_AgeNe_a_program_to_calculate_Ne_and_Nb_in_age-structured_populations_/3551643?backTo=/collections/Calculating_i_N_i_sub_e_sub_and_i_N_i_sub_e_sub_i_N_i_in_age-structured_populations_a_hybrid_Felsenstein-Hill_approach/3304059))
- `agene/agene.xlsx` : life tables characteristics (identical to `GENETIC_DIVERSITY_DATA.xlsx` in sheet SLiM)
- `agene/cogediv_sensbility.txt` : example of input of `AgeNe`
- `agene/output_sensibility.txt` : example of output of `AgeNe`
- `agene/agene_output.Rdata` : species estimate of variance in reproductive success for each life tables characteristcis.

:bar_chart: Scripts :
- `Data/agene.R` : creates input and run `AgeNe` from life tables. 
- `agene_analysis.R` : estimate of Ne/N from `AgeNe` simulations and correlates with species genetic diversity.

### SLiM

:file_folder: Files :
- `Input*/` : input directories containing species input files for `SLiM` simulations.

:bar_chart: Scripts :
- `forward_simulation_slim.R` : creates input for `SLiM` simulations
- `forward_slim/snakefile` : snakefile to run `SLiM` simulations with script below :
```ruby
SLiM/build/slim -t -m -d K={params.K} -d L={params.L} -d iter={wildcards.itera} -d mu={params.mu} -d rec={params.rec} -d register_each={params.register_each} {input.input}
```
- `forward_slim/job` : job scripts to run on IFB cluster.
- `Slim_analysis.R` : estimate of Ne/N from `SLiM` simulations and correlates with species genetic diversity.

## Simulated life tables from alternatives age-specific survival and fecundity curves

:file_folder: Files :
- `sim_lifetime/AgeNe.exe` : AgeNe algorithm from Waples, R. S., Do, C., and Chopelet, J. (2011). Calculating Ne and Ne/N in age-structured populations : A hybrid Felsenstein-Hill approach. Ecology, 92(7) :1513–1522 (download [here](https://figshare.com/articles/dataset/Supplement_1_AgeNe_a_program_to_calculate_Ne_and_Nb_in_age-structured_populations_/3551643?backTo=/collections/Calculating_i_N_i_sub_e_sub_and_i_N_i_sub_e_sub_i_N_i_in_age-structured_populations_a_hybrid_Felsenstein-Hill_approach/3304059))
- `sim_lifetime/cogediv_sensbility.txt` : example of input of `AgeNe`
- `sim_lifetime/output_sensibility.txt` : example of output of `AgeNe`
- `simulation_lifetime_linear_0.01.Rdata` : Ne/N of the 16 simulated species for all the combinations of c and f values from linear age-fecundity model. 
- `simulation_lifetime_exp_0.01_noscale.Rdata` : Ne/N of the 16 simulated species for all the combinations of c and f values from exponential age-fecundity model. 
- `simulation_lifetime_power_0.01_noscale.Rdata` : Ne/N of the 16 simulated species for all the combinations of c and f values from power-law age-fecundity model. 
- `simulation_lifetime_poly_0.01_first.Rdata` : Ne/N of the 16 simulated species for all the combinations of c and f values from polynomial age-fecundity model. 

:bar_chart: Scripts :
- `sim_lifetime.R` : generates simulated life tables and estimate Ne/N with `AgeNe`
- `sim_lifetime_analysis.R` : estimate slope between adult lifespan and Ne/N estimations

## Figures

:bar_chart: Scripts :
- `main_figure.R` : export Figure 1 and Figure 3 from the main text.

## :arrow_right: Running all scripts

`Master_Script.R` run all scripts and render final notebook.

## :wrench: Tools needed

* [fastp v0.20.0](https://github.com/OpenGene/fastp)
* [jellyfish](https://github.com/gmarcais/Jellyfish)
* [genomescope](https://github.com/schatzlab/genomescope)
* [SLiM](https://messerlab.org/slim/)
* [AgeNE](https://figshare.com/articles/dataset/Supplement_1_AgeNe_a_program_to_calculate_Ne_and_Nb_in_age-structured_populations_/3551643?backTo=/collections/Calculating_i_N_i_sub_e_sub_and_i_N_i_sub_e_sub_i_N_i_in_age-structured_populations_a_hybrid_Felsenstein-Hill_approach/3304059)
* [snakemake](https://github.com/snakemake/snakemake)

## :red_circle: Youtube video

For french-speakers, [a youtube video](https://www.youtube.com/watch?v=98pTKuRNgAE&list=PL_rJBQvKDsY--gDXJ5d21QWQboa9OoGmh&index=3) summarising in 10 minutes the article.