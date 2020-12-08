# Life tables drive genetic diversity in marine fishes

Scripts and files to generate results and output of the article 

> Barry P., Broquet T., Gagnaire P.-A., 2021. Life tables drive genetic diversity in marine fishes. Evolution Letters

## Sampling information

`sampling_map.R` : display sampling locations of all individuals

`plot_age_surv_fec.R` : show  sex-specific age-specific survival and fecundity, cumulative survival for the 16 species retrieving from the litterature.

## Pre-processing fasta files

```ruby
fastp -i {INPUT.R1} -I {INPUT.R2} -o {OUTPUT.FASTP_R1} -O {OUTPUT.FASTP_R2} --trim_poly_g --correction --low_complexity_filter --html {OUTPUT.REPORT_HTML} --json {OUTPUT.REPORTJSON} --report_title {SAMPLE} --thread 8 --dont_overwrite
```

`fastp.R` : display quality statistics of individual paired-end raw sequences

## Estimate individual genome-wide genetic diversity

```ruby
jellyfish count -C -m 21 -s 1000000 -t 8 {INPUT.FASTQ} -o {OUTPUT.JELLYFISH}
jellyfish histo -t 10 {INPUT.JELLYFISH} > {OUTPUT.HISTO}
R --vanilla --slave --args {OUTPUT.HISTO} 21 150 OutputGenomeScope 1000000 Summary {SAMPLE} < ~/GenomeScope/genomescope_cluster.R
```

`genomescope.R` : display individual estimated and standard deviation of genome-wide genetic diversity, genome length, genome unique length, genome repeat length and model fit estimated by `GenomeScope`

`historical_contingencies.R` : heatmap clustering of intraspecific variance in genetic diversity.

## Simple determinants of genetic diversity

`lfh_diversity.R` : correlates nine life history traits (body size, trophic level, fecundity, propagule size, age at first maturity, lifespan, adult lifespan, hermaphroditism, brooding strategy) with median species genetic diversity.

## Estimation $\frac{N_{e}}{N}$ and correlation with genetic diversity 

### AgeNe

`forward_simulation_slim.R` : creates input for `SLiM` simulations

`agene_analysis.R` : estimate of $\frac{N_{e}}{N}$ from `AgeNe` simulations and correlates with species genetic diversity.

### SLiM

`forward_simulation_slim.R` : creates input for `SLiM` simulations

```ruby
SLiM/build/slim -t -m -d K={params.K} -d L={params.L} -d iter={wildcards.itera} -d mu={params.mu} -d rec={params.rec} -d register_each={params.register_each} {input.input}
```

`slim_analysis.R` : estimate of $\frac{N_{e}}{N}$ from `SLiM` simulations and correlates with species genetic diversity.

## Simulated life tables from alternatives age-specific survival and fecundity curves

`sim_lifetime.R` : generates simulated life tables and estimate $\frac{N_{e}}{N}$ with `AgeNe`

`sim_lifetime_analysis.R` : estimate slope between adult lifespan and $\frac{N_{e}}{N}$ estimations

## Figures



## Running all scripts

`Master_Script.R` run all scripts and render 

## Tools needed

* [fastp](https://github.com/OpenGene/fastp)
* [jellyfish](https://github.com/gmarcais/Jellyfish)
* [genomescope](https://github.com/schatzlab/genomescope)
* [SLiM](https://messerlab.org/slim/)
* [AgeNE](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/10-1796.1)
* [snakemake](https://github.com/snakemake/snakemake)

## Youtube video

For french-speakers, [a youtube video] (https://www.youtube.com/watch?v=98pTKuRNgAE&list=PL_rJBQvKDsY--gDXJ5d21QWQboa9OoGmh&index=3) summarising in 10 minutes the article.