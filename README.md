# MutationTimeR
MutationTimeR is an R package to time somatic mutations relative to clonal and subclonal copy number states and calculate
the relative timing of copy number gains.
It has been used by the PCAWG consortium to calculate mutation times of 2,778 whole genome sequencing samples. http://dx.doi.org/10.1101/161562	

## Installation
MutationTimeR runs in most current `R` versions. You can install and load MutationTimeR using
```{R}
devtools::instal_github("gerstung-lab/MutationTimeR")
```

## Running MutationTimeR
MutationTimeR requires a `vcf` object containing point mutations (SNVs, MNVs, or indels) and a `GRanges` object with the copy number segments of the sample.

```{R}
library("MutationTimeR")
vcf <- readVcf("myvcf.vcf") # Point mutations, needs `info` columns t_alt_count t_ref_count
bb <- GRanges(, major_cn= , minor_cn=, clonal_frequency=purity) # Copy number segments, needs columns  major_cn, minor_cn and clonal_frequency of each segment
```

To run MutationTimeR simply use
```{R}
mt <- mutationTime(vcf, bb)
```

## Annotation of point mutations
MutationTimer produces two main outputs. The first is the probability for each individual point mutation from the original `vcf` object to belong to different copy number levels:
```{R}
head(param$D)
```

These can be converted into basic clonal states (`early clonal/late clonal/clonal/subclonal/NA`)
```{R}
table(classifyMutations(mt$D))
```

## Timing of copy number gains
The proportion of point mutations in different copy number states, in turn, defines the molecular timing of gains. This can be calculated using: 
```{R}
bb$timing_param <- mt$P
bbToTime(bb)
```

<img src="MutationTime.png" alt="MutationTime.R output" width="50%"/>
