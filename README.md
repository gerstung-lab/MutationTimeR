# MutationTimeR
MutationTimeR is an R package to time somatic mutations relative to clonal and subclonal copy number states and calculate
the relative timing of copy number gains. Time is measured as a fraction of point mutations; this is termed _mutation time_. Mutation time 
is proportional to real time if the number of mutations acquired per bp and per year is constant. 

MutationTimeR has been used by the PCAWG consortium to calculate mutation times of 2,778 whole genome sequencing samples. 
Please see M. Gerstung, C. Jolly, I. Leshchiner, S. Dentro, S. Gonzalez _et al._, [The Evolutionary History of 2,658 Cancers](https://doi.org/10.1038/s41586-019-1907-7), _Nature_. *578*, pages 122-128(2020).	

## Installation
MutationTimeR runs in most current `R` versions. You can install and load MutationTimeR using

```r
devtools::install_github("mg14/mg14")
devtools::install_github("gerstung-lab/MutationTimeR")
```

## Input data
MutationTimeR requires a `vcf` object containing point mutations (SNVs, MNVs, or indels) and a `GRanges` object with the copy number segments of the sample.


```r
library("MutationTimeR")
#vcf <- readVcf("myvcf.vcf") # Point mutations, needs `geno` entries `AD` and `DP` or `info` columns t_alt_count t_ref_count.
#bb <- GRanges(, major_cn= , minor_cn=, clonal_frequency=purity) # Copy number segments, needs columns  major_cn, minor_cn and clonal_frequency of each segment
#clusters <- data.frame(cluster= , n_ssms=, proportion=) # Optional data.frame with subclonal cluster locations (VAF proportion) and size (number of variants n_ssms)
```

Here we cut this short by

```r
data(MutationTimeR)
```

The `vcf` object containing all point mutations is:

```r
vcf
```

```
## class: ExpandedVCF 
## dim: 5499 1 
## rowRanges(vcf):
##   GRanges with 4 metadata columns: REF, ALT, QUAL, FILTER
## info(vcf):
##   DataFrame with 2 columns: t_alt_count, t_ref_count
## info(header(vcf)):
##                Number Type    Description     
##    t_alt_count 1      Integer Tumour alt count
##    t_ref_count 1      Integer Tumour ref count
## geno(vcf):
##   SimpleList of length 3: AD, DP, FT
## geno(header(vcf)):
##       Number Type    Description                                             
##    AD 2      Integer Allelic depths (number of reads in each observed allele)
##    DP 1      Integer Total read depth                                        
##    FT 1      String  Variant filters
```

The `GRanges` with allele-specific copy number data

```r
bb
```

```
## GRanges object with 80 ranges and 3 metadata columns:
##        seqnames              ranges strand |  major_cn  minor_cn
##           <Rle>           <IRanges>  <Rle> | <integer> <integer>
##    [1]        1     17500-144915331      * |         1         1
##    [2]        1 144915332-147959385      * |         2         1
##    [3]        1 147959386-149718742      * |         2         1
##    [4]        1 149718743-249250620      * |         2         1
##    [5]        2     17500-242979500      * |         1         1
##    ...      ...                 ...    ... .       ...       ...
##   [76]       20   25006082-25168069      * |         1         1
##   [77]       20   25168070-29532764      * |         1         0
##   [78]       20   29532765-63025519      * |         1         1
##   [79]       21          1-48129894      * |         1         1
##   [80]       22   12200000-51304565      * |         1         0
##        clonal_frequency
##               <numeric>
##    [1]             0.54
##    [2]             0.54
##    [3]             0.54
##    [4]             0.54
##    [5]             0.54
##    ...              ...
##   [76]             0.54
##   [77]             0.54
##   [78]             0.54
##   [79]             0.54
##   [80]             0.54
##   -------
##   seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

Lastly, the prevalence and number of subclonal is

```r
clusters
```

```
##   cluster n_ssms proportion
## 1       0   5518       0.54
## 2       1    845       0.23
```

## Running MutationTimeR
To run MutationTimeR simply use


```r
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
```

## Annotation of point mutations
MutationTimer produces two main outputs. The first is the probability for each individual point mutation from the original `vcf` object to belong to different copy number levels:


```r
head(mt$V)
```

```
## DataFrame with 6 rows and 15 columns
##       MutCN MutDeltaCN     MajCN     MinCN  MajDerCN  MinDerCN       CNF
##   <numeric>  <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
## 1         1          0         1         1         1         1      0.54
## 2         1          0         1         1         1         1      0.54
## 3         1          0         1         1         1         1      0.54
## 4         1          0         1         1         1         1      0.54
## 5         1          0         1         1         1         1      0.54
## 6         1          0         1         1         1         1      0.23
##            CNID            pMutCN     pGain              pSingle
##   <IntegerList>         <numeric> <numeric>            <numeric>
## 1             1 0.999999992268111         0    0.999999992268111
## 2             1   0.9999753674765         0      0.9999753674765
## 3             1 0.999937495120268         0    0.999937495120268
## 4             1 0.999999348396561         0    0.999999348396561
## 5             1 0.999999960786733         0    0.999999960786733
## 6             1 0.999795151536053         0 0.000204848463947305
##                   pSub pAllSubclones         pMutCNTail         CLS
##              <numeric>        <List>          <numeric>    <factor>
## 1 7.73188863450679e-09                0.976558542742796 clonal [NA]
## 2  2.4632523499745e-05                0.653374217923694 clonal [NA]
## 3 6.25048797320028e-05                0.336557312762548 clonal [NA]
## 4 6.51603439185877e-07                0.926037697462508 clonal [NA]
## 5 3.92132672834029e-08                0.940677539480767 clonal [NA]
## 6    0.999795151536053               0.0814516989630261   subclonal
```

These probabilities are the basis for the following simple clonal states (`early clonal/late clonal/clonal/subclonal/NA`)


```r
table(mt$V$CLS)
```

```
## 
## clonal [early]  clonal [late]    clonal [NA]      subclonal 
##            198             97           4356            848
```

To add this annotation to the vcf use

```r
vcf <- addMutTime(vcf, mt$V)
```

## Timings of copy number gains
The proportion of point mutations in different copy number states, in turn, defines the molecular timing of gains. These are stored in the following `DataFrame`:


```r
head(mt$T)
```

```
## DataFrame with 6 rows and 10 columns
##                                                                        timing_param
##                                                                              <List>
## 1                                                    1:1:0.27:...,2:1:0.115:...,...
## 2 1:1:0.21259842519685:...,1:2:0.425196850393701:...,2:1:0.0905511811023622:...,...
## 3 1:1:0.21259842519685:...,1:2:0.425196850393701:...,2:1:0.0905511811023622:...,...
## 4 1:1:0.21259842519685:...,1:2:0.425196850393701:...,2:1:0.0905511811023622:...,...
## 5                                                    1:1:0.27:...,2:1:0.115:...,...
## 6                                                    1:1:0.27:...,2:1:0.115:...,...
##                type              time           time.lo           time.up
##            <factor>         <numeric>         <numeric>         <numeric>
## 1                NA                NA                NA                NA
## 2 Mono-allelic Gain                 1 0.485179575172229                 1
## 3 Mono-allelic Gain                NA                NA                NA
## 4 Mono-allelic Gain 0.773534222763912 0.666501932404523 0.815455093388117
## 5                NA                NA                NA                NA
## 6                NA                NA                NA                NA
##    time.2nd time.2nd.lo time.2nd.up time.star n.snv_mnv
##   <numeric>   <numeric>   <numeric>  <factor> <integer>
## 1        NA          NA          NA        NA       148
## 2        NA          NA          NA       ***         8
## 3        NA          NA          NA        NA         0
## 4        NA          NA          NA       ***       128
## 5        NA          NA          NA        NA       304
## 6        NA          NA          NA        NA        22
```
The relevant columns are `time` with 95% confidence intervals `time.lo` and `time.hi` as well as the counterparts `time2/time2.lo/time2.hi` for the second gain in cases where one allele has more than one gained copy.
The field `time.star` indicates a tiers: `***` indicates gains +1. `**` gains +2, which are found to be slightly less reliable and need certain assumptions about their temporal sequence. `*` are subclonal gains which are hit an miss.

The DataFrame can be added to the copy number `GRanges` object for convenience.


```r
mcols(bb) <- cbind(mcols(bb),mt$T)
```

## Plot output
Timing annotated VCF and copy number can be plotted using the following command:


```r
plotSample(vcf,bb)
```

<img src="MutationTimeR.png" width="576" />

This shows the observed and expected variant allele frequencies of point mutations on the top. This is very useful to spot inconsistencies with purity and copy number configuration. As a rule of thumb the states (horizontal bars) should run
right through the middle of the clouds of point mutations. Colours indicate the timing category: Blue = clonal [other], purple = clonal [late], green = clonal [early], red = subclonal. 

The middle plot shows the copy number as stacked barplots. Subclonal CN is indicated by fractional bars. Dark grey is major, light grey minor allele.

The bottom plot shows the estimated mutation time of primary and secondary gains (shaded). Boxes denote 95% CIs. The histogram at the right shows the distribution of timing events. Blue = mono-allelic gains (N:1), pink = CN-LOH/gain+loss (N:0) and green = bi-allelic gains (N:2).
