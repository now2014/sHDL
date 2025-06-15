# Stratified high-definition likelihood inference of heritability enrichment <img src="logo.png" align="right" height=140/>

<br>

## Introduction

The sHDL software provides a powerful tool for analyzing the genetic architecture of complex traits through the distribution of heritability across the genome. By integrating state-of-the-art statistical techniques with genomic functional annotations, sHDL offers a comprehensive approach to understanding the genetic basis of complex traits.

sHDL builds upon the success of the stratified linkage disequilibrium score regression (sLDSC) method by introducing a new stratified high-definition likelihood model. This advanced statistical model enhances the estimation efficiency for heritability enrichment parameters while reducing bias, as demonstrated through simulations and real-data analyses. Compared to sLDSC, sHDL offers 1.4 to 7.4-fold higher efficiency in estimating heritability enrichment parameters.

One of the key features of sHDL is its ability to incorporate diverse genomic annotations, including gene expressions specific to cell types and cancer epigenetic information. This enables the identification of specific genomic regions and cell types associated with complex traits, leading to a deeper understanding of their genetic architecture.

Through extensive validation and comparison with sLDSC, sHDL has been shown to discover novel signals and provide more accurate estimations of heritability enrichment parameters. This makes sHDL a robust and versatile tool for researchers in the field of quantitative genetics, offering new insights into the polygenic contributions to complex traits.

## Installation

To install the latest version of `sHDL` package via Github, run the following code in `R`:

```R
# install.packages('remotes')
remotes::install_github("now2014/sHDL", ref='main')
```

## Quick vignette

Download the LD reference panel from the [HDL](https://github.com/zhenin/HDL) software:

```bash
wget -c -t 1 \
  https://www.dropbox.com/s/fuvpwsf6r8tjd6c/UKB_array_SVD_eigen90_extraction.tar.gz?dl=0 \
  --no-check-certificate -O /path/to/UKB_array_SVD_eigen90_extraction.tar.gz

# checking the md5 = ff3fadd7ea08bd29759b6c652618cd1f
md5sum /your/path/to/UKB_array_SVD_eigen90_extraction.tar.gz

# extraction
tar -xvf /your/path/to/UKB_array_SVD_eigen90_extraction.tar.gz
```

Run sHDL

```R
remotes::install_github("zhenin/HDL/HDL")
data(gwas1.example, package="HDL")

library(sHDL)

## The GWAS summary statistics for birth weight loaded from HDL package.
M <- nrow(gwas1.example)
set.seed(1234)
D <- rbinom(M, 1, 0.01) # random D vector
names(D) <- gwas1.example$SNP

## The path to the directory where linkage disequilibrium (LD) information is stored.
LD.path <- "/your/path/to/UKB_array_SVD_eigen90_extraction"

## To speed up the test, we set a large lam.cut value.
res.sHDL <- sHDL(D, gwas1.example, LD.path, nthreads=4, stepwise=TRUE, lam.cut=10, Dr.path=NULL, mode="memory")
print(as.data.frame(res.sHDL))
```

## More tutorials

**For more detail** please see tutorials in [https://now2014.github.io/sHDL/](https://now2014.github.io/sHDL/).

## Citation

If you use the `sHDL` software, please cite:

Lan A. and Shen X. Modeling the genetic architecture of complex traits via stratified high-definition likelihood. (2024).

[Ning, Z., Pawitan, Y. & Shen, X. High-definition likelihood inference of genetic correlations across human complex traits. *Nat Genet* (2020)](https://www.nature.com/articles/s41588-020-0653-y).

## Contact

Reporting bugs by opening a new issue on this [Github page](https://github.com/now2014/sHDL/issues).

Sending email to authors:  [Ao Lan](mailto:lanao@mail2.sysu.edu.cn) or [Xia Shen](mailto:shenx@fudan.edu.cn).

