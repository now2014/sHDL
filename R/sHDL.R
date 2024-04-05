#' Stratified High-Definition Likelihood (sHDL) Inference of Heritability Enrichment
#'
#' The function return estimation and standard error of heritability enrichment fold of trait based on GWAS summary statistics.
#'
#' @param D A vector of genomic annotations with vector names of SNP IDs.
#' @param gwas.df A data frame including GWAS summary statistics of genetic variants for a trait.
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS.
#' @param LD.path Path to the directory where linkage disequilibrium (LD) information is stored, compatible with HDL LD reference panels.
#' @param output.file Where the log and results should be written. If you do not specify a file, the log will be printed on the console.
#' @param nthreads Number of threads to use, default \code{nthreads = 1}.
#' @param stepwise Whether to estimate enrichment fold by estimating heritability and intercept first, default \code{stepwise = FALSE}. If \code{fix.h2} and \code{fix.intercept} are not NULL, stepwise will be overridden.
#' @param fill.missing.N If the sample size is missing in the GWAS summary statistics, fill.missing.N is used to fill the missing values.
#' @param lim Tolerance limitation to ensure the positive-definiteness of covariance matrices, default \code{lim = exp(-18)}.
#' @param lam.cut Eigenvalue cutoff for LD matrices, default \code{lam.cut = NULL}, which means no cutoff.
#' @param Dr.path Path to the directory where the Dr matrices are stored, default \code{Dr.path = "./Dr"}.
#' @param overwrite Whether to overwrite the existing Dr matrices, default \code{overwrite = FALSE}.
#' @param verbose Whether to print the log on the console, default \code{verbose = FALSE}.
#' @param fix.h2 Whether to fix the heritability to \code{fix.h2} or estimate the heritability, default \code{fix.h2 = NULL}, which means estimate the heritability.
#' @param fix.intercept Whether to fix the intercept to \code{fix.intercept} or estimate the intercept, default \code{fix.intercept = NULL}, which means estimate the intercept.
#' @param maxit Maximum number of iterations, default \code{maxit = 1000}.
#' @param pgtol Tolerance for convergence, default \code{pgtol = 1e-3}.
#' @param start.v A vector of starting values \code{c(fold, h2, intercept)} for optimization.
#' @param mode Whether to store Dr to disk or memory, default \code{mode = "disk"}. If \code{mode = "disk"}, \code{Dr} is stored to disk (path returned only) and lam are not returned. If \code{mode = "memory"}, \code{Dr} and \code{lam} are returned.
#' @param pattern Chromosome and picece pattern of LD files, default is \code{".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"}.
#' @note Users can download the precomputed eigenvalues and eigenvectors of LD correlation matrices for European ancestry population. The download link can be found at https://github.com/zhenin/HDL/wiki/Reference-panels
#' These are the LD matrices and their eigen-decomposition from 335,265 genomic British UK Biobank individuals. Three sets of reference panel are provided:
#' 1) 1,029,876 QCed UK Biobank imputed HapMap3 SNPs. The size is about 33 GB after unzipping. Although it takes more time, using the imputed panel provides more accurate estimates of genetic correlations.
#' Therefore if the GWAS includes most of the HapMap3 SNPs, then we recommend using the imputed reference panel.
#' 2) 769,306 QCed UK Biobank imputed HapMap2 SNPs. The size is about 18 GB after unzipping.If one of your GWAS includes most of the HapMap 2 SNPs, but many SNPs (more than 1%) in the above HapMap 3 reference panel are absent,
#' then this HapMap2 panel is more proper to be used for HDL.
#' 3) 307,519 QCed UK Biobank Axiom Array SNPs. The size is about 7.5 GB after unzipping.
#'
#' @return A list is returned with:
#' \itemize{
#' \item{time }{The time used for optimization.}
#' \item{fold }{The estimated heritability enrichment fold.}
#' \item{h2 }{The estimated SNP-based heritability.}
#' \item{intercept }{The estimated intercept.}
#' \item{fold.se }{The standard error of the estimated heritability enrichment fold.}
#' \item{h2.se }{The standard error of the estimated heritability.}
#' \item{intercept.se }{The standard error of the estimated intercept.}
#' \item{fold.p }{P-value based on Wald test for the estimated heritability enrichment fold.}
#' \item{h2.p }{P-value based on Wald test for the estimated heritability.}
#' \item{intercept.p }{P-value based on Wald test for the estimated intercept.}
#' \item{converged }{Whether the optimization converges.}
#' \item{message }{The message returned by \code{\link{optim}}.}
#' }
#'
#' @author Ao Lan, Xia Shen
#'
#' @references
#' Lan A and Shen X (2024). Modeling the Genetic Architecture of Complex Traits via Stratified High-Definition Likelihood.
#'
#' @importFrom dplyr filter
#' @importFrom parallel makeCluster stopCluster clusterExport
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads
#'
#' @examples
#'\dontrun{
#' ## The GWAS summary statistics for birth weight loaded from HDL package.
#' data(gwas1.example, package='HDL')
#' M <- nrow(gwas1.example)
#' set.seed(1234)
#' D <- rbinom(M, 1, 0.01) # random D vector
#' names(D) <- gwas1.example$SNP
#'
#' ## The path to the directory where linkage disequilibrium (LD) information is stored.
#' LD.path <- "path/to/UKB_array_SVD_eigen90_extraction"
#' 
#' ## To speed up the test, we set a large lam.cut value.
#' res.sHDL <- sHDL(D, gwas1.example, LD.path, nthreads=1, stepwise=TRUE, lam.cut=10, Dr.path="./Dr", verbose=F)
#' as.data.frame(res.sHDL)
#' system("rm -rf ./Dr") # remove the Dr directory
#' }
#' @export
#'

sHDL <-function(
  D, gwas.df, LD.path, output.file = NULL, nthreads = 1, stepwise = TRUE,
  fill.missing.N = NULL, lim = exp(-18), lam.cut = NULL,
  Dr.path = "./Dr", overwrite = FALSE, verbose = FALSE,
  fix.h2 = NULL, fix.intercept = NULL, maxit=1000,
  pgtol=1e-3, start.v = c(1, 0.1, 1), mode=c("disk", "memory"),
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"){

  mode <- match.arg(mode)
  if(!is.null(output.file)){
    log.file <- paste0(output.file, ".log")
    if(file.exists(log.file)) unlink(log.file)
  }else{
    log.file <- ""
  }

  if(nthreads > 1){
    clust <- makeCluster(nthreads, outfile=log.file)
  }else{
    clust <- NULL
  }

  ## min-max normalization on D
  minv <- min(D, na.rm = T)
  maxv <- max(D, na.rm = T)
  D <- (D - minv)/(maxv - minv)
  D <- D[setdiff(names(D), D[duplicated(names(D))])] ## remove duplicates

  message("Formating GWAS summary statistics ...\n")
  gwas.df <- format.gwas(gwas.df, LD.path, fill.missing.N, log.file, pattern)

  N <- median(gwas.df$N, na.rm=TRUE)
  z <- gwas.df$Z
  z <- matrix(z, ncol=1)
  rownames(z) <- gwas.df$SNP

  message("Converting z (D) to zr (Dr) ...\n")
  ref.data <- sHDL:::sHDL.reduct.dim(LD.path, z=z, D=D, lam.cut=lam.cut,
    Dr.path=Dr.path, overwrite=overwrite, mode=mode,
    nthreads=nthreads, pattern=pattern)
  M <- sum(unlist(lapply(ref.data, function(x) x$M)))
  Md <- sum(unlist(lapply(ref.data, function(x) x$Md)))

  message("Optimizing ...\n")
  res <- sHDL.optim(
    ref.data, N, start.v, output.file, log.file,
    stepwise, fix.h2, fix.intercept, lim, verbose, clust,
    lwr=NULL, upr=NULL, maxit=maxit, pgtol=pgtol
  )
  return(res)
}