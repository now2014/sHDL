#' Maximize the log likelihood to infer the heritability enrichment fold.
#'
#' @param ref.data A list of data matched to LD reference generated by \code{\link{sHDL.reduct.dim}}.
#' @param start.v A vector of starting values \code{c(fold, h2, intercept)} for optimization.
#' @param N The sample size of the GWAS.
#' @param output.file Where the log and results should be written.
#' @param log.file Where the log should be written.
#' @param stepwise Whether to estimate enrichment fold by estimating heritability and intercept first, default \code{stepwise = FALSE}. If \code{fix.h2} and \code{fix.intercept} are not \code{NULL}, \code{stepwise} will be overridden.
#' @param fix.h2 Whether to fix the heritability to \code{fix.h2} or estimate the heritability, default \code{fix.h2 = NULL}, which means estimate the heritability.
#' @param fix.intercept Whether to fix the intercept to \code{fix.intercept} or estimate the intercept, default \code{fix.intercept = NULL}, which means estimate the intercept.
#' @param lim Tolerance limitation to ensure the positive-definiteness of covariance matrices, default \code{lim = exp(-18)}.
#' @param verbose Whether to print the log on the console, default \code{verbose = FALSE}.
#' @param clust A cluster object generated by \code{\link{parallel::makeCluster}}.
#' @param lwr Lower bounds for \code{c(fold, h2, intercept)}. Default is \code{NULL}, which means \code{c(0, 0, 0.1)}.
#' @param upr Upper bounds for \code{c(fold, h2, intercept)}. Default is \code{NULL}, which means \code{c(M / Md, 1, 5)}, where \code{Md} is the sum of annotation weights and \code{M} is the total number of SNPs.
#' @param maxit Maximum number of iterations, default maxit = 1000.
#' @param pgtol Tolerance for convergence, default pgtol = 1e-3.
#' @return A list is returned with:
#' \itemize{
#' \item{time }{Time elapsed in seconds for optimization.}
#' \item{fold }{The estimated heritability enrichment fold.}
#' \item{h2 }{The estimated SNP-based heritability.}
#' \item{intercept }{The estimated intercept.}
#' \item{fold.se }{The standard error of the estimated heritability enrichment fold.}
#' \item{h2.se }{The standard error of the estimated heritability.}
#' \item{intercept.se }{The standard error of the estimated intercept.}
#' \item{fold.p }{P-value based on Wald test for the estimated heritability enrichment fold.}
#' \item{h2.p }{P-value based on Wald test for the estimated heritability.}
#' \item{intercept.p }{P-value based on Wald test for the estimated intercept.}
#' \item{stepwise }{Whether the optimization is done in a stepwise manner.}
#' \item{converged }{Whether the optimization converges.}
#' \item{message }{The message returned by \code{\link{optim}}.}
#' }
#' @export
#'
#' @importFrom dplyr filter distinct
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads

sHDL.optim <- function(
  ref.data, N, start.v=c(1, 0.1, 1), output.file=NULL, log.file="", stepwise=FALSE,
  fix.h2=NULL, fix.intercept=NULL, lim=exp(-18), verbose=FALSE, clust=NULL,
  lwr=NULL, upr=NULL, maxit=1000, pgtol=1e-3
  ){
  t0 <- Sys.time()
  ref.lam.M.Md <- lapply(ref.data, sHDL:::read.lam)
  ## ref.data: list of zr, Dr, lam
  if(!is.null(fix.intercept) & !is.null(fix.h2)) stepwise <- FALSE
  if(file.exists(log.file) & log.file!="") unlink(log.file)
  Md <- sum(unlist(lapply(ref.lam.M.Md , function(x) x$Md)))
  M <- sum(unlist(lapply(ref.lam.M.Md , function(x) x$M)))
  if(!is.null(fix.intercept) & length(start.v)==3) start.v <- start.v[-3]
  if(!is.null(fix.h2)) start.v <- start.v[-2]
  if(is.null(lwr)) lwr <- c(0, 0, 0.1)
  if(is.null(upr)) upr <- c(M / Md, 1, 5)

  
  if(stepwise){
    # optim h2 & intercept first
    lam.all <- unlist(lapply(ref.lam.M.Md, function(x) x$lam))
    zr.all <- unlist(lapply(ref.data, function(x) x$zr))
    opt.first <- optim(
      start.v[-1], sHDL:::log.lik.HDL,
      lam=lam.all, zr=zr.all, N=N, M=M, fix.intercept=fix.intercept, lim=lim,
      method ="L-BFGS-B", lower=lwr[-1], upper=upr[-1], 
      hessian=TRUE, control=list(maxit=maxit, fnscale=-1, pgtol=pgtol)
    )

    se <- tryCatch(
      sqrt(diag(solve(-opt.first$hessian))),
      error = function(e){
        warn.msg <- "Hessian matrix is NOT invertible.\n"
        sHDL:::log.msg(warn.msg, log.file, type="warning")
        NA
      }
    )

    h2 <- opt.first$par[1]
    h2.se <- se[1]
    if(is.null(fix.intercept)){
      intercept <- opt.first$par[2]
      intercept.se <- se[2]
      fix.intercept <- intercept
    }else{
      intercept <- fix.intercept
      intercept.se <- NA
    }
    fix.h2 <- h2
    start.v <- start.v[1]
    lwr <- lwr[1]
    upr <- upr[1]
  }


  opt <- optim(
    start.v, log.lik.wg,
    ref.data=ref.data, Md=Md, M=M, N=N,
    verbose=verbose, method ="L-BFGS-B", lower=lwr, upper=upr,
    log.file=log.file, h2=fix.h2, intercept=fix.intercept, lim=lim, clust=clust,
    hessian=TRUE, control=list(maxit=maxit, fnscale=-1, pgtol=pgtol)
  )

  time <- as.numeric(Sys.time() - t0, units="secs")
  sHDL:::log.msg(sprintf("Optimization done in %.3f seconds.\n", time), log.file)
  covergence <- opt$convergence
  msg <- opt$message

  ## convert hessian to SE & p-value
  se <- tryCatch(
    sqrt(diag(solve(-opt$hessian))),
    error = function(e){
      warn.msg <- "Hessian matrix is NOT invertible.\n"
      sHDL:::log.msg(warn.msg, log.file, type="warning")
      NA
    }
  )
  se <- as.numeric(se)
  fold <- opt$par[1]
  fold.se <- se[1]
  

  if(is.null(fix.h2)){
    h2 <- opt$par[2]
    h2.se <- se[2]
  }else if(!stepwise){
    h2 <- fix.h2
    h2.se <- NA
  }

  if(is.null(fix.intercept)){
    intercept <- opt$par[3]
    intercept.se <- se[3]
  }else if(!stepwise){
    intercept <- fix.intercept
    intercept.se <- NA
  }

  fold.p <- pchisq(((fold-1)/fold.se)^2, df = 1, lower.tail = FALSE)
  intercept.p <- pchisq((intercept/intercept.se)^2, df = 1, lower.tail = FALSE)
  h2.p <- pchisq((h2/h2.se)^2, df = 1, lower.tail = FALSE)

  res <- list(
    time=time, fold=fold, h2=h2, intercept=intercept,
    fold.se=fold.se, h2.se=h2.se, intercept.se=intercept.se,
    fold.p=fold.p, h2.p=h2.p, intercept.p=intercept.p, stepwise=stepwise,
    converged=covergence==0, optim.message=msg
  )
  if(!is.null(output.file)){
    write.table(res, output.file, sep="\t", row.names=F, quote=F, col.names=T)
  }
  if(!is.null(log.file)){
    write.table(unlist(res), log.file, sep="\t", quote=F, row.names = T, col.names = F, append = TRUE)
  }
  if(!is.null(clust)){
    tmp <- tryCatch({
      suppressWarnings(stopCluster(clust))
    }, error = function(e){
      NULL
    })
  }
  return(res)
}