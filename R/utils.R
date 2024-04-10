#' @title Format GWAS data
#' @noRd
#' @keywords internal
format.gwas <- function(gwas.df, LD.path, fill.missing.N = NULL, log.file = "",
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"){

  bim.files <- sHDL:::list.LD.ref.files(LD.path, suffix='\\.bim', pattern=pattern)
  bim <- do.call(rbind, lapply(
    bim.files, read.table, header=F,
    col.names=c('CHR', 'SNP', 'CM', 'POS', 'A1', 'A2')
  ))
  Mref <- nrow(bim)
  
  gwas.df <- filter(gwas.df, SNP %in% bim$SNP)

  gwas.df$A1 <- toupper(as.character(gwas.df$A1))
  gwas.df$A2 <- toupper(as.character(gwas.df$A2))

  if(!("Z" %in% colnames(gwas.df))){
    if(("b" %in% colnames(gwas.df)) && ("se" %in% colnames(gwas.df))){
      if(abs(median(gwas.df$b) - 1) < 0.1){
        cat("Taking log(b) in GWAS 1 because b is likely to be OR in stead of log(OR). \n")
        if(log.file != ""){
          cat("Taking log(b) in GWAS 1 because b is likely to be OR in stead of log(OR). \n", file = log.file, append = T)
        }
        gwas.df$Z <- log(gwas.df$b) / gwas.df$se
      } else{
        gwas.df$Z <- gwas.df$b / gwas.df$se
      }
    } else{
      error.message <- "Z is not available in GWAS 1, meanwhile either b or se is missing. Please check."
      if(log.file != ""){
        cat(error.message, file = log.file, append = T)
      }
      stop(error.message)

    }
  }
  gwas.df$Z[!is.finite(gwas.df$Z)] <- NA
  k1.0 <- length(unique(gwas.df$SNP))

  if(is.null(fill.missing.N)){
    gwas.df <- filter(gwas.df, !is.na(Z), !is.na(N))

  } else if(fill.missing.N == "min"){
    gwas.df <- filter(gwas.df, !is.na(Z))
    gwas.df$N[is.na(gwas.df$N)] <- min(gwas.df$N, na.rm = T)
  } else if(fill.missing.N == "max"){
    gwas.df <- filter(gwas.df, !is.na(Z))
    gwas.df$N[is.na(gwas.df$N)] <- max(gwas.df$N, na.rm = T)
  } else if(fill.missing.N == "median"){
    gwas.df <- filter(gwas.df, !is.na(Z))
    gwas.df$N[is.na(gwas.df$N)] <- median(gwas.df$N, na.rm = T)
  } else{
    error.message <- "If given, the argument fill.missing.N can only be one of below: 'min', 'max', 'median'."
    if(log.file != ""){
      cat(error.message, file = log.file, append = T)
    }
    stop(error.message)
  }

  k1 <- length(unique(gwas.df$SNP))
  k1.percent <- paste("(",round(100*k1 / Mref, 2), "%)", sep="")

  if(k1 < Mref*0.99){
    warn.message <- "More than 1% SNPs in reference panel are missed in the GWAS. This may generate bias in estimation. Please make sure that you are using correct reference panel.  \n"
    if(log.file != ""){
      cat(warn.message, file = log.file, append = T)
    }else{
      warning(warn.message)
    }
  }

  if(log.file != ""){
    cat(k1, "out of", Mref, k1.percent, "SNPs in reference panel are available in the GWAS."," \n", file = log.file, append = T)
  }else{
    cat(k1, "out of", Mref, k1.percent, "SNPs in reference panel are available in the GWAS."," \n")
  }
  N <- median(gwas.df$N, na.rm=TRUE)

  ## remove duplicates
  dup.snps <- gwas.df$SNP[duplicated(gwas.df$SNP)]
  gwas.df <- filter(gwas.df, !(SNP %in% dup.snps))

  ## match direction of z-scores
  ref.A1 <- bim$A1
  names(ref.A1) <- bim$SNP
  ref.A1 <- ref.A1[gwas.df$SNP]
  idx.direc <- ifelse(ref.A1 == gwas.df$A1, 1, -1)
  gwas.df$Z <- gwas.df$Z * idx.direc
  return(gwas.df)
}

#' @title List LD reference files
#' @noRd
#' @keywords internal
list.LD.ref.files <- function(LD.path, suffix = ".rda", full.names = TRUE,
                              pattern = ".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"){
  segfiles <- list.files(LD.path, pattern = suffix, full.names = F)
  segfiles <- segfiles[grepl(pattern = pattern, segfiles)]
  chrom <- gsub(pattern, "\\1", segfiles)
  piece <- gsub(pattern, "\\2", segfiles)
  chrom <- as.numeric(chrom)
  piece <- as.numeric(piece)
  idx <- order(chrom, piece)
  segfiles <- segfiles[idx]
  if(full.names) segfiles <- paste0(LD.path, "/", segfiles)
  return(segfiles)
}

#' @title read lam, M, Md from Dr.file
#' @noRd
#' @keywords internal
read.lam <- function(refd){
  lam <- refd$lam
  M <- refd$M
  Md <- refd$Md
  if(is.null(lam) | is.null(Md) | is.null(M)){
    load(refd$Dr)
  }
  return(list(lam=lam, M=M, Md=Md))
}

#' @title log likelihood function without enrichment fold parameter
#' @noRd
#' @keywords internal
log.lik.HDL <- function(param, lam, zr, N, M, lim=exp(-18), fix.intercept=NULL){
  h2 <- param[1]
  if(!is.null(fix.intercept)){
    intercept <- fix.intercept
  }else{
    intercept <- param[2]
  }
  alpha <- intercept * lam + N/M * h2 * lam**2
  alpha <- pmax(alpha, lim)
  logdet <- sum(log(alpha))
  quadra <- sum(zr**2/alpha)
  lnL <- -0.5 * (logdet + quadra)
  return(lnL)
}

#' @title log likelihood function for piece
#' @noRd
#' @keywords internal
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads
log.lik <- function(refd, h20, h2d, intercept,
  N, M, Md, lim=exp(-18)){
  zr <- refd$zr
  Dr <- refd$Dr
  lam <- refd$lam
  if(is.null(lam)){
    e1 <- new.env()
    load(Dr, envir = e1)
    Dr <- e1$Dr
    lam <- e1$lam
  }
  alpha <- intercept * lam + N*h20/M * lam**2
  alpha <- pmax(alpha, lim)

  if(length(Dr) == 1){
    Dr <- as.vector(Dr)
    sigma <- alpha + N*h2d/Md * Dr
    logdet <- log(sigma)
    quadra <- zr**2/sigma
    lnL <- -0.5 * (logdet + quadra)
    return(lnL)
  }

  # use RhpcBLASctl to control the number of threads
  blas_set_num_threads(1)
  omp_set_num_threads(1)
  sigma <- diag(alpha) + N*h2d/Md * Dr
  L <- tryCatch({
    chol(sigma)
  }, error = function(e){
    # cat("Error in chol(sigma) \n")
    # cat("min(alpha): ", min(alpha), "\n")
    # cat("h20: ", h20, "\n")
    # cat("h2d: ", h2d, "\n")
    # cat("intercept: ", intercept, "\n")
    NULL
  })
  if(is.null(L)) return(-1e18)
  logdet <- 2 * sum(log(diag(L)))
  quadra <- crossprod(zr,
    forwardsolve(L, backsolve(L, zr, transpose=T), upper.tri=T))
  logdet <- as.vector(logdet)
  quadra <- as.vector(quadra)
  lnL <- -0.5 * (logdet + quadra)
  # cat(h20, h2d, logdet, quadra, lnL, '\n')
  return(lnL)
}

#' @title log likelihood function for whole genome
#' @noRd
#' @keywords internal
log.lik.wg <- function(param, ref.data, Md, M, N,
  log.file="", h2=NULL, intercept=NULL, fold=NULL,
  verbose=FALSE, lim=exp(-18), clust=NULL){
  t0 <- Sys.time()

  fold <- param[1]
  if(is.null(h2)) h2 <- param[2]
  if(is.null(intercept)) intercept <- param[3]

  h2d <- Md*(fold-1)/(M - Md)*h2
  h20 <- h2 - h2d

  if(!is.null(clust)){
    lnL <- parLapply(cl=clust, ref.data, sHDL:::log.lik,
      N=N, M=M, Md=Md, h20=h20, h2d=h2d, intercept=intercept, lim=lim)
  }else{
    lnL <- lapply(ref.data, sHDL:::log.lik,
      N=N, M=M, Md=Md, h20=h20, h2d=h2d, intercept=intercept, lim=lim)
  }

  lnL <- sum(unlist(lnL))
  time <- as.numeric(Sys.time() - t0, units="secs")
  if(verbose) cat(sprintf(
    "fold: %.3f h2: %.3f intercept: %.3f lnL: %.3f time: %.3f \n",
    fold, h2, intercept, lnL, time
  ), file=log.file, append=T)
  return(lnL)
}

#' @title D normalization
#' @noRd
#' @keywords internal
normD <-function(
  D, LD.path, log.file="", norm.method=c("minmax", "scaled", "none"),
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"){

  method <- match.arg(norm.method)
  bim.files <- sHDL:::list.LD.ref.files(LD.path, suffix='\\.bim', pattern=pattern)
  bim <- do.call(rbind, lapply(
    bim.files, read.table, header=F,
    col.names=c('CHR', 'SNP', 'CM', 'POS', 'A1', 'A2')
  ))
  ref.snps <- bim$SNP
  M <- nrow(bim)
  dup.snps <- names(D)[duplicated(names(D))]
  D <- D[setdiff(names(D), dup.snps)] ## remove duplicates

  int.snps <- intersect(ref.snps, names(D))
  dD <- rep(0, M)
  names(dD) <- ref.snps
  dD[int.snps] <- D[int.snps]
  dD[is.na(dD)] <- 0
  D <- dD
  minv <- min(D)

  if(method == "minmax"){
    maxv <- max(D)
    D <- (D - minv)/(maxv - minv)
    Md <- sum(D > 0)
  }else if(method == "scaled"){
    if(minv < 0){
      D <- D - minv
    }
    Md <- sum(D > 0)
    if(Md == M){
      D <- D - min(D)
    }
    D <- Md / sum(D) * D
  }else if(method=="none"){
    if(minv < 0) warn.msg <- "The annotation weights contain negative values, which may cause bias in the estimation.\n"
    if(log.file != "" && !is.null(warn.msg)){
      cat(warn.msg, file = log.file, append = T)
    }else if(!is.null(warn.msg)){
      warning(warn.msg)
    }
    Md <- sum(D != 0)
    msg <- sprintf("No normalization applied on %d (%.3f%%) annotated variants. The theoretical upper boundary for enrichment fold is M / Md = %.3f.\n",
      Md, Md/M, M/sum(D))
    if(log.file != ""){
      cat(msg, file = log.file, append = T)
    }else{
      cat(msg)
    }
    return(D)
  }else{
    stop("Unknown normalization method.")
  }

  msg <- sprintf("Applied `%s` weight nomalization on %d (%.3f%%) annotated variants.\n",
    method, Md, Md/M*100)
  msg <- sprintf("%sThe theoretical upper boundary for enrichment fold is M / Md = %.3f", msg, M/sum(D))
  if(sum(D)/M > 0.9 && method == "scaled"){
    warn.msg <- paste0(msg, ", which is very close to 1, please consider the `minmax` option or reducing the dense of annotations.\n")
    msg <- NULL
  }else if(sum(D)/M > 0.9){
    warn.msg <- paste0(msg, ", which is very close to 1, please consider reducing the dense of annotations.\n")
    msg <- NULL
  }else{
    msg <- paste0(msg, ".\n")
    warn.msg <- NULL
  }
  
  if(log.file != ""){
    if(!is.null(msg)) cat(msg, file = log.file, append = T)
    if(!is.null(warn.msg)) cat(warn.msg, file = log.file, append = T)
  }else{
    if(!is.null(msg)) cat(msg)
    if(!is.null(warn.msg)) warning(warn.msg)
  }
  return(D)
}