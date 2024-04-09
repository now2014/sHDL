#' Reduce the dimension of covariance matrix by converting z (D) to zr (Dr).
#'
#' @param LD.file Path to the \code{.rda} file where the Eigen decomposition of LD matrix is stored.
#' @param D A vector of genomic annotations with vector names of SNP IDs.
#' @param z A matrix of Z-scores with rownames of SNP IDs. Supporting multiple columns for multiple traits.
#' @param lam.cut Eigenvalue cutoff for LD matrices, default lam.cut = NULL, which means no cutoff.
#' @param Dr.path Path to the directory where the Dr matrices are stored, default Dr.path = NULL, which means do not store Dr to disk.
#' @param overwrite Whether to overwrite the existing Dr matrices, default overwrite = FALSE.
#' @param mode Whether to store Dr to disk or memory, default \code{mode = "disk"}. If \code{mode = "disk"}, \code{Dr} is stored to disk (path returned only) and lam are not returned. If \code{mode = "memory"}, \code{Dr} and \code{lam} are returned.
#' @param nthreads Number of threads to use for matrix operations, default \code{nthreads = 1}.
#' @param pattern Chromosome and picece pattern of LD files, default is \code{".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*"}.
#' @param norm.method The normalization method, either \code{"minmax"} (default), \code{"scaled"} or \code{"none"}. If \code{"minmax"}, the annotation weight vector D is normalized to [0, 1]. If \code{"scaled"}, the sum of normalized vector D is scaled to the number of annotated SNPs. If \code{"none"}, the annotation weight vector D is not normalized.
#' @return A list is returned with:
#' \itemize{
#' \item{Dr }{The reduct RDR matrix.}
#' \item{zr }{The reduct z-score vector.}
#' \item{lam }{The eigenvalues of LD matrix.}
#' \item{Md }{The sum of annotation weights of SNPs in given LD matrix.}
#' \item{M }{The number of SNPs in given LD matrix.}
#' }
#' @export
#' 
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads

sHDL.reduct.dim <- function(LD.path, z=NULL, D=NULL, lam.cut=NULL,
  Dr.path=NULL, overwrite=FALSE, mode=c("disk", "memory"), nthreads=1,
  pattern=".*chr(\\d{1,2})\\.(\\d{1,2})[_\\.].*", norm.method=c("minmax", "scaled", "none")){
  # use RhpcBLASctl to control the number of thread
  blas_set_num_threads(nthreads)
  omp_set_num_threads(nthreads)
  if(!is.null(Dr.path) && Dr.path!="" && !dir.exists(Dr.path))
    dir.create(Dr.path, recursive = TRUE)
  mode <- match.arg(mode)
  LD.files <- sHDL:::list.LD.ref.files(LD.path, suffix=".rda", pattern=pattern)

  if(!is.null(D)) D <- sHDL:::normD(D, LD.path, norm.method=norm.method, pattern=pattern)
  ref.data <- list()
  for(i in seq_along(LD.files)){
    LD.file <- LD.files[i]
    if(!is.null(Dr.path)) Dr.file <- paste0(Dr.path, "/", basename(LD.file))
    if(!is.null(Dr.path) && file.exists(Dr.file) && !overwrite && mode=="disk"){
      Dr <- Dr.file
      Md <- NULL
      lam <- NULL
    }else if(!is.null(Dr.path) && file.exists(Dr.file) && !overwrite && mode=="memory"){
      load(Dr.file) # Dr, lam, Md, M
    }else{
      Dr <- NULL
      Md <- NULL
      lam <- NULL
    }

    if(is.null(z) && !is.null(Dr)){
      ref.data[[i]] <- list(Dr=Dr, zr=NULL, lam=lam, Md=Md, M=M)
      next
    }

    load(LD.file) # V, lam, LDsc
    if(exists('U') && !exists('V')){
      V <- U
      rm('U')
    }
    if(!is.null(lam.cut)){
      idx <- lam > lam.cut
      if(sum(idx)==0) idx[1] <- TRUE
      V <- V[,idx]
      lam <- lam[idx]
    }
    snps.ref <- row.names(V)
    M <- length(snps.ref)

    if(!is.null(z)){
      z[is.na(z)] <- 0
      z[!is.finite(z)] <- 0
      Z <- matrix(0, nrow=M, ncol=ncol(z))
      row.names(Z) <- snps.ref
      int.snps <- intersect(snps.ref, row.names(z))
      Z[int.snps,] <- z[int.snps,]
      zr <- crossprod(V, Z)
    }else{
      zr <- NULL
    }

    if(!is.null(D) && is.null(Dr)){
      dD <- numeric(M)
      names(dD) <- snps.ref
      int.snps <- intersect(snps.ref, names(D))
      dD[int.snps] <- D[int.snps]
      Dr <- diag(lam) %*% t(V) %*% diag(dD) %*% V %*% diag(lam)
      Md <- sum(dD)

      if(!is.null(Dr.path) && (overwrite || !file.exists(Dr.file))) save(Dr, lam, Md, M, file=Dr.file)

      if(mode=="disk"){
        Dr <- Dr.file
      }
    }
    rm('V')
    ref.data[[i]] <- list(Dr=Dr, zr=zr, lam=lam, Md=Md, M=M)
  }
  return(ref.data)
}