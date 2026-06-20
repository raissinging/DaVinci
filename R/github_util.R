require(rsvd)
require(tripack)
require(ggpubr)



rotateSVD=function(svdres){
  upos=svdres$u
  uneg=svdres$u
  upos[upos<0]=0
  uneg[uneg>=0]=0
  uneg=-uneg
  sumposu=colSums(upos)
  sumnegu=colSums(uneg)


  for(i in 1:ncol(svdres$u)){
    if(sumnegu[i]>sumposu[i]){
      svdres$u[,i]=-svdres$u[,i]
      svdres$v[,i]=-svdres$v[,i]
    }
  }
  svdres
}#rotateSVD




tscale <- function(mat){
  return(t(scale(t(mat))))
}#tscale



normalize_barcode_sums_to_median <- function(gbm){
  bc_sums <- Matrix::colSums(gbm)
  median_sum <- median(bc_sums)
  return(sweep(gbm,2,median_sum/bc_sums, '*'))
}#normalize_barcode_sums_to_median




#' Normalize the spatial RNA-seq
#' @export
gene_normalization <- function(gene.exp, frac.thr = 0.95, MT.remove = T, median.norm = T){
  #use_gene <- which(apply(gene.exp==0,1,sum)/ncol(gene.exp) < frac.thr)

  use_gene <- which(Matrix::rowSums(gene.exp==0)/ncol(gene.exp) < frac.thr)
  if (is.data.frame(gene.exp)){
    gene.exp <- as.matrix(gene.exp[use_gene,])
  }else{
    gene.exp <- gene.exp[use_gene,]
  }#else

  #remove MT- genes
  if (MT.remove){
    use_gene <- which(!grepl("MT-", rownames(gene.exp)))
    gene.exp <- gene.exp[use_gene,]

    use_gene <- which(!grepl("mt-", rownames(gene.exp)))
    gene.exp <- gene.exp[use_gene,]
  }#if

  #check standard deviation
  gene.exp <- gene.exp[which(apply(gene.exp,1,sd)>0),]

  #remove all zero cells
  gene.exp <- gene.exp[, which( Matrix::colSums(gene.exp)>0)]

  if (median.norm){
    gene.exp <- normalize_barcode_sums_to_median(gene.exp)
  }#median.norm

  gene.exp.log <- log2(gene.exp+1)

  #tscale
  #gene.exp <- tscale(gene.exp[apply(gene.exp,1,sd)>0,])
  return(list(gene.exp = gene.exp, gene.exp.log = gene.exp.log))
}#gene_normalization





#http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/
#gene.exp <- as.matrix(gene.exp)
#tf = t( t(gene.exp) / Matrix::colSums(gene.exp))
#tf <- log1p(tf*1e5)
#idf = log(1 + ncol(gene.exp) / Matrix::rowSums(gene.exp))
#gene.exp <- tf*idf

#https://satijalab.org/seurat/archive/v3.1/atacseq_integration_vignette
#directly based on the signac/R/preprocessing.R and https://stuartlab.org/signac/reference/runtfidf

#' Normalize the spatial ATAC-seq
#' @export 
LSI_normalization <- function(atac, method = "1", scale.factor = 1e5, peak.filter.thr = NULL){

  #filter out no count peaks
  #############################################
  include.index <- which(Matrix::rowSums(atac)>0)
  atac <- atac[include.index,]


  #filter out peaks based on fractions
  #############################################
  if (!is.null(peak.filter.thr)){
    include.index <- which(apply(atac==0,1,sum)/ncol(atac) < peak.filter.thr)
    atac <- atac[include.index,]
  }#if


  atac <- atac[which(apply(atac,1,sd)>0),]

  #binarize the matrix
  #############################################
  atac[atac >0] <- 1

  #TF calculation
  #############################################
  if (method == 4 ){
    tf <- atac
  }else{
    tf <- Matrix::t(Matrix::t(atac)/Matrix::colSums(atac))
  }#else


  #IDF calculation
  #############################################
  idf <- ncol(atac)/Matrix::rowSums(atac)

  if (method==2){
    idf <- log1p(idf)
  }else if (method==3){
    tf <- log1p(tf*scale.factor)
    idf <- log1p(idf)
  }#else if

  #final product
  #############################################
  idf.diag <- Matrix::sparseMatrix(i = 1:length(idf),
                                   j = 1:length(idf),
                                   x = idf,
                                   dims = c(length(idf), length(idf)))
  #idf.diag <- diag(idf)
  output <- idf.diag %*%tf

  if (method==1){
    output <- log1p(output*scale.factor)
  }#if


  #set NA to 0
  output[is.na(output)] <- 0
  rownames(output) <- rownames(atac)
  colnames(output) <- colnames(atac)

  return(output)
}#LSI_normalization




#wrapper to combine everything
#' Wrapper to prepare the input
#' @export
preprocess <- function(mat, 
                       coor, 
                       type = "rna", 
                       graph.opt = "Tri.mesh", 
                       frac.thr = 0.95, 
                       MT.remove = T, 
                       median.norm = T, 
                       LSI.method = "3",  
                       LSI.scale.factor = 1e5, 
                       LSI.peak.filter.thr = NULL){

  #make sure aligned
  if (ncol(mat)!=nrow(coor)){
    stop("Dimensions don't align.")
  }#if

  if (!all(colnames(mat)==rownames(coor))){
    stop("Spot names don't match")
  }#if


  #processing the mat
  ################################################
  if (type %in% c("rna", "protein", "rna_raw")){

    tmp <- gene_normalization(mat, 
                              frac.thr = frac.thr, 
                              MT.remove = MT.remove, 
                              median.norm = median.norm)

    if (type=="rna_raw"){
      mat <- tmp$gene.exp
    }else{
      mat <- tmp$gene.exp.log  
    }#else
    

  }else if (type=="atac"){

    mat <- LSI_normalization(mat, 
                             method = LSI.method, 
                             scale.factor = LSI.scale.factor, 
                             peak.filter.thr = LSI.peak.filter.thr)

  }#else

  message(paste0(type, " normalization finishes."))

  #align the mat and coor again
  ################################################
  sids.included <- intersect(rownames(coor), colnames(mat))
  coor <- coor[sids.included,]
  mat <- mat[,sids.included]


  #processing the coor
  ################################################
  temp <- L_generate(coor, 
                     opt = graph.opt)
  L <- temp$L

  message(paste0("Graph generation finishes."))

  return(list(mat = Matrix::as.matrix(mat), 
              coor = coor, 
              L = L, 
              L.visual = temp$L.visual))
}#preprocess










#' Generate the Laplacian input
#' @import tripack
#' @import ggpubr
#' @export
L_generate <- function(coor, 
                       opt = "grid", 
                       dist.thr = 2, 
                       num.of.neighbor = NULL,
                       verbose = F){

  if (is.null(rownames(coor))){
    stop("Coor should have unique row names")
  }else if (!all(sort(colnames(coor))==c("array_col", "array_row"))){
    stop("colnames should be named as 'array_col' and 'array_row'")
  }#else if

  if (!opt %in% c("grid", "Tri.mesh", "nearest.neighbor")){
    stop("Opt is not valid")
  }else if (opt == "grid"){
    if (verbose){
      message(paste0("Your distance cutff is ", dist.thr))
    }#if verbose
  }#opt


  #initalize
  x <- c()
  y <- c()
  x.end <- c()
  y.end <- c()
  A <- matrix(0, nrow = nrow(coor), ncol = nrow(coor))
  rownames(A) <- rownames(coor)
  colnames(A) <- rownames(coor)

  if (opt=="grid"){
    coor.dist <- dist(coor)
    coor.dist <- as.matrix(coor.dist)

    for (i in 1: (nrow(coor.dist)-1)){
      temp <- which(coor.dist[i, (i+1):ncol(coor.dist)] <= dist.thr)
      if (length(temp)>0){
        temp.names <- colnames(coor.dist)[(i+1):ncol(coor.dist)]
        temp.names <- temp.names[temp]

        A[rownames(coor.dist)[i], temp.names] <- 1
        A[temp.names, rownames(coor.dist)[i]] <- 1
        x <- c(x, rep(coor[i,1], length(temp)))
        y <- c(y, rep(coor[i,2], length(temp)))
        x.end <- c(x.end, coor[temp.names,1])
        y.end <- c(y.end, coor[temp.names,2])
      }#if
    }#for i


  }else if (opt == "Tri.mesh"){
    tri.mesh.res <- tripack::tri.mesh(x = coor[,"array_row"], y = coor[,"array_col"])
    neighbor.info <- tripack::neighbours(tri.mesh.res)

    for (i in 1: length(neighbor.info)){
      A[i , neighbor.info[[i]] ] <- 1
      A[neighbor.info[[i]], i] <- 1
      x <- c(x, rep(coor[i,1], length(neighbor.info[[i]])))
      y <- c(y, rep(coor[i,2], length(neighbor.info[[i]])))
      x.end <- c(x.end, coor[neighbor.info[[i]],1])
      y.end <- c(y.end, coor[neighbor.info[[i]],2])
    }#for i


  }else if (opt == "nearest.neighbor"){
    if (is.null(num.of.neighbor)){
      stop("Number of neighbors is not set")
    }#if

    coor.dist <- dist(coor)
    coor.dist <- as.matrix(coor.dist)

    if (num.of.neighbor > nrow(coor)){
      num.of.neighbor <- nrow(coor)
      message(paste0("The number of neighbors has been adjusted to be ", num.of.neighbor))
    }#if

    for (i in 1:nrow(coor)){
      tt <- setdiff(1:nrow(coor), i)
      temp <- order( coor.dist[i, tt] , decreasing = F)[1:num.of.neighbor]
      temp <- tt[temp]

      A[i, temp] <- 1
      A[temp, i] <- 1
      x <- c(x, rep(coor[i,1], length(temp)))
      y <- c(y, rep(coor[i,2], length(temp)))
      x.end <- c(x.end, coor[temp,1])
      y.end <- c(y.end, coor[temp,2])
    }#for i

  }#else if


  D <-  diag(apply(A,1,sum))
  L <- D-A

  #L.visual: visualization
  dat.point <- data.frame(x = coor[,"array_row"], y=coor[,"array_col"])
  dat.edge <- data.frame(x=x,y=y, xend = x.end, yend = y.end)

  L.visual <- ggplot()+
              geom_point(data = dat.point, aes(x=x, y=y))+
              geom_segment(data = dat.edge, aes(x=x,y=y, xend =xend, yend = yend))+
              theme_pubclean()
  return(list(L = L, L.visual = L.visual))
}#L_generate





#' Main function with a fixed L4 parameter
#' @export
manifoldDecomp=function(Y, 
                        L, 
                        k,
                        svdres=NULL, 
                        L1=NULL, 
                        L2=NULL, 
                        L4 = NULL,
                        max.iter=200, 
                        tol=5e-6, 
                        trace=F,
                        B=NULL, 
                        B.random = F,
                        scale=1,  
                        adaptive.frac=0.05, 
                        adaptive.iter=30,
                        verbose = T,
                        random.seed = 123){

  round2=function(x){signif(x,4)}
  getT=function(x){-quantile(x[x<0], adaptive.frac)}



  pos.adj=3

  ng=nrow(Y)
  ns=ncol(Y)

  Bdiff=Inf
  BdiffTrace=double()
  BdiffCount=0
  message("****")

  if(is.null(svdres)){

    message("Computing SVD")
    set.seed(random.seed)
    svdres=rsvd(Y, k = k)

    svdres=rotateSVD(svdres)

    #  show(svdres$d[k])
  }#svdres


  #L1
  #######################################################
  if(is.null(L1)){
    L1=svdres$d[k]*scale
    if(!is.null(pos.adj)){
      L1=L1/pos.adj
    }#if
  }#if is.null(L1)


  #L2
  #######################################################
  if(is.null(L2)){
    L2=svdres$d[k]*scale
  }#if is.null(L2)

  #L1=svdres$d[k]/2*scale
  print(paste0("L1 is set to ",L1))
  print(paste0("L2 is set to ",L2))


  #B
  #######################################################
  if(is.null(B)){
    #initialize B with svd
    if (verbose){
      message("Init")
    }#if
    
    B=t(svdres$v[1:ncol(Y), 1:k]%*%diag(sqrt(svdres$d[1:k])))
  }else{
    
    if (verbose){
      message("B given")
    }#if verbose
    
  }#B initialization


  if (B.random) {
    
    if (verbose){
      message("using random start")
    }#if    
    
    set.seed(random.seed)
    B = t(apply(B, 1, sample))
  }#is.null


  right <- L4*t(L)+L4*L
  right.shur <- rcpp_shur(right)
  message("shur done")


  #updates
  #######################################################
  for ( i in 1:max.iter){

    #print(i)

    #Z update
    ######################################
    Zraw=Z=(Y%*%t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    if(i>=adaptive.iter && adaptive.frac>0){


      cutoffs=apply(Zraw,2, getT)

      for(j in 1:ncol(Z)){
        Z[Z[,j]<cutoffs[j],j]=0
      }#for j
    }else{
      Z[Z<0]=0
    }#else


    #B update
    ######################################
    oldB=B
    #B=solve(crossprod(Z)+L2*diag(k))%*%(t(Z)%*%Y)
    left <- 2*crossprod(Z)+2*L2*diag(k)
    total <- 2*t(Z) %*% Y

    B <- sylvester_pre(left, right.shur$U, right.shur$S, total)


    #update error
    ######################################
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    err0=sum((Y-Z%*%B)^2)+sum((Z)^2)*L1+sum(B^2)*L2
    if(trace){
      message(paste0("iter",i, " errorY= ",erry<-round2(mean((Y-Z%*%B)^2)), ", Bdiff= ",round2(Bdiff), ", Bkappa=", round2(kappa(B))))
    }

    #check for convergence
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
    }else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }#else if

    if(Bdiff<tol &&i>40){
      message(paste0("converged at  iteration ", i))
      break
    }#if

    if( BdiffCount>5&&i>40){
      message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
      break
    }#if

  }#for i

  rownames(B)=colnames(Z)=paste("LV",1:k)
  Zproject=Z%*%solve(crossprod(Z)+L2*diag(k))
  return(list(B=B, Z=Z, Zraw=Zraw, Zproject=Zproject,L1=L1, L2=L2, L4 = L4, k = k))
}#manifoldDecomp






#' Variance explained by latent variables
#' @export
VarianceExplained <- function(Y, Z, B, option = "simple", normalize = F){
  if (option == "simpleboth"){
    if (normalize){
      Z.norm <- apply(Z,2,normF)
      B.norm <- apply(B, 1, normF)
      Z <- sweep(Z,2,Z.norm, "/")
      B <- sweep(B,1,B.norm, "/")
    }#if

    res <- diag(t(Z)%*%Y%*%t(B))
  }else if (option == "simple"){
    res <- diag(B%*%t(Y)%*%Y%*%t(B))

  }else if (option == "firstfew"){
    res <- c()
    for (ii in 1:nrow(B)){
      print(ii)
      Vk <- t(B[1:ii, ,drop=F])
      Xk <- Y%*%Vk%*% solve(t(Vk) %*% Vk)%*%t(Vk)
      #?Y%*%Vk = Zk
      res[ii] <- sum(diag(t(Xk) %*% Xk))
    }#for ii
    res <- res-c(0, res[1:(length(res)-1)])

  }else if (option == "firstfewboth"){
    res <- c()
    for (ii in 1:nrow(B)){
      print(ii)
      Vk <- t(B[1:ii, ,drop=F])
      Zk <- Z[,1:ii, drop =F]
      Xk <- Zk%*% solve(t(Vk) %*% Vk)%*%t(Vk)
      #?Y%*%Vk = Zk
      res[ii] <- sum(diag(t(Xk) %*% Xk))
    }#for ii
    res <- res-c(0, res[1:(length(res)-1)])
  }#else if

  names(res) <- rownames(B)
  return(res)
}#VarianceExplained




recon.error <- function(Y, Y.re, index){
  temp <- (Y-Y.re)^2
  return(mean(temp[index]))
}#recon.error

