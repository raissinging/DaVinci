# https://github.com/jmonlong/Hippocamplus/blob/master/content/post/2018-06-09-ClusterEqualSize.Rmd
#' Same Size Clustering
#'
#' This is a wrapper for several implementation that classify samples into
#' same size clusters, the details please see [this blog](http://jmonlong.github.io/Hippocamplus/2018/06/09/cluster-same-size/).
#' The source code is modified based on code from the blog.
#'
#' @param mat a data/distance matrix.
#' @param diss if `TRUE`, treat `mat` as a distance matrix.
#' @param clsize integer, number of sample within a cluster.
#' @param algo algorithm.
#' @param method method.
#'
#' @return a vector.
#' @export
#'
#' @examples
#' set.seed(1234L)
#' x <- rbind(
#'   matrix(rnorm(100, sd = 0.3), ncol = 2),
#'   matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2)
#' )
#' colnames(x) <- c("x", "y")
#'
#' y1 <- same_size_clustering(x, clsize = 10)
#' y11 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, diss = TRUE)
#'
#' y2 <- same_size_clustering(x, clsize = 10, algo = "hcbottom", method = "ward.D")
#'
#' y3 <- same_size_clustering(x, clsize = 10, algo = "kmvar")
#' y33 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, algo = "kmvar", diss = TRUE)
same_size_clustering <- function(mat, diss = FALSE, clsize = NULL,
                                 algo = c("nnit", "hcbottom", "kmvar"),
                                 method = c(
                                   "maxd", "random", "mind", "elki",
                                   "ward.D", "average", "complete", "single"
                                 )) {
  stopifnot(is.numeric(clsize))

  algo <- match.arg(algo)
  method <- match.arg(method)
  do.call(algo, args = list(mat = mat, diss = diss, clsize = clsize, method = method))
}#same_size_clustering




#minimized the difference in set size
minimal.diff <- function(n, base_size){
  k <- n %/% base_size
  base_size <- n %/% k
  remainder <- n %%k
  set_sizes <- c(rep(base_size, k-remainder), rep(base_size+1, remainder))
  return(set_sizes)
}#minimal.diff



nnit <- function(mat,
                 clsize = NULL,
                 diss = FALSE,
                 method = "maxd") {
  
  stopifnot(is.logical(diss))

  #default, can go below clsize
  #clsize.rle <- rle( as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))) )
  #clsize <- clsize.rle$lengths
  
  #make sure the cluster size is at least clsize
  #clsize.rle <- rle( as.numeric(cut(1:nrow(mat), floor( nrow(mat) / clsize))) )
  #clsize <- clsize.rle$lengths
  
  #equvalant, but easier to understand
  clsize <- minimal.diff(nrow(mat), clsize)
  
  lab <- rep(NA, nrow(mat))
  if (isFALSE(diss)) {
    dmat <- as.matrix(dist(mat))
  } else {
    dmat <- mat
  }
  cpt <- 1
  while (sum(is.na(lab)) > 0) {
    lab.ii <- which(is.na(lab))
    dmat.m <- dmat[lab.ii, lab.ii, drop = F]
    ii <- switch(method,
                 maxd = which.max(rowSums(dmat.m)),
                 mind = which.min(rowSums(dmat.m)),
                 random = sample.int(nrow(dmat.m), 1),
                 stop("unsupported method in 'nnit'!")
    )
    lab.m <- rep(NA, length(lab.ii))
    lab.m[head(order(dmat.m[ii, ]), clsize[cpt])] <- cpt
    lab[lab.ii] <- lab.m
    cpt <- cpt + 1
  }
  if (any(is.na(lab))) {
    lab[which(is.na(lab))] <- cpt
  }
  lab
}#nnit


kmvar <- function(mat,
                  clsize = NULL,
                  diss = FALSE,
                  method = "maxd") {
  stopifnot(is.logical(diss))

  k <- ceiling(nrow(mat) / clsize)
  if (isFALSE(diss)) {
    km.o <- kmeans(mat, k)
    # distance to centers
    centd <- lapply(1:k, function(kk) {
      euc <- t(mat) - km.o$centers[kk, ]
      sqrt(apply(euc, 2, function(x) sum(x^2)))
    })
    centd <- matrix(unlist(centd), ncol = k)
  } else {
    message("PAM algorithm is applied when input distance matrix.")
    pam.o <- cluster::pam(mat, k, diss = TRUE)
    # medoids
    # distance to medoids
    centd <- mat[, pam.o$id.med, drop = FALSE]
  }

  labs <- rep(NA, nrow(mat))
  clsizes <- rep(0, k)

  ptord <- switch(method,
                  maxd = order(-apply(centd, 1, max)),
                  mind = order(apply(centd, 1, min)),
                  random = sample.int(nrow(mat)),
                  elki = order(apply(centd, 1, min) - apply(centd, 1, max)),
                  stop("unsupported method in 'kmvar'!")
  )

  for (ii in ptord) {
    bestcl <- which.max(centd[ii, ])
    labs[ii] <- bestcl
    clsizes[bestcl] <- clsizes[bestcl] + 1
    if (clsizes[bestcl] >= clsize) {
      centd[, bestcl] <- NA
    }
  }
  return(labs)
}#kmvar


#' @importFrom stats as.dist
hcbottom <- function(mat,
                     clsize = NULL,
                     diss = FALSE,
                     method = "ward.D") {
  stopifnot(is.logical(diss))

  method <- match.arg(method, choices = c("ward.D", "average", "complete", "single"))
  if (isFALSE(diss)) {
    dmat <- as.matrix(dist(mat))
  } else {
    dmat <- mat
  }
  clsize.rle <- rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))))
  clsizes <- clsize.rle$lengths
  cpt <- 1
  lab <- rep(NA, nrow(mat))
  for (clss in clsizes[-1]) {
    lab.ii <- which(is.na(lab))
    hc.o <- hclust(as.dist(dmat[lab.ii, lab.ii]), method = method)
    clt <- 0
    ct <- length(lab.ii) - clss
    while (max(clt) < clss) {
      cls <- cutree(hc.o, ct)
      clt <- table(cls)
      ct <- ct - 1
    }
    cl.sel <- which(cls == as.numeric(names(clt)[which.max(clt)]))
    lab[lab.ii[head(cl.sel, clss)]] <- cpt
    cpt <- cpt + 1
  }
  lab[is.na(lab)] <- cpt
  lab
}#hcbottom








#' Tile the slice for scalable implementation
#' @import sf
#' @export
tile_the_slice <- function(coord, 
                           random.seed = 1, 
                           L2_number = 1000, 
                           tile.minimum = 3){
  
  if (is.null(rownames(coord))){
    stop("Input doesn't have rownames.")
  }#if

  #pre-set
  L1_size <- 7000
  #L2_number <- 1000 #about 1000 spots in the problem
  message("L2_number is not recommended to be higher than 2000. Higher the number is, the more granularity you expect.")

  #derived
  L1_number <- nrow(coord)/L1_size #number
  L2_size <- nrow(coord) %/% L2_number #balanced.cluster.size

  if (L2_size < tile.minimum){
    L2_size <- tile.minimum
  }#if 

  message(paste0("Each tile includes ", L2_size, " spots."))


  #tile the slice
  #################################
  bbox <- sf::st_bbox(c(xmin = min(coord[,1]), 
                        ymin = min(coord[,2]), 
                        xmax = max(coord[,1]), 
                        ymax = max(coord[,2]) ))
  
  xy_ratio <- (bbox$xmax-bbox$xmin)/(bbox$ymax-bbox$ymin)
  ny <- ceiling(sqrt(L1_number/xy_ratio))
  nx <- ceiling(ny*xy_ratio)
  grid <- sf::st_make_grid(bbox, n = c(nx, ny), square = T)
  pos <- sf::st_as_sf(as.data.frame(coord), coords = c("array_row","array_col"))


  #adjust the tile: size requirement
  #################################
  point_in_grid <- sf::st_intersects(grid, pos)
  point_counts <- lengths(point_in_grid)
  exclude_index <- sort(which(point_counts < tile.minimum), decreasing = T)
  

  while (length(exclude_index)>0){
    message("Tile is being adjusted to include enough points.")
    
    for (cell in exclude_index){
      
      if (point_counts[cell]==0){
        grid <- grid[-cell]
      }else{
        
        adjacency_matrix <- sf::st_touches(grid)
        adjacent_cells <- adjacency_matrix[[cell]]
        
        #merge with the first cell
        target_cell <- setdiff( adjacent_cells, exclude_index)[1] 
        
        if (is.na(target_cell)){
          
          message("No cell survives. Need to consider 2nd level neighbors.")
          
          new_adjacent_cells <- c()
          for (ii in 1:length(adjacent_cells)){
            new_adjacent_cells <- c(new_adjacent_cells, adjacency_matrix[[adjacent_cells[ii]]])
          }#for ii
          
          target_cell <- setdiff(unique(new_adjacent_cells), exclude_index)[1]
        }#if 
        
        
        #merged step
        merged_polygon <- sf::st_union(grid[c(cell, target_cell)])
        grid[target_cell] <- merged_polygon
        grid <- grid[-cell]
        
      }#else
      
    }#for cell

    point_in_grid <- sf::st_intersects(grid, pos)
    point_counts <- lengths(point_in_grid)
    exclude_index <- sort(which(point_counts < tile.minimum), decreasing = T)
    
  }#while length(exclude_index)>0
  
  #validation
  #table(point_counts)
  
  #adjust the tile: avoid colinearity
  #################################
  collinearity.index <- c()

  for (ii in 1:length(point_in_grid)){

    num.x <- length(unique(coord[point_in_grid[[ii]],1]))
    num.y <- length(unique(coord[point_in_grid[[ii]],2]))
    if ( (num.x==1) | (num.y==1) ){
        collinearity.index <- c(collinearity.index, ii)
    }#if 
  }#for ii


  while ( length(collinearity.index)>0 ){
    
    message("Tile is being adjusted to remove collinearity.")

    for (cell in collinearity.index){

        adjacency_matrix <- sf::st_touches(grid)
        adjacent_cells <- adjacency_matrix[[cell]]
        
        #merge with the first cell
        target_cell <- setdiff( adjacent_cells, collinearity.index)[1] 
        
        if (is.na(target_cell)){
          
          message("No cell survives. Need to consider 2nd level neighbors.")
          
          new_adjacent_cells <- c()
          for (ii in 1:length(adjacent_cells)){
            new_adjacent_cells <- c(new_adjacent_cells, adjacency_matrix[[adjacent_cells[ii]]])
          }#for ii
          
          target_cell <- setdiff(unique(new_adjacent_cells), collinearity.index)[1]
        }#if 
                
        #merged step
        merged_polygon <- sf::st_union(grid[c(cell, target_cell)])
        grid[target_cell] <- merged_polygon
        grid <- grid[-cell]
        
    }#for cell



    #recalcualte the colinearity.idnex
    point_in_grid <- sf::st_intersects(grid, pos)
    collinearity.index <- c()
    for (ii in 1:length(point_in_grid)){

      num.x <- length(unique(coord[point_in_grid[[ii]],1]))
      num.y <- length(unique(coord[point_in_grid[[ii]],2]))
      if ( (num.x==1) | (num.y==1) ){
          collinearity.index <- c(collinearity.index, ii)
      }#if 
    }#for ii

  }#while length(collinearity.index) > 0




  #calculate the L1 label
  #################################
  ids_L1 <- unlist(lapply(sf::st_intersects(pos, grid), function(sublist) sublist[1]))
  names(ids_L1) <- rownames(coord)


  #L2 fine grid
  #################################
  ids_L2 <- ids_L1
  ids_L1_unique <- sort(unique(ids_L1))

  for (ii in 1:length(ids_L1_unique)){
    temp <- which(ids_L1==ids_L1_unique[ii])
    if (length(temp)<L2_size){
      ids_L2[temp] <- 1
    }else{
      ids_L2[temp] <- swk(coord[temp,], 
                          method = "balanced", 
                          L2norm = F, 
                          balanced.cluster.size = L2_size, 
                          random.seed = random.seed, 
                          verbose = F)
    }#else

  }#for ii

  #combine together
  #################################
  grid_ids <- paste0(ids_L1, "_", ids_L2)
  mapping <- 1:length(unique(grid_ids))
  names(mapping) <- unique(grid_ids)
  grid_ids <- mapping[grid_ids]

  names(grid_ids) <- rownames(coord)
  return(grid_ids)
}#tile_the_slice







#' Scalable implementation
#' 
#' @param L2_number Number of meta spot expected. L2_number is not recommended to be higher than 3000 if you want to get results really fast. 
#' 
#' @export
manifoldDecomp.scalable <- function(gene.exp, coor, L2_number = 4000, k.arg = 20, L4.arg = 50, max.iter = 10, frac.thr = 0.95){
  
  if (!all(colnames(gene.exp)==rownames(coor))){
      stop("gene.exp and coor don't align.")
  }#if

  #step 0: filter spot without normalization
  #############################################################################
  use_gene <- which(Matrix::rowSums(gene.exp==0)/ncol(gene.exp) < frac.thr)
  if (is.data.frame(gene.exp)){
    gene.exp <- as.matrix(gene.exp[use_gene,])
  }else{
    gene.exp <- gene.exp[use_gene,]
  }#else

  #check standard deviation
  gene.exp <- gene.exp[which(apply(gene.exp,1,sd)>0),]

  #remove all zero cells
  gene.exp <- gene.exp[, which( Matrix::colSums(gene.exp)>0)]
  
  #remember to filter coor accordingly
  coor <- coor[colnames(gene.exp),]



  #step 1: generate a partition
  #############################################################################
  if (nrow(coor) < L2_number*k.arg){
    L2_number <- nrow(coor) %/% k.arg
    message( paste0("L2_number is reset to ", L2_number, " to fulfill the rank requirement."))
  }#if

  grid_ids <- tile_the_slice(coor, 
                             random.seed = 1, 
                             L2_number = L2_number, 
                             tile.minimum = k.arg+1)
  message("Tiling finishes.")
  
  #mean aggregation per tile
  ##################################
  unique_grid_ids <- unique(grid_ids)
  tile_list <- lapply(unique_grid_ids, function(x){which(grid_ids==x)})
  tile_aggregation <- lapply(tile_list, function(x){apply(gene.exp[,names(x),drop=F],1,mean)} )
  tile.mat <- do.call(cbind, tile_aggregation)
  colnames(tile.mat) <- paste0("Tile_", unique_grid_ids)
  names(tile_list) <- paste0("Tile_", unique_grid_ids)
    
  
  #keep the coordinates for each tile
  ##################################
  coor_aggregation <- lapply(tile_list, function(x){apply(coor[names(x),,drop=F],2,mean)} )
  coor.tile <- do.call(rbind, coor_aggregation)
  colnames(coor.tile) <- c("array_col", "array_row")
  rownames(coor.tile) <- paste0("Tile_",seq_along(coor.tile[,1]))
  message("Aggregation per tile finishes.")

  
  #step 2: tile-level decomposition, get Z
  ##################################
  #tile.mat, coor.tile
  fin.tile <- DaVinci::preprocess(tile.mat, coor.tile, type = "rna", graph.opt = "Tri.mesh",  frac.thr = 0.95)
  
  gene.exp.tile <- fin.tile$mat
  L.tile <- fin.tile$L

  #modify gene.exp and coor accordingly
  gene.exp <- gene.exp[rownames(gene.exp.tile),]


  message("Tile-level decomposition starts.")  
  ICAp.res.tile0 <- manifoldDecomp_adaptive(gene.exp.tile, L.tile, k = k.arg, L4 = L4.arg, L4_adaptive = 2, to_drop = T, save.complete = T)
  
  
  #step 3: within-tile update
  ##################################
  gene.exp.list <- list()
  L.list <- list()
  L1.grid <- list()
  L2.grid <- list()
  shur0.grid <- list()
  density.grid <- unlist(lapply(tile_list, length))
  B.list <- list()
  
  for (ii in 1:length(tile_list)){
    #print(ii)
    if (density.grid[ii] > 3){
      gene.exp.inuse <- as.matrix(gene.exp[,names(tile_list[[ii]])])
      
      coor.inuse <- coor[names(tile_list[[ii]]), ,drop = F]
      temp <- L_generate(coor.inuse, opt = "Tri.mesh")
      
      gene.exp.list[[ii]] <- gene.exp.inuse
      L.list[[ii]] <- temp$L
      
      ICAp.res.inuse <- impact_adaptive(ICAp.res.tile0$Z, gene.exp.inuse, query.L = temp$L, query.L4 = L4.arg, to_drop = F, scale=1, max.iter = 200, cor.thr = 0.8, verbose = F)
      
      #Slide.LvPlot(coor.inuse, LVs= ICAp.res.inuse$B, gene.verbose = F, plot.all = T)
      
      L1.grid[[ii]] <- ICAp.res.inuse$L1
      L2.grid[[ii]] <- ICAp.res.inuse$L2
      shur0.grid[[ii]] <- ICAp.res.inuse$shur0
      B.list[[ii]] <- ICAp.res.inuse$B
      
    }else{
      
      gene.exp.list[[ii]] <- NA
      L.list[[ii]] <- NA
      L1.grid[[ii]] <- NA  
      L2.grid[[ii]] <- NA
      shur0.grid[[ii]] <- NA
      B.list[[ii]] <- NA
    }#else
    
  }#for ii
  
  
  #step 4: construct new B to update Z
  ##################################
  #global B - averaging B from tile
  B.agg.pre <- lapply(B.list, function(x){apply(x,1,mean)} )
  B.aggregate <- do.call(cbind, B.agg.pre)
  
  
  
  
  #repeat step 2,3,4 till convergence
  #############################################################################
  error.relative <- normF(ICAp.res.tile0$B - B.aggregate)/normF(ICAp.res.tile0$B)
  error.accum <- c(error.relative)
  
  count <- 0
  message("Iteration starts.")
  
  while ( (error.relative >=1e-3) & (count <= max.iter) ){
    
    count <- count+1
    print(count)
    
    #check correlation
    #B.to_report <- do.call(cbind, B.list)
    #corr.list[[count]] <- cor(t(B.to_report), t(gt.B[,colnames(B.to_report)]))
    
    
    #ptm <- proc.time()
    ICAp.res.tile <- manifoldDecomp_adaptive(gene.exp.tile, L.tile, B = B.aggregate, k = k.arg,  svdres = NA, L1 = ICAp.res.tile0$L1, L2 = ICAp.res.tile0$L2, L4 = L4.arg, L4_adaptive = 2, to_drop = T, save.complete = T)
    #print(proc.time()-ptm)
    
    #visualization check
    #pheatmap::pheatmap(cor(t(ICAp.res.tile$B)), display_numbers = T, fontsize = 20)
    #Slide.LvPlot(coor.tile, LVs= ICAp.res.tile$B, gene.verbose = F, plot.all = T)
    
    
    
    #tile-wise update
    ###############################################################################
    B.list <- list()
    
    
    for (ii in 1:length(tile_list)){
      #print(ii)
      if (density.grid[ii] > 3){
        
        ICAp.res.inuse <- impact_adaptive(ICAp.res.tile$Z, gene.exp.list[[ii]] , query.L = L.list[[ii]], query.L1 = L1.grid[[ii]], query.L2 = L2.grid[[ii]], query.L4 = L4.arg, query.shur0 = shur0.grid[[ii]], to_drop = F, scale=1, max.iter = 200, cor.thr = 0.8, verbose =F)
        
        #Slide.LvPlot(coor.inuse, LVs= ICAp.res.inuse$B, gene.verbose = F, plot.all = T)
        B.list[[ii]] <- ICAp.res.inuse$B
        
      }else{        
        B.list[[ii]] <- NA
      }#else
    }#for ii
    
    
    #global B - averaging B from tile
    B.agg.pre <- lapply(B.list, function(x){apply(x,1,mean)} )
    B.aggregate <- do.call(cbind, B.agg.pre)
    
    #visualization check
    #Slide.LvPlot(coor.tile, LVs= B.aggregate, gene.verbose = F, plot.all = T)
    
    
    error.relative <- normF(ICAp.res.tile$B - B.aggregate)/normF(ICAp.res.tile$B)
    message(paste0("Error:", error.relative))
    error.accum <- c(error.accum, error.relative)
    
  }#while
  
  
  #return the result
  #############################################################################
  
  #tile level: coor.tile, gene.exp.tile, B.aggregate, ICAp.res.tile
  
  #spot level: B.spot, B.list
  B.spot <- do.call(cbind, B.list)[,rownames(coor)]
  
  
  colnames(B.aggregate) <- rownames(coor.tile)

  return(list(coor.tile = coor.tile, gene.exp.tile = gene.exp.tile, B.aggregate = B.aggregate, dav.res.tile = ICAp.res.tile, B.list = B.list, B.spot = B.spot, error = error.accum))
  
  
}#manifoldDecomp.scalable


