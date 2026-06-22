
#' @details This function requires `Seurat`. Make sure it is installed.
GenePlot <- function(dat, features, pt.size = 2, ratio = 1, use.myratio = F){
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for this function but is not installed. Please install it.")
  }
  coord <- Seurat::GetTissueCoordinates(object = dat@images[[1]])
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  message(paste0("The aspect ratio is ", myratio))
  
  if (use.myratio){
    ratio <- myratio
  }#if
  
  SpatialFeaturePlot(dat, features = features, pt.size.factor = pt.size)+ggplot2::theme(aspect.ratio = ratio) 
}#GenePlot




#input should be the gene, L and ICAp.res
reciprocal_default <- function(Y.list, Lap, ICAp.res.list, num.of.reference.slice = NULL){
  
  num.of.slice <- length(Y.list)
  if (is.null(num.of.reference.slice)){
    num.of.reference.slice <- length(ICAp.res.list)  
  }#if
  #else it will be first several slices based on the provided numbers
  message(paste0("The number of refernce slices is: " , num.of.reference.slice))
  message(paste0("The number of total slices is: ", num.of.slice))

  #calculate the projection and construct a 2d list
  ##############################################################################
  proj <- list()
  
  for (ii in 1:num.of.reference.slice){
    proj.ii.on.jj <- list()
    
    for (jj in 1:num.of.slice){
    
      if (ii == jj){
        proj.ii.on.jj[[jj]] <- ICAp.res.list[[jj]]
        
      }else{
          proj.ii.on.jj[[jj]] <- impact(ICAp.res.list[[ii]]$Z, Y.list[[jj]], Lap[[jj]], ICAp.res.list[[jj]]$L4, query.shur0 = ICAp.res.list[[jj]]$shur0, query.ICAp.res = ICAp.res.list[[jj]], scale= 1, max.iter = 200, cor.thr = 0.8)
          
      }#else
      
    }#for jj 
    
    proj[[ii]] <- proj.ii.on.jj
  }#for ii
  
  return(proj)
  
}#reciprocal_default





p.intersect <- function(x){
  
  res <- x[[1]]
  if (length(x)>1){
    
    for (ii in 2:length(x)){
      res <- intersect(res, x[[ii]])
    }#for ii
  }#if
  
  return(res)
}#p.intersect



#' apply the existing integration result on the new dataset
#' @export
reciprocal_with_Z <- function(Y.list, Lap, ICAp.res.list, ICAp.res.list.reference, Z.names){
  
  num.of.slice <- length(Y.list)
  num.of.reference.slice <- length(ICAp.res.list.reference)
  message(paste0("The number of refernce slices is: " , num.of.reference.slice))
  message(paste0("The number of total slices is: ", num.of.slice))
  
  
  #first construct Z
  ##########################################################
  slice.index <- as.numeric(lapply(strsplit(Z.names, "_"), function(x){x[1]}))
  LV.names <- as.character(lapply(strsplit(Z.names, "_"), function(x){x[2]}))
  
  Z.inuse.list <- list()
  for (ii in 1:length(Z.names)){
      Z.inuse.list[[ii]] <- ICAp.res.list.reference[[slice.index[ii]]]$Z[,LV.names[ii]]
  }#for ii
  
  gene.names <- p.intersect(lapply(Z.inuse.list, names))
    
  Z.inuse <- do.call(cbind, lapply(Z.inuse.list, function(x){x[gene.names]}))
  rownames(Z.inuse) <- gene.names
  colnames(Z.inuse) <- Z.names
  
  #calculate the projection
  ##########################################################
  proj <- list()
  
  for (ii in 1:num.of.slice){
    proj[[ii]] <- impact(Z.inuse, Y.list[[ii]], Lap[[ii]], ICAp.res.list[[ii]]$L4, query.shur0 = ICAp.res.list[[ii]]$shur0, query.ICAp.res = ICAp.res.list[[ii]], scale= 1, max.iter = 200, cor.thr = 0.8)
  }#for ii
  
  LVs.inuse <- do.call(cbind, lapply(proj, function(x){x$B}))
  
  col.names <- c()
  for (ii in 1:num.of.slice){
    col.names <- c(col.names, paste0(ii, "_", colnames(Y.list[[ii]])) )
  }#for ii
  
  rownames(LVs.inuse) <- Z.names
  colnames(LVs.inuse) <- col.names
  
  return(LVs.inuse)
}#reciprocal_with_Z





reciprocal_concat <- function(proj, option = "L2norm"){

  num.of.slice <- length(proj)
  
  #row names and colnames
  row.names <- c()
  col.names <- c()
  for (ii in 1:num.of.slice){
    row.names <- c(row.names, paste0(ii, "_", rownames(proj[[ii]][[ii]]$B)) )
    col.names <- c(col.names, paste0(ii, "_", colnames(proj[[ii]][[ii]]$B)) )
  }#for ii
   
  
  #merge
  dim.list <- list()
  
  for (ii in 1:num.of.slice){
    if (option == "default"){
      dim.list[[ii]] <- do.call(cbind, lapply(proj[[ii]], function(x){x$B}) )
    }else if (option == "L2norm"){
      dim.list[[ii]] <- do.call(cbind, lapply( lapply(proj[[ii]], function(x){x$B}), L2Norm))
    }else if (option == "L2norm.joint"){
      dim.list[[ii]] <- L2Norm(do.call(cbind, lapply(proj[[ii]], function(x){x$B}) ))
    }#else if
  }#for ii
  
  #construct the LVs
  ###########################################
  LVs <- do.call(rbind, dim.list)

  rownames(LVs) <- row.names
  colnames(LVs) <- col.names
  
  return(LVs)
}#reciprocal_concat



#construct the B.list and global LVs with different options
#' Finalize the integration by filtering out duplicated LVs based on similarity and strategy
#' @export
reciprocal_deco <- function(LVs, 
                            LVs.filter.thr = 0.9, 
                            mod = "all"){
    
  #filter by correlation
  ###########################################
   
    
    if (mod == "all"){
      
      #use as many LVs < LVs.filter.thr as possible to generate the final LVs in use
      cor.res <- cor(t(LVs))
      
      diag(cor.res) <- 0
      cor.res[lower.tri((cor.res))] <- 0

      cor.to_drop <- which(cor.res > LVs.filter.thr, arr.ind = T)
      cor.to_drop.cor <- which(cor.res > LVs.filter.thr, arr.ind = F)
      cor.to_drop.cor <- cor.res[cor.to_drop.cor]
      
      cor.to_drop <- as.data.frame(cor.to_drop)
      cor.to_drop$cor <- cor.to_drop.cor
      cor.to_drop <- cor.to_drop[order(cor.to_drop$cor, decreasing = T),]

      LV.to_drop <- c()
      while(nrow(cor.to_drop)>0){
          LV.index <- max(cor.to_drop[1,1:2])
          LV.to_drop <- c(LV.to_drop, LV.index)
          tmp <- which(cor.to_drop[,1]!=LV.index & cor.to_drop[,2]!=LV.index)
          cor.to_drop <- cor.to_drop[tmp,]
      }#while

      LVs <- LVs[-LV.to_drop,]


      #not efficient
      if (F){
          #iteratively remove the correlations
          while (max(cor.res)> LVs.filter.thr){
        #remove the larger LV index
        LV.to_drop <- max(which(cor.res==max(cor.res[upper.tri(cor.res)]), arr.ind = T))
        LVs <- LVs[-LV.to_drop,]
        
        cor.res <- cor.res[-LV.to_drop, -LV.to_drop]
        #cor.res <- cor(t(LVs))
        #diag(cor.res) <- 0
        #cor.res[lower.tri((cor.res))] <- 0
      }#while 
   
      }#if F
      
         
    }else if (mod == "common"){
      
      #only selected the most highly correlated ones
      cor.res <- cor(t(LVs))
      
      #converted to a connected graph based on the threshold
      adj_matrix <- cor.res > LVs.filter.thr
      adj_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = F)
      components_list <- igraph::components(adj_graph)
      components_selected <- which(components_list$csize>1)
      
      #only pick up one from correlated blocks
      include.names <- c()
      for (ii in 1:length(components_selected)){
        #select the one with smaller index
        include.names <- c(include.names, names(components_list$membership)[which(components_list$membership==components_selected[ii])][1] )
      }#for ii
      
      LVs <- LVs[include.names,]
      
    }#else if mod == "common"
    
  return(LVs)
}#reciprocal_deco










#' Horizontal Integration
#'
#' @param LVs.filter.thr The correlation threshold for LV filtering.
#' @param mod The deafult value is "all" which takes all LVs to integrate, and the alternative value is "common", which only includes highly correlated LVs to integrate.
#' 
#' @import harmony
#' @export
Horizontal.Integration <- function(Y.list, 
                                   L.list, 
                                   coor.list = NULL, #can be removed
                                   dav.res.list = NULL, 
                                   k.args = rep(12, length(Y.list)), 
                                   LVs.filter.thr = 0.8, 
                                   mod = "all", 
                                   remove.LV1 = F, 
                                   L2.option = "default", 
                                   LV.filter = T,
                                   batch.correction = "harmony" #NULL, harmony, combat
                                   ){
  
  if ((length(Y.list) > 30) & (mod == "all")){
    message("There are more than 30 slices. 'mod' is recommended to be set as 'common'")
  }#if
  
  #single slice
  ##########################################################
  if (is.null(dav.res.list)){
    dav.res.list <- list()
  
    for (ii in 1:length(Y.list)){
      dav.res.list[[ii]] <- manifoldDecomp_adaptive(Y.list[[ii]], L.list[[ii]], k = k.args[ii], L4 = 50, L4_adaptive = 2, to_drop = T, save.complete = T)
    }#for ii
    
  }#if
  
  #horizontal integration  
  ##########################################################
  message("Reciprocal projection starts.")
  
  proj <- reciprocal_default(Y.list, L.list, dav.res.list)
  message("Projection finishes.")

  LVs <- reciprocal_concat(proj, option = L2.option)
  message("Concatenation finishes.")


  message("Select the principal axis.")
  if (LV.filter){
    LVs <- reciprocal_deco(LVs, LVs.filter.thr = LVs.filter.thr, mod = mod)
  }#if LV.filter
  message("Selection finishes.")
  

  #remove LV 1
  #if (remove.LV1){
  #  LVs <- LVs[which(!grepl("LV 1$", rownames(LVs))),]  
  #}#if remove.L1
  
  if (remove.LV1){
      LVs <- LVs[which(!rownames(LVs) %in% remove.LV1),]
  }#if 
  
  
  #batch correction
  ##########################################################
   if (is.null(batch.correction)){

    LVs_embeddings <- t(LVs)

  }else if (batch.correction=="harmony"){

    message("Batch correction starts.")
    .meta <- data.frame(section = unlist(lapply(strsplit(colnames(LVs), "_"), function(x){x[[1]]})))
    harmony_object <- harmony::RunHarmony(data_mat = t(LVs), meta_data = .meta,
                                          vars_use = "section", return_object = TRUE)
    LVs_embeddings <- t(harmony_object$getZcorr())
    rownames(LVs_embeddings) <- colnames(LVs)
    colnames(LVs_embeddings) <- rownames(LVs)
    message("Batch correction finishes.")

  }else if (batch.correction == "combat"){
    message("Batch correction starts.")
    LVs_embeddings <- sva::ComBat(LVs,
                                  unlist(lapply(strsplit(colnames(LVs), "_"), function(x){x[[1]]}))
                                  )
    LVs_embeddings <- t(LVs_embeddings)
    rownames(LVs_embeddings) <- colnames(LVs)
    message("Batch correction finishes.")

  }#else if
  

  slice_id <- unlist(lapply(strsplit(rownames(LVs_embeddings), "_"), function(x){x[1]}))
  
  #extract the loadings corrresponding to LVs_embeddings
   ######################################################
   loading.uid <- colnames(LVs_embeddings)
   loading.ii <- as.numeric(unlist(lapply(strsplit(colnames(LVs_embeddings), "_"), function(x){x[1]})))

  if (stringr::str_count(rownames(LVs)[1], "_") ==1){
    loading.jj <- unlist(lapply(strsplit(colnames(LVs_embeddings), "_"), function(x){x[2]}))
  }else{
    loading.jj <- unlist(lapply(strsplit(colnames(LVs_embeddings), "_"), function(x){ paste0(x[-1], collapse = "_")}))
  }#else

   
   #extract from dav.res.list
   ##################################
   loading <- list()
   for (ii in 1:length(loading.uid)){
      loading[[ii]] <- dav.res.list[[loading.ii[ii]]]$Z[,loading.jj[ii]]

   }#for ii

   names(loading) <- loading.uid 

  return(list(LVs = LVs, 
              LVs_embeddings = LVs_embeddings, 
              slice_id = slice_id, 
              loading = loading))
  
  #LVs - the intermediate resutls after reciprocal_concat
  #LVs_embeddings - the product
  #slice_id - sample id per slice
  
}#Horizontal.Integration





# harmony first, then filter LVs based on the correlation
# instead of filtering LVs first and then running harmony
Horizontal.Integration.first <- function(Y.list, 
                                         L.list, 
                                         coor.list = NULL, #can be removed
                                         dav.res.list = NULL, 
                                         k.args = rep(12, length(Y.list)), 
                                         LVs.filter.thr = 0.8, 
                                         mod = "all", 
                                         remove.LV1 = F, 
                                         L2.option = "default", 
                                         LV.filter = T,
                                         batch.correction = "harmony" #NULL, harmony, combat
                                         ){
  
  if ((length(Y.list) > 30) & (mod == "all")){
    message("There are more than 30 slices. 'mod' is recommended to be set as 'common'")
  }#if
  
  #single slice
  ##########################################################
  if (is.null(dav.res.list)){
    dav.res.list <- list()
  
    for (ii in 1:length(Y.list)){
      dav.res.list[[ii]] <- manifoldDecomp_adaptive(Y.list[[ii]], L.list[[ii]], k = k.args[ii], L4 = 50, L4_adaptive = 2, to_drop = T, save.complete = T, verbose = F)
    }#for ii
    
  }#if
  
  #horizontal integration  
  ##########################################################
  message("Reciprocal projection starts.")
  
  #ptm <- proc.time()
  proj <- reciprocal_default(Y.list, L.list, dav.res.list)
  #print(proc.time()-ptm)
  message("Projection finishes.")
  
  LVs <- reciprocal_concat(proj, option = L2.option)
  message("Concatenation finishes.")

  #batch correction
  ##########################################################
  if (is.null(batch.correction)){

    LVs_embeddings <- t(LVs)

  }else if (batch.correction == "harmony"){
    message("Batch correction starts.")
    .meta <- data.frame(section = unlist(lapply(strsplit(colnames(LVs), "_"), function(x){x[[1]]})))
    harmony_object <- harmony::RunHarmony(data_mat = t(LVs), meta_data = .meta,
                                          vars_use = "section", return_object = TRUE)
    LVs_embeddings <- t(harmony_object$getZcorr())
    rownames(LVs_embeddings) <- colnames(LVs)
    colnames(LVs_embeddings) <- rownames(LVs)
    message("Batch correction finishes.")

  }else if (batch.correction == "combat"){
    message("Batch correction starts.")
    LVs_embeddings <- sva::ComBat(LVs,
                                  unlist(lapply(strsplit(colnames(LVs), "_"), function(x){x[[1]]}))
                                  )
    LVs_embeddings <- t(LVs_embeddings)
    rownames(LVs_embeddings) <- colnames(LVs)
    message("Batch correction finishes.")

  }#else if
  

  slice_id <- unlist(lapply(strsplit(rownames(LVs_embeddings), "_"), function(x){x[1]}))



  #Deco
  ##########################################################
  message("Select the principal axis.")
  
  if (LV.filter){
    #ptm <- proc.time()
    LVs_embeddings <- reciprocal_deco(t(LVs_embeddings), 
                                      LVs.filter.thr = LVs.filter.thr, 
                                      mod = mod)
    #print(proc.time()-ptm)
    LVs_embeddings <- t(LVs_embeddings)
  }#if LV.filter

  message("Selection finishes.")
  
    
  #remove LV 1
  #if (remove.LV1){
  #  LVs_embeddings <- LVs_embeddings[,which(!grepl("LV 1$", colnames(LVs_embeddings)))]  
  #}#if remove.L1
  
  if (remove.LV1){
      LVs_embeddings <- LVs_embeddings[,which(!colnames(LVs_embeddings) %in% remove.LV1)]
  }#if   
   
   #extract the loadings corrresponding to LVs_embeddings
   ######################################################
   loading.uid <- colnames(LVs_embeddings)
   loading.ii <- as.numeric(unlist(lapply(strsplit(colnames(LVs_embeddings), "_"), function(x){x[1]})))

   if (stringr::str_count(rownames(LVs)[1], "_") ==1){
    loading.jj <- unlist(lapply(strsplit(colnames(LVs_embeddings), "_"), function(x){x[2]}))
  }else{
    loading.jj <- unlist(lapply(strsplit(colnames(LVs_embeddings), "_"), function(x){ paste0(x[-1], collapse = "_")}))
  }#else
  
    
   #extract from dav.res.list
   ##################################
   loading <- list()
   for (ii in 1:length(loading.uid)){
      loading[[ii]] <- dav.res.list[[loading.ii[ii]]]$Z[,loading.jj[ii]]

   }#for ii

   names(loading) <- loading.uid 

  return(list(LVs = LVs, 
              LVs_embeddings = LVs_embeddings, 
              slice_id = slice_id, 
              loading = loading))
  
  #LVs - the intermediate resutls after reciprocal_concat
  #LVs_embeddings - the product
  #slice_id - sample id per slice
  
}#Horizontal.Integration.first







#' Concatenate multiple gene expression matrix together 
#' @export
p.concatenate <- function(gene.exp.list, slice_id = NULL){
   
  #figure out the gene names
  #######################################################
  gene.name.list <- lapply(gene.exp.list, rownames)
  gene.name.inuse <- p.intersect(gene.name.list)
  
  gene.exp.list.joint <- lapply(gene.exp.list, function(x){x[gene.name.inuse,]})
  
  #concatenate
  #######################################################
  res <- p.cbind(gene.exp.list.joint)
  
  #rename the column names
  #######################################################
  slice_id_numeric <- as.character(rep(1:length(gene.exp.list), times = unlist(lapply(gene.exp.list, ncol)) ))
  if (!is.null(slice_id)){
     if (!all(slice_id == slice_id_numeric)){
        stop("Need to check slice_id carefully.")
     }else{
       slice_id <- slice_id_numeric
     }#else
  }#if
  
  colnames(res) <- paste0(slice_id, "_", colnames(res))
  return(res)
  
}#p.concatenate






p.cbind <- function(x){
  if (length(x)==1){
    return(x[[1]])
  }else if (length(x)==2){
    return(cbind(x[[1]], x[[2]]))
  }else if (length(x)>2){
    res <- cbind(x[[1]], x[[2]])
    for (ii in 3:length(x)){
      res <- cbind(res, x[[ii]])
    }#for ii
    return(res)
  }#else if
}#p.cbind


p.rbind <- function(x){
  if (length(x)==1){
    return(x[[1]])
  }else if (length(x)==2){
    return(rbind(x[[1]], x[[2]]))
  }else if (length(x)>2){
    res <- rbind(x[[1]], x[[2]])
    for (ii in 3:length(x)){
      res <- rbind(res, x[[ii]])
    }#for ii
    return(res)
  }#else if
}#p.rbind



