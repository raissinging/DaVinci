#' Adaptive Leiden/Louvain clustering
#'
#' When the number of clusters is set, run louvain/leiden clustering in an adaptive manner.
#' 
#' @export
leiden_adaptive <- function(nn, 
                            num.of.cluster = 5, 
                            resolution.start = 0.5, 
                            adaptive.size = 2, 
                            method = "louvain", 
                            verbose = T, 
                            full = F,
                            remove.singleton = NULL){
  
  
  #partition <- leiden::leiden(nn, resolution_parameter = resolution.start)
  if (method == "leiden"){
    partition <- Seurat:::RunLeiden(nn, resolution.parameter = resolution.start, method = "matrix") 
  }else if (method =="louvain"){
    partition <- Seurat:::RunModularityClustering(nn, resolution = resolution.start, algorithm = 1, print.output = F) 
  }#else if

 

    #correct for singleton
    if (!is.null(remove.singleton)){
        
        p.summary <- table(partition)

        p.singleton <- as.numeric(names(p.summary)[p.summary<3])
        message("Resolve singleton")

        if (length(p.singleton) > 1){
          #singleton exists
          p.cluster <- as.numeric(names(p.summary)[p.summary>=3])
          
          p.cluster.center <- c()
          #calculate the mean positions for p.cluster
          for (ii in 1:length(p.cluster)){
              temp.index <- which(partition==p.cluster[ii])
              p.cluster.center <- rbind(p.cluster.center, apply(remove.singleton[temp.index,], 2, mean))
          }#for ii
            
          #create a nearest neighbor index
          prebuild.index <- BiocNeighbors::buildIndex(p.cluster.center, BNPARAM = BiocNeighbors::AnnoyParam() )
          
          #assign singleton based on the distance
          p.singleton.index <- which(partition %in% p.singleton)
          p.singleton.to_query <- remove.singleton[p.singleton.index,]
          
          to_query.res <- FindNN(query = p.singleton.to_query, number.of.NN = 1, prebuild.index = prebuild.index)
          
          partition[p.singleton.index] <- p.cluster[to_query.res$index]
          
        }#if 

    }#if 

#wrapper, directly provide results for a given resolution
  if (is.null(num.of.cluster)){

     pp <- as.character(partition)
    if (!is.null(rownames(nn))){
        names(pp) <- rownames(nn)
    }#if

    if (full){
      return(list(partition = pp, reso = resolution.start))
    }else{
      return(pp)
    }#else

  }#if



  if (length(unique(partition)) < num.of.cluster){
    reso.left <- resolution.start
    reso.right <- resolution.start*adaptive.size
    reso <- reso.right
  }else if (length(unique(partition)) > num.of.cluster){
    reso.left <- resolution.start/adaptive.size
    reso.right <- resolution.start
    reso <- reso.left
  }else{
    reso <- resolution.start
    
    pp <- as.character(partition)
    if (!is.null(rownames(nn))){
        names(pp) <- rownames(nn)
    }#if

    if (full){
      return(list(partition = pp, reso = reso))
    }else{
      return(pp)
    }#else
    
  }#else
  
  
    
    while(length(unique(partition)) != num.of.cluster){
 
      if (reso < 0){
        stop("The resolution parameter is not valid.")
      }#if
      
      if (abs(reso.left-reso.right) < 1e-4){
        if (verbose){
          message("Re-initalize.") 
        }#if verbose
        
        return( leiden_adaptive(nn, num.of.cluster, resolution.start = resolution.start+0.1, adaptive.size, method, verbose, full) )
        
      }else{
      
        if (method == "leiden"){
          #partition <- leiden::leiden(nn, resolution_parameter = reso)
          partition <- Seurat:::RunLeiden(nn, resolution.parameter = reso, method = "matrix")  
        }else if (method == "louvain"){
          partition <- Seurat:::RunModularityClustering(nn, resolution = reso, algorithm = 1, print.output = F)  
        }#else if
        
        #correct for singleton
        if (!is.null(remove.singleton)){
               p.summary <- table(partition)

                p.singleton <- as.numeric(names(p.summary)[p.summary<3])
                message("Resolve singleton")

                if (length(p.singleton) > 1){
                  #singleton exists
                  p.cluster <- as.numeric(names(p.summary)[p.summary>=3])
                  
                  p.cluster.center <- c()
                  #calculate the mean positions for p.cluster
                  for (ii in 1:length(p.cluster)){
                      temp.index <- which(partition==p.cluster[ii])
                      p.cluster.center <- rbind(p.cluster.center, apply(remove.singleton[temp.index,], 2, mean))
                  }#for ii
                    
                  #create a nearest neighbor index
                  prebuild.index <- BiocNeighbors::buildIndex(p.cluster.center, BNPARAM = BiocNeighbors::AnnoyParam() )
                  
                  #assign singleton based on the distance
                  p.singleton.index <- which(partition %in% p.singleton)
                  p.singleton.to_query <- remove.singleton[p.singleton.index,]
                  
                  to_query.res <- FindNN(query = p.singleton.to_query, number.of.NN = 1, prebuild.index = prebuild.index)
                  
                  partition[p.singleton.index] <- p.cluster[to_query.res$index]
                  
                }#if 

        }#if 



        if ( length(unique(partition)) < num.of.cluster){
          if (verbose){
            message("Explode")
          }#if verbose
          
          if (reso == reso.left){
            reso.left <- reso.left
            reso <- (reso.left+reso.right)/2
            reso.right <- reso.right
            
          }else if (reso < reso.right){
            reso.left <- reso
            reso <- (reso+reso.right)/2
            reso.right <- reso.right
          }else if (reso == reso.right){
            reso <- reso.right*adaptive.size
            reso.left <- reso.right
            reso.right <- reso
          }#else if
          
        }else if ( length(unique(partition)) > num.of.cluster){
          if (verbose){
            message("Shrink")
          }#if verbose          
          
          if (reso == reso.left){
            reso <- reso/adaptive.size
            reso.right <- reso.left
            reso.left <- reso
            
          }else if (reso < reso.right){
            reso.right <- reso
            reso <- (reso.left+reso)/2
            reso.left <- reso.left
            
          }else if (reso == reso.right){
            reso.left <- reso.left
            reso <- (reso.left+reso.right)/2
            reso.right <- reso.right
            
          }#else if
          
        }else{
          
          pp <- as.character(partition)
          if (!is.null(rownames(nn))){
              names(pp) <- rownames(nn)
          }#if
          
          if (full){
              return(list(partition = pp, reso = reso))
          }else{
              return(pp)
          }#else
        
        }#else
        
        if (verbose){
          cat(paste0(reso, " ", reso.left, " ", reso.right, "\n"))  
        }#if verbose  
        
      }#else
      
    }#while
    
 
}#leiden_adaptive






  
self_extract <- function(proj, ids, opt.in = "Z"){
    
    res <- c()
    
    ii.list <- as.numeric(unlist(lapply(strsplit(ids, "_"), function(x){x[1]})))
    jj.list <- as.numeric(unlist(lapply(strsplit(ids, " "), function(x){x[2]})))
    
    if (opt.in=="Z"){
      #extract B      
      for (ii in 1:length(ii.list)){
        res <- rbind(res, proj[[ii.list[ii]]]$B[jj.list[ii],])
      }#for ii
      rownames(res) <- ids
      colnames(res) <- colnames(proj[[1]]$B)
      
    }else if (opt.in == "B"){
      #extract Z
      for (ii in 1:length(ii.list)){
        res <- rbind(res, proj[[ii.list[ii]]]$Z[,jj.list[ii]])
      }#for ii
      rownames(res) <- ids
      colnames(res) <- rownames(proj[[1]]$Z)

    }#else if
    
    return(res)
    
}#self_extract
  
  

#' Consensus Decomposition
#'
#' Figure out the consensus of latent variables
#' 
#' @export
self_deco <- function(proj,  
                      LVs.filter.thr = 0.9, 
                      freq = 1, 
                      opt = "B",
                      verbose = F){
  
  num.of.slice <- length(proj)
  
  #row names and colnames
  row.names <- c()
  LVs <- c()
  
  if (opt == "B"){
    
    col.names <- colnames(proj[[1]]$B)
    
    for (ii in 1:num.of.slice){
      LVs <- rbind(LVs, proj[[ii]]$B)
      row.names <- c(row.names, paste0(ii, "_", rownames(proj[[ii]]$B)) )
      if (!all(col.names==colnames(proj[[ii]]$B))){
        warning("column samples are mismatched at ", ii)
        break
      }#if
    }#for ii
    

  }else if (opt == "Z"){
    col.names <- rownames(proj[[1]]$Z)
    
    for (ii in 1:num.of.slice){
      LVs <- rbind(LVs, t(proj[[ii]]$Z))
      row.names <- c(row.names, paste0(ii, "_", colnames(proj[[ii]]$Z)) )
      if (!all(col.names==rownames(proj[[ii]]$Z))){
        warning("column samples are mismatched at ", ii)
        break
      }#if
    }#for ii

  }#else if opt == "Z"
  
  
  #construct the LVs
  ###########################################
  
  rownames(LVs) <- row.names
  colnames(LVs) <- col.names
    
  
  #filter by correlation
  ###########################################
  
      #in case there are 0-sd ones       
      include.index <- which(apply(LVs, 1, sd) > 1e-7)
      LVs <- LVs[include.index,]


      #only selected the most highly correlated ones
      cor.res <- cor(t(LVs))
      
      #converted to a connected graph based on the threshold
      adj_matrix <- cor.res > LVs.filter.thr
      adj_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = F)
      components_list <- igraph::components(adj_graph)
      components_selected <- which(components_list$csize>freq)
      components_selected <- components_selected[order(components_list$csize[components_selected] , decreasing = T)]
      
      if (verbose){
        print(sort(components_list$csize, decreasing = T))
      }#if verbose      
      
      #output the module size
      components.size <- components_list$csize[components_selected]
          
      #only pick up one from correlated blocks
      include.names <- c()
      for (ii in 1:length(components_selected)){
        #select the one with smaller index
        include.names <- c(include.names, names(components_list$membership)[which(components_list$membership==components_selected[ii])][1] )
      }#for ii
      
      LVs <- LVs[include.names,]



#extract the corresponding ones
  if (opt == "B"){
    LVs.pair <- self_extract(proj, rownames(LVs), opt.in = "B")
  }else if (opt == "Z"){
    LVs.pair <- self_extract(proj, rownames(LVs), opt.in = "Z")
  }#else if
  
  return(list(LVs=LVs, 
              LVs.pair = LVs.pair, 
              components.size = components.size))
}#self_deco






#radius based, nearest neighbor based
#continuous smoothing or voting
#nearest neighbor: label imputation
#nearest neighbor: value smoothing
refinement <- function(y, 
                       coor, 
                       neighbor.option="KNN", 
                       neighbor.arg = 6, 
                       radius.arg = 10, 
                       tasks = "discrete", 
                       self = F, 
                       random.seed = 1){

  set.seed(random.seed)
  #get the neighbor index per spot
  #neighbor.index <- list()

  if (neighbor.option == "KNN"){
    nn.res <- Seurat::FindNeighbors(coor, k.param = neighbor.arg+1, return.neighbor = T, compute.SNN = F, verbose = F)
    nn.idx <- Seurat::Indices(nn.res)
    neighbor.index <- lapply(seq_len(nrow(coor)), function(i){
      idx <- nn.idx[i, ]
      idx[idx != i]   #drop self
    })

  }else if (neighbor.option == "Radius"){
    #decide neighbors based on the radius
    coor.dist <- as.matrix(dist(coor))
    diag(coor.dist) <- Inf

    neighbor.index <- apply(coor.dist, 1, function(x){which(x < radius.opt)})
  }#else if


  if (tasks == "discrete"){

    y.table <- lapply(neighbor.index, function(x){table(y[x])})
    y.smooth <- unlist( lapply(y.table, function(x){ if (max(x)>neighbor.arg/2){names(which.max(x))}else{NA} } )  )
    y.smooth[is.na(y.smooth)] <- y[is.na(y.smooth)]

  }else if (tasks == "continuous"){

    if (self){
      y.smooth <- (unlist(lapply(neighbor.index, function(x){mean(y[x])}))+y)/2
    }else{
      y.smooth <- unlist(lapply(neighbor.index, function(x){mean(y[x])}))
    }#else

  }#else if

  return(y.smooth)

}#refinement




#' Smooth out by LV
#'
#' Smooth out by LV
#' 
#' @export
refinement.batch <- function(data, 
                             coor, 
                             neighbor.option="KNN", 
                             neighbor.arg = 6, 
                             radius.arg = 10, 
                             tasks = "discrete", 
                             self = F){
  if (ncol(data) < nrow(data)){
    data <- t(data)
  }#if

  for (ii in 1:nrow(data)){
    data[ii, ] <- refinement(data[ii,], coor = coor, neighbor.option = neighbor.option, neighbor.arg = neighbor.arg, radius.arg =radius.arg, tasks = tasks, self = self)
  }#for ii

  return(data)
}#refinement.batch







#' Proximity-dependent gene discovery
#'
#' Proximity-dependent gene discovery
#' 
#' @export
Proximity.Dependency <- function(mat, 
                                 coor, 
                                 niche, 
                                 Query = list("0+3"=c("0","3")), 
                                 Anchor = list("1+2"=c("1","2")), 
                                 n_neighbors = 5,
                                 meth = "p",
                                 fn = mean){
  
  #create index for Anchor
  ###################################
  anchor.list <- list()
  
  for (ii in 1:length(Anchor)){
    include.index <- which(niche %in% Anchor[[ii]])
    
    prebuild.index.modal <- BiocNeighbors::buildIndex( coor[names(niche)[include.index],], BNPARAM = BiocNeighbors::AnnoyParam() )  
    anchor.list[[ii]] <- prebuild.index.modal
  }#for ii
  names(anchor.list) <- names(Anchor)
  
  
  #query
  ###################################
  query.list <- list()
  for (ii in 1:length(Query)){
    
    #print(ii)
    to_query <- coor[ names(niche)[which(niche %in% Query[[ii]])], ]
    
    temp <- list()
    for (jj in 1:length(Anchor)){
      modal.neighbor <- FindNN(query=to_query, number.of.NN = n_neighbors, prebuild.index = anchor.list[[jj]])  
      
      buffer <- apply(modal.neighbor$distance,1, fn)
      names(buffer) <- rownames(to_query)
      
      temp[[jj]] <- buffer
    }#for jj
    names(temp) <- names(Anchor)
    
    query.list[[ii]] <- temp
  }#for ii
  
  names(query.list) <- names(Query)
  
  
  #correlation test
  ###################################
  cor.list <- list()
  
  for (ii in 1:length(query.list)){
    
    cor.temp.list <- list()
    
    for (jj in 1:length(anchor.list)){
      
      dd <- query.list[[ii]][[jj]]
      yy <- mat[,names(dd)]
      
      
      calc <- apply(yy, 1, function(x){
        cor.test(dd, x, method = meth)
      })
      
      p.val <- sapply(calc, "[[", "p.value")
      cor.val <- sapply(calc, "[[", "estimate")
      q.val <- p.adjust(p.val, method = "BH")
      calc.summary <- data.frame(gene = rownames(mat),
                                 cor.val = cor.val,
                                 p.val = p.val,
                                 q.val = q.val
                                 )
      rownames(calc.summary) <- rownames(mat)
      calc.summary <- calc.summary[order(calc.summary$q.val, decreasing = F),]
      cor.temp.list[[jj]] <- calc.summary
      
    }#for jj
    
    names(cor.temp.list) <- names(anchor.list)
    cor.list[[ii]] <- cor.temp.list
    
  }#for ii
  
  names(cor.list) <- names(query.list)
  
  return(list(Query.dist = query.list, cor.list = cor.list))
  
}#Proximity.Dependency




#' Proximity-dependent gene visualization
#'
#' Proximity-dependent gene visualization
#' 
#' @export
Proximity.Dependency.scatter <- function(gene.name, 
                                         mat,
                                         dist,
                                         x.lab = "Distance to all the other niches",
                                         y.lab = "Expression levels"){
  
  dd <- dist
  yy <- mat[gene.name, names(dd)]
  
  to_plot <- data.frame(x = dd, y = yy)
  
  pp <- ggplot(to_plot, aes(x=x,y=y))+
    geom_point(size = 2)+
    theme_pubr(base_size = 30)+
    geom_smooth(method = "lm")+
    stat_cor(size=8, label.x = 5, label.y = 0.2)+
    xlab(x.lab)+
    ylab(y.lab)+
    ggtitle(gene.name)
  
  return(pp)
  
}#Proximity.Dependency.scatter





#Other parameters
#save.label <- "rna" #rna, peak as input        #preprocess
#normaliza.version <- "gene" #should be used when processing the data     #preprocess
# - rna, adt -> gene
# - atac -> gene, LSI_1, LSI_2, LSI_3

#clustering input parameters
# - L2.opt <- "Yes" #Yes, No - whether L2norm before finding the nearest neighbors or mclust
# - k.opt <- 40 #neighboring paramters for louvain

#' Cross-sample (Horizontal) integration wrapper
#'
#' Cross-sample (Horizontal) integration
#' 
#' @export
Horizontal.Integration.Assemble <- function(
    dataset.opts = NULL,
    Y.list,
    L.list, 
    coor.list,
    k.arg.list = c(4, 6, 8, 10, 11, 12, 13, 15, 16, 18, 20, 22, 24), #the parameters for self-contrastive learning step
    L4.arg = 50, #parameters
    input.opt = "FineTune", #OnlyDeco, FineTune
    h.opt = "first",   #default, first
    mod.opt = "all",   #all, common
    L2.in = "default", #default, L2norm, L2norm.joint
    random.seed = 123,
    smooth = F
    ){

  
  #sanity check
  #########################################################
  if (is.null(dataset.opts)){
    stop("Please include a vector of character names corresponding to the slices.")
  }#if
  
  
  #self-contrastive learning step
  #########################################################
  exhaustive.list <- list()

  #ptm <- proc.time()
  for (ii in 1:length(Y.list)){
    
    message(paste0("Working on Sample ", ii, " / ", length(Y.list)))
    proj <- list()
    s0 <- NULL  
    
    
    count <- 0
    for (k.arg in k.arg.list){
      message("Working on k=", k.arg)
      
      if (is.null(s0)){

        ICAp.res <- manifoldDecomp_adaptive(Y.list[[ii]], 
                                            L.list[[ii]], 
                                            k = k.arg, 
                                            L4 = L4.arg, 
                                            L4_adaptive = 2, 
                                            to_drop = T, 
                                            save.complete = T,
                                            verbose = F, 
                                            random.seed = random.seed)  
      }else{

        ICAp.res <- manifoldDecomp_adaptive(Y.list[[ii]], 
                                            L.list[[ii]], 
                                            k = k.arg, 
                                            L4 = L4.arg, 
                                            L4_adaptive = 2, 
                                            to_drop = T, 
                                            save.complete = T, 
                                            shur0 = s0, 
                                            verbose = F, 
                                            random.seed = random.seed)  
      }#else
      
      #saveRDS(ICAp.res, file = paste0(output.folder, "/", save.label, "@", normalize.version, "@k=", k.arg,".RDS") )
      count <- count+1
      proj[[count]] <- ICAp.res
      
      if (is.null(s0)){
        s0 <- ICAp.res$shur0
      }#if
      
    }#for k.arg

    exhaustive.list[[ii]] <- proj  
  }#for ii
  #print(proc.time()-ptm)



    #self contrastive learning step
    #######################################################################################
  if (length(k.arg.list) == 1){
      
      #create the object needed for the downstream analysis
      embed.list <- list()
      embed.list.finetune <- list()
      dav.res.list <- list()

      for (ii in 1:length(Y.list)){
          embed.list[[ii]] <- exhaustive.list[[ii]][[1]]  
          embed.list.finetune[[ii]] <- exhaustive.list[[ii]][[1]]  
          dav.res.list[[ii]] <- exhaustive.list[[ii]][[1]]
      }#for ii
      
  }else{
    
    #self contrastive learning
    #############################
    embed.list <- list()
    for (ii in 1:length(exhaustive.list)){

      embed.list[[ii]] <- self_deco(exhaustive.list[[ii]], 
                                    LVs.filter.thr = 0.9, 
                                    freq = 1, 
                                    opt = "B")
    }#for ii
    

    #run again to finetune
    #############################
    
    embed.list.finetune <- list()
    for (ii in 1:length(embed.list)){
      
      embed.list.finetune[[ii]] <- manifoldDecomp_adaptive(Y.list[[ii]], 
                            L.list[[ii]], 
                            k = nrow(embed.list[[ii]]$LVs), 
                            B = embed.list[[ii]]$LVs, 
                            L4 = L4.arg, 
                            L4_adaptive = 2, 
                            to_drop = T, 
                            save.complete = T, 
                            shur0 = exhaustive.list[[ii]][[1]]$shur0, 
                            verbose = F, 
                            random.seed = random.seed)
            
    }#for ii


    #Finetune after self-contrastive learning
    if (input.opt == "OnlyDeco"){
    dav.res.list <- list()
    
    for (ii in 1:length(embed.list.finetune)){
      
      ICAp.res <- embed.list.finetune[[ii]]
      #construct the pseudo ICAp.res objects
      tmp <- list(Z = t(embed$LVs.pair), 
                  B = embed$LVs, 
                  L4 = ICAp.res$L4, 
                  shur0 = ICAp.res$shur0, 
                  L1 = ICAp.res$L1, 
                  L2 = ICAp.res$L2)
      #cor.res <- cor(t(embed$LVs), t(ICAp.res$B))
      
      dav.res.list[[ii]] <- tmp
    }#for ii
    
    }else if (input.opt == "FineTune"){
      
      dav.res.list <- embed.list.finetune
      
    }#else if
  }#else



  if (length(dataset.opts)==1){

      #no need to integrate here, directly return here
      mat <- t(dav.res.list[[1]]$B)
      mat <- as.matrix(mat)
      mat.slice.id <- rep(dataset.opts[1], nrow(mat))
      #no integration.res
      return(list(Y.list = Y.list, 
                  coor.list = coor.list, 
                  L.list = L.list, 
                  embed.list = embed.list, 
                  embed.list.finetune = embed.list.finetune, 
                  dataset.opts = dataset.opts, 
                  mat.slice.id = mat.slice.id, 
                  mat = mat,
                  exhaustive.list = exhaustive.list))

  }else{
    #Horizontal integration
    ########################################################

    if (h.opt == "default"){
    
      integration.res <- Horizontal.Integration(Y.list, 
                                                L.list, 
                                                coor.list, 
                                                dav.res.list = dav.res.list, 
                                                LVs.filter.thr = 0.8, 
                                                mod = mod.opt, 
                                                remove.LV1 = F, 
                                                L2.option = L2.in)
      
    }else if (h.opt == "first"){
          
      integration.res <- Horizontal.Integration.first(Y.list, 
                                                      L.list, 
                                                      coor.list, 
                                                      dav.res.list = dav.res.list, 
                                                      LVs.filter.thr = 0.8, 
                                                      mod = mod.opt, 
                                                      remove.LV1 = F, 
                                                      L2.option = L2.in)
      
    }#else
    
    
    #cleanup the output
    #########################################################
    embed <- integration.res$LVs_embeddings
    #make sure the sample names can align
    sids <- unlist(lapply(strsplit( rownames(embed), "_"), function(x){paste0(x[-1], collapse = "_")}))
    mat.slice.id <- unlist(lapply(strsplit(rownames(embed),"_"), function(x){x[1]}))
    mat.slice.id <- dataset.opts[as.numeric(mat.slice.id)]

    rownames(embed) <- sids             
    mat <- as.matrix(embed)

    
    
    #smoothing first no matter used or not
    #########################################################
    if (!smooth){
        return(list(Y.list = Y.list, 
                    coor.list = coor.list, 
                    L.list = L.list, 
                    embed.list = embed.list, 
                    embed.list.finetune = embed.list.finetune, 
                    dataset.opts = dataset.opts, 
                    mat.slice.id = mat.slice.id, 
                    integration.res = integration.res, 
                    mat = mat,
                    exhaustive.list = exhaustive.list))
    }else{
      mat.smooth <- mat

    for (ii in 1:length(dataset.opts)){
      print(ii)
      
      subpart.index <- which(mat.slice.id==dataset.opts[ii])
      
      subpart <- mat[subpart.index,]
      subpart.smooth <- refinement.batch(subpart, 
                                         as.matrix(coor.list[[ii]])[rownames(subpart),], 
                                         neighbor.option = "KNN", 
                                         neighbor.arg = 8, 
                                         tasks = "continuous")
      
      mat.smooth[subpart.index,] <- (subpart+t(subpart.smooth))/2
      
    }#for ii
    
      return(list(Y.list = Y.list, 
                  coor.list = coor.list, 
                  L.list = L.list, 
                  embed.list = embed.list, 
                  embed.list.finetune = embed.list.finetune, 
                  dataset.opts = dataset.opts, 
                  mat.slice.id = mat.slice.id, 
                  integration.res = integration.res, 
                  mat = mat, 
                  mat.smooth = mat.smooth,
                  exhaustive.list = exhaustive.list))

    }#else

  }#else

}#Horizontal.Integration.Assemble




#Other parameters
#save.label <- "rna" #rna, peak as input        #preprocess
#normaliza.version <- "gene" #should be used when processing the data     #preprocess
# - rna, adt -> gene
# - atac -> gene, LSI_1, LSI_2, LSI_3

#clustering input parameters
# - L2.opt <- "Yes" #Yes, No - whether L2norm before finding the nearest neighbors or mclust
# - k.opt <- 40 #neighboring paramters for louvain

#' Cross-sample (Horizontal) integration wrapper
#'
#' Cross-sample (Horizontal) integration
#' 
#' @export
Horizontal.Integration.Assemble.Parallel <- function(
    dataset.opts = NULL,
    Y.list,
    L.list, 
    coor.list,
    k.arg.list = c(4, 6, 8, 10, 11, 12, 13, 15, 16, 18, 20, 22, 24), #the parameters for self-contrastive learning step
    L4.arg = 50, #parameters
    input.opt = "FineTune", #OnlyDeco, FineTune
    h.opt = "first",   #default, first
    mod.opt = "all",   #all, common
    L2.in = "default", #default, L2norm, L2norm.joint
    random.seed = 123,
    smooth = F,
    workers = 2, # number of cores/workers for parallelization
    ram = 15, # in GB max RAM usage 
    sequential = F,
    multi = "multicore" # multicore (mac and linux), multisession (for windows)
    ){ 

  
  #sanity check
  #########################################################
  if (is.null(dataset.opts)){
    stop("Please include a vector of character names corresponding to the slices.")
  }#if

  if (sequential){
    #sequential execution
    future::plan(future::sequential)
  }else{
    #parallel execution
    if (multi == "multisession"){
      future::plan(future::multisession, workers = workers)
    }else{
      future::plan(future::multicore, workers = workers)
    }
    pkgs_to_load <- if (multi == "multisession") c("DaVinci", "progressr") else NULL
  }

  options(future.globals.maxSize = ram * 1024^3)

  # progress bar settings 
  progressr::handlers("progress")

  options(
    progressr.enable = TRUE,
    progressr.enable_after = 0,      
    progressr.delay_stdout = FALSE   
  )

  progressr::handlers(progressr::handler_txtprogressbar(
    file = stdout(),
    enable_after = 0
  ))

  #self-contrastive learning step
  #########################################################
  
  exhaustive.list <- progressr::with_progress({
    
    p <- progressr::progressor(steps = length(Y.list) * length(k.arg.list))

    future.apply::future_lapply(seq_along(Y.list), function(ii) {
      
      message(paste0("Working on Sample ", ii, " / ", length(Y.list)))
      proj <- list()
      s0 <- NULL  
      
      
      count <- 0
      for (k.arg in k.arg.list){
        message("Working on k=", k.arg)
        
        if (is.null(s0)){

          ICAp.res <- manifoldDecomp_adaptive(Y.list[[ii]], 
                                              L.list[[ii]], 
                                              k = k.arg, 
                                              L4 = L4.arg, 
                                              L4_adaptive = 2, 
                                              to_drop = T, 
                                              save.complete = T,
                                              verbose = F, 
                                              random.seed = random.seed)  
        }else{

          ICAp.res <- manifoldDecomp_adaptive(Y.list[[ii]], 
                                              L.list[[ii]], 
                                              k = k.arg, 
                                              L4 = L4.arg, 
                                              L4_adaptive = 2, 
                                              to_drop = T, 
                                              save.complete = T, 
                                              shur0 = s0, 
                                              verbose = F, 
                                              random.seed = random.seed)  
        }#else
        
        count <- count+1
        proj[[count]] <- ICAp.res

        p(sprintf("Sample %d: k=%d", ii, k.arg))
        
        if (is.null(s0)){
          s0 <- ICAp.res$shur0
        }#if
        
      }#for k.arg

      return(proj)
    }, future.seed = TRUE, future.packages = pkgs_to_load)

  })


    #self contrastive learning step
    #######################################################################################
  if (length(k.arg.list) == 1){
      
      #create the object needed for the downstream analysis
      embed.list <- list()
      embed.list.finetune <- list()
      dav.res.list <- list()

      for (ii in 1:length(Y.list)){
          embed.list[[ii]] <- exhaustive.list[[ii]][[1]]  
          embed.list.finetune[[ii]] <- exhaustive.list[[ii]][[1]]  
          dav.res.list[[ii]] <- exhaustive.list[[ii]][[1]]
      }#for ii
      
  }else{
    
    #self contrastive learning
    #############################
    
    embed.list <- future.apply::future_lapply(seq_along(exhaustive.list), function(ii) {
      self_deco(exhaustive.list[[ii]], 
                LVs.filter.thr = 0.9, 
                freq = 1, 
                opt = "B")
    }, future.seed = TRUE, future.packages = pkgs_to_load)
    

    #run again to finetune
    #############################
    
    embed.list.finetune <- future.apply::future_lapply(seq_along(embed.list), function(ii) {
      manifoldDecomp_adaptive(Y.list[[ii]], 
                              L.list[[ii]], 
                              k = nrow(embed.list[[ii]]$LVs), 
                              B = embed.list[[ii]]$LVs, 
                              L4 = L4.arg, 
                              L4_adaptive = 2, 
                              to_drop = T, 
                              save.complete = T, 
                              shur0 = exhaustive.list[[ii]][[1]]$shur0, 
                              verbose = F, 
                              random.seed = random.seed)
    }, future.seed = TRUE, future.packages = pkgs_to_load)


    #Finetune after self-contrastive learning
    if (input.opt == "OnlyDeco"){
    
    dav.res.list <- future.apply::future_lapply(seq_along(embed.list.finetune), function(ii) {
      ICAp.res <- embed.list.finetune[[ii]]
      #construct the pseudo ICAp.res objects
      list(Z = t(embed$LVs.pair), 
           B = embed$LVs, 
           L4 = ICAp.res$L4, 
           shur0 = ICAp.res$shur0, 
           L1 = ICAp.res$L1, 
           L2 = ICAp.res$L2)
    }, future.seed = TRUE, future.packages = pkgs_to_load)
    
    }else if (input.opt == "FineTune"){
      
      dav.res.list <- embed.list.finetune
      
    }#else if
  }#else



  if (length(dataset.opts)==1){

      #no need to integrate here, directly return here
      mat <- t(dav.res.list[[1]]$B)
      mat <- as.matrix(mat)
      mat.slice.id <- rep(dataset.opts[1], nrow(mat))
      #no integration.res
      return(list(Y.list = Y.list, 
                  coor.list = coor.list, 
                  L.list = L.list, 
                  embed.list = embed.list, 
                  embed.list.finetune = embed.list.finetune, 
                  dataset.opts = dataset.opts, 
                  mat.slice.id = mat.slice.id, 
                  mat = mat,
                  exhaustive.list = exhaustive.list))

  }else{
    #Horizontal integration
    ########################################################

    if (h.opt == "default"){
    
      integration.res <- Horizontal.Integration(Y.list, 
                                                L.list, 
                                                coor.list, 
                                                dav.res.list = dav.res.list, 
                                                LVs.filter.thr = 0.8, 
                                                mod = mod.opt, 
                                                remove.LV1 = F, 
                                                L2.option = L2.in)
      
    }else if (h.opt == "first"){
          
      integration.res <- Horizontal.Integration.first(Y.list, 
                                                      L.list, 
                                                      coor.list, 
                                                      dav.res.list = dav.res.list, 
                                                      LVs.filter.thr = 0.8, 
                                                      mod = mod.opt, 
                                                      remove.LV1 = F, 
                                                      L2.option = L2.in)
      
    }#else
    
    
    #cleanup the output
    #########################################################
    embed <- integration.res$LVs_embeddings
    #make sure the sample names can align
    sids <- unlist(lapply(strsplit( rownames(embed), "_"), function(x){paste0(x[-1], collapse = "_")}))
    mat.slice.id <- unlist(lapply(strsplit(rownames(embed),"_"), function(x){x[1]}))
    mat.slice.id <- dataset.opts[as.numeric(mat.slice.id)]

    rownames(embed) <- sids             
    mat <- as.matrix(embed)

    
    
    #smoothing first no matter used or not
    #########################################################
    if (!smooth){
        return(list(Y.list = Y.list, 
                    coor.list = coor.list, 
                    L.list = L.list, 
                    embed.list = embed.list, 
                    embed.list.finetune = embed.list.finetune, 
                    dataset.opts = dataset.opts, 
                    mat.slice.id = mat.slice.id, 
                    integration.res = integration.res, 
                    mat = mat,
                    exhaustive.list = exhaustive.list))
    }else{
      mat.smooth <- mat

    smooth_updates <- future.apply::future_lapply(seq_along(dataset.opts), function(ii) {
      print(ii)
      
      subpart.index <- which(mat.slice.id==dataset.opts[ii])
      subpart <- mat[subpart.index,]
      subpart.smooth <- refinement.batch(subpart, 
                                         as.matrix(coor.list[[ii]])[rownames(subpart),], 
                                         neighbor.option = "KNN", 
                                         neighbor.arg = 8, 
                                         tasks = "continuous")
      
      return(list(idx = subpart.index, 
                  val = (subpart+t(subpart.smooth))/2))
    }, future.seed = TRUE, future.packages = pkgs_to_load)
    
    for (update in smooth_updates) {
      mat.smooth[update$idx, ] <- update$val
    }
    
      return(list(Y.list = Y.list, 
                  coor.list = coor.list, 
                  L.list = L.list, 
                  embed.list = embed.list, 
                  embed.list.finetune = embed.list.finetune, 
                  dataset.opts = dataset.opts, 
                  mat.slice.id = mat.slice.id, 
                  integration.res = integration.res, 
                  mat = mat, 
                  mat.smooth = mat.smooth,
                  exhaustive.list = exhaustive.list))

    }#else

  }#else

  future::plan(sequential) # Explicitly close multisession workers by switching plan

}#Horizontal.Integration.Assemble










#' Cross-modality (Vertical) integration wrapper
#'
#' Cross-modality (Vertical) integration
#' 
#' @export
Vertical.Integration.Assemble <- function(
    dataset.opts,
    Y1.list,
    Y2.list,
    L.list, 
    coor.list,
    k.arg.list = c(4, 6, 8, 10, 11, 12, 13, 15, 16, 18, 20, 22, 24), #the parameters for self-contrastive learning step
    L4.arg = 50, #parameters
    input.opt = "FineTune", #OnlyDeco, FineTune
    h.opt = "first",   #default, first
    mod.opt = "all",   #all, common
    L2.in = "default", #default, L2norm, L2norm.joint
    random.seed = 123,
    smooth = F,
    n_neighbors = 40,    #vertical integration parameters - optimized
    n_neighbors_large = 50,  #vertical integration parameters - optimized
    L2_norm = "column",  #vertical integration parameters - optimized
    library.size.1 = NULL,
    library.size.2 = NULL,
    library.type.1 = "rna", 
    library.type.2 = "atac",
    baseline.1 = NULL,
    baseline.2 = NULL,
    verbose = F
   ){
    
  
  #sanity check
  #########################################################
  if (is.null(dataset.opts)){
    stop("Please include a vector of character names corresponding to the subfolders.")
  }#if
  
  
  
  #first run through Horizontal.Integration.Assembele
  #########################################################
  #Modality 1
  message("Working on Modality 1")
  ptm <- proc.time()
  HI.1.res <- Horizontal.Integration.Assemble(
                dataset.opts = dataset.opts,
                Y.list  = Y1.list,
                L.list  = L.list,
                coor.list = coor.list,
                k.arg.list = k.arg.list,
                L4.arg = L4.arg,
                input.opt = input.opt,
                h.opt = h.opt,
                mod.opt = mod.opt,
                L2.in = L2.in,
                random.seed = random.seed,
                smooth = smooth
  )
  if (verbose){
      print(proc.time()-ptm)
      #4455s
  }#if
    
  #Modality 2
  message("Working on Modality 2")
  ptm <- proc.time()
  HI.2.res <- Horizontal.Integration.Assemble(
                dataset.opts = dataset.opts,
                Y.list  = Y2.list,
                L.list  = L.list,
                coor.list = coor.list,
                k.arg.list = k.arg.list,
                L4.arg = L4.arg,
                input.opt = input.opt,
                h.opt = h.opt,
                mod.opt = mod.opt,
                L2.in = L2.in,
                random.seed = random.seed,
                smooth = smooth
  )
  if (verbose){
      print(proc.time()-ptm)
      #5427s
  }#if verbose

  #vertical integration
  #########################################################
  #make sure names aligned
  ##########################
  input1 <- HI.1.res$mat
  input2 <- HI.2.res$mat

  message("Make sure names aligned.")
  #all(paste0(HI.2.res$mat.slice.id, "@", rownames(input1)) == paste0(HI.1.res$mat.slice.id, "@", rownames(input2)))
    
  if (!is.null(library.size.1)){
      if (all(paste0(HI.1.res$mat.slice.id, "@", rownames(input1)) == names(library.size.1))){
          message("Modality 1 aligned.")
      }#if
  }#if

  if (!is.null(library.size.2)){
      if (all(paste0(HI.2.res$mat.slice.id, "@", rownames(input2)) == names(library.size.2))){
          message("Modality 2 aligned.")
      }#if
  }#if

  mat.slice.id <- HI.1.res$mat.slice.id

  #make sure the barcodes are unique
  if (sum(duplicated(rownames(input1))) >0){
    rownames(input1) <- paste0(HI.1.res$mat.slice.id, "@", rownames(input1))
    if (sum(duplicated(rownames(input1))) ==0){
        message("Rename modality 1 to make sure barcodes are unique.")
    }else{
        stop("Modality 1 barcodes are not unique.")
    }#else
  }#if

  if (sum(duplicated(rownames(input2))) >0){
    rownames(input2) <- paste0(HI.2.res$mat.slice.id, "@", rownames(input2))
    if (sum(duplicated(rownames(input2))) ==0){
        message("Rename modality 2 to make sure barcodes are unique.")
    }else{
        stop("Modality 2 barcodes are not unique.")
    }#else
  }#if


  #start the vertical integation
  ##########################
  integration.res <- BiModalIntegration(
                            modal.1 = input1, 
                            modal.2 = input2, 
                            n_neighbors = n_neighbors, 
                            n_neighbors_large = n_neighbors_large, 
                            L2norm = L2_norm,
                            library.size.1 = library.size.1,
                            library.size.2 = library.size.2,
                            library.type.1 = library.type.1,
                            library.type.2 = library.type.2,
                            baseline.1 = baseline.1,
                            baseline.2 = baseline.2
                            )
  
  mat <- integration.res$snn.mat
  
  #assign names
  rownames(mat) <- names(integration.res$weight.1)
  colnames(mat) <- names(integration.res$weight.1)

  return(list(mat = mat, 
              mat.slice.id = mat.slice.id, 
              HI.1.res = HI.1.res, 
              HI.2.res = HI.2.res, 
              integration.res = integration.res))
    
}#Vertical.Integration.Assemble

#sids <- rownames(input1)
#mat.slice.id <- unlist(lapply(strsplit( sids,"_"), function(x){x[1]}))











#c(4, 6, 8, 9, 10, 11, 12, 14, 16, 18, 20)
#(1:15)/10
#L2.opt <- "Yes" #Yes, No - whether L2norm before finding the nearest neighbors or mclust
#k.opt <- 40 #neighboring paramters for louvain

#' Niche report
#'
#' Cluster the latent variables to get the niche labels
#' 
#' @export
Niche.report <- function(
  mat, 
  method = "louvain", 
  L2.opt = T,                #key parameters
  zscore = F, 
  mclust.model = "EEE", 
  num.cluster.args = NULL, 
  num.neighbors = 40,        #key parameters
  resolution.args = NULL, 
  random.seed = 1,
  verbose = F
  ){
  
  message("mclust, leiden and louvain have been preset.")
  
  par.list <- list()  
  set.seed(random.seed)
  
  #flip the coordinates if needed
  if (ncol(mat) > nrow(mat)){
    mat <- t(mat)
    message("There are ", nrow(mat), " spots in the latent space (n = ", ncol(mat), ").")
  }#if
  
  
  if (is.null(num.cluster.args) & is.null(resolution.args)){
    stop("You need to set at least one of the two parameters: num.cluster.args or resolution.args, in concordnace with the method argument.") 
  }#if 
  
  if (method=="mclust" & is.null(num.cluster.args)){
    stop("You are using mclust, but the number of clusters is not specified. You can specify a range of values, e.g., c(4, 6, 8, 9, 10, 11, 12, 14, 16, 18, 20) or a single value.")
  }#if
  
  if (method=="louvain" & is.null(num.cluster.args) & is.null(resolution.args) ){
    stop("You are using louvain algorithm, but neither of the number of clusters nor the resolution parameter is set.")
  }#if
  
  if (method=="leiden" & is.null(num.cluster.args) & is.null(resolution.args) ){
    stop("You are using leiden algorithm, but neither of the number of clusters nor the resolution parameter is set.")
  }#if
  
  
  
  #main body starts here
  ###############################################################
  if (method=="mclust" & !is.null(num.cluster.args)){
    
    #call mclust
    res.pre <- Cluster.Finetune.pre(mat,method = "mclust", L2norm = L2.opt, mclust.num.args = num.cluster.args)
    args <- res.pre$args
    par.list <- res.pre$par.list
    
  }else{
    
    #leiden and louvain
    #construct the graph
    
    if (L2.opt){
      snn.res <- Seurat::FindNeighbors(L2Norm(mat, MARGIN = 2), k.param = num.neighbors, return.neighbor = F, compute.SNN = T, verbose = F)
      #snn.res.smooth <- Seurat::FindNeighbors(L2Norm(mat.smooth, MARGIN = 2), k.param = num.neighbors, return.neighbor = F, compute.SNN = T, verbose = F)
      
    }else{
      snn.res <- Seurat::FindNeighbors(mat, k.param = num.neighbors, return.neighbor = F, compute.SNN = T, verbose = F) 
      #snn.res.smooth <- Seurat::FindNeighbors(mat.smooth, k.param = num.neighbors, return.neighbor = F, compute.SNN = T, verbose = F)  
    }#else
    
    if (is.null(num.cluster.args) & !is.null(resolution.args)){
      
      #resolution
      res.pre <- Cluster.Finetune.pre(snn.res$snn, method = method, ld.resolution.args = resolution.args)
      
      args <- res.pre$args
      par.list <- res.pre$par.list
      
    }else{
      
      #number of clusters
      par.list <- list()
      
      for (pp in 1:length(num.cluster.args)){
        par.list[[pp]] <- leiden_adaptive(snn.res$snn, num.of.cluster = num.cluster.args[pp], resolution.start = 0.4, adaptive.size = 2, method = method, verbose = verbose)
        names(par.list[[pp]]) <- rownames(mat)
      }#for pp
      
      args <- num.cluster.args
      
    }#else
    
  }#else - using louvain or leiden
  
  return(list(args = args, niche.list = par.list))
  
}#Niche.report








#' Niche finetune
#'
#' Apply local smoothing on the output of Niche.report().
#' 
#' @export
Niche.report.finetune <- function(par.list, dataset.opts, coor.list, mat.slice.id, neighbor.arg = 8){
  
  #refine every partition
  par.list.refine <- par.list
  
  for (ii in 1:length(dataset.opts)){
    
    subpart.index <- which(mat.slice.id==dataset.opts[ii])
    
    for (jj in 1:length(par.list.refine)){
      
      p.part <- par.list.refine[[jj]][subpart.index]
      p.part <- refinement(p.part, as.matrix(coor.list[[ii]])[names(p.part),], neighbor.option = "KNN", neighbor.arg = neighbor.arg, tasks = "discrete")
      par.list.refine[[jj]][subpart.index] <- p.part
      
    }#for jj
    
  }#for ii
  
  return(par.list.refine)
  
}#Niche.report.finetune




