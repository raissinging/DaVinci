#' Exhaustive decomposition
#'
#' Tiles the slice, decomposes on tile aggregates, then transfers back to pixel level
#' using impact_adaptive.
#' @export
exhaustive <- function(input.folder,        # path to folder containing {slice.id}.RDS files
                       output.dir,          # path to output folder; a subfolder per slice is created
                       k.arg      = 30,     # number of latent variables
                       L2_number  = 2000,   # target number of tiles per section
                       UMI.thr    = -Inf,   # minimum UMI per pixel; -Inf keeps all pixels
                       max.iter   = 0,      # iterative refinement passes; 0 = single-pass, no iterations
                       normalize  = T,      # run gene_normalization(); set F if data is pre-normalized
                       slice.id   = NULL){  # section ID (e.g. "43775"); NULL = process all .RDS in input.folder

  #if no slice.id, run all files in the input folder
  if (is.null(slice.id)){
    ids <- sub("\\.RDS$", "", list.files(input.folder, pattern = "\\.RDS$"))
    for (id in ids){
      exhaustive(input.folder = input.folder,
                 output.dir   = output.dir,
                 k.arg        = k.arg,
                 L2_number    = L2_number,
                 UMI.thr      = UMI.thr,
                 max.iter     = max.iter,
                 normalize    = normalize,
                 slice.id     = id)
    }#for id
    return(invisible(NULL))
  }#if

  L4.arg <- 50

  slice.out <- paste0(output.dir, slice.id, "/")
  dir.create(slice.out, recursive = T, showWarnings = F)

  #load
  input <- readRDS(paste0(input.folder, slice.id, ".RDS"))

  keep.index     <- which(Matrix::colSums(input$gene.exp) > UMI.thr)
  input$gene.exp <- input$gene.exp[, keep.index]
  input$coor     <- input$coor[keep.index, ]

  if (normalize){
    temp     <- gene_normalization(input$gene.exp, frac.thr = 0.95, MT.remove = T, median.norm = T)
    gene.exp <- temp$gene.exp
  }else{
    gene.exp <- input$gene.exp
  }#else
  coor <- input$coor[colnames(gene.exp), ]
  colnames(coor) <- c("array_row", "array_col")

  #tile the section
  grid_ids        <- tile_the_slice(coor, random.seed = 1, L2_number = L2_number)
  unique_grid_ids <- unique(grid_ids)
  tile_list       <- lapply(unique_grid_ids, function(x){ which(grid_ids == x) })

  tile_aggregation   <- lapply(tile_list, function(x){ apply(gene.exp[, names(x), drop = F], 1, mean) })
  tile_mat           <- do.call(cbind, tile_aggregation)
  colnames(tile_mat) <- paste0("Tile_", unique_grid_ids)
  names(tile_list)   <- paste0("Tile_", unique_grid_ids)

  coor_aggregation    <- lapply(tile_list, function(x){ apply(coor[names(x), , drop = F], 2, mean) })
  coor.tile           <- do.call(rbind, coor_aggregation)
  colnames(coor.tile) <- c("array_col", "array_row")
  rownames(coor.tile) <- paste0("Tile_", seq_along(coor.tile[, 1]))

  gene.exp.tile <- tile_mat
  L.tile        <- L_generate(coor.tile, opt = "Tri.mesh")$L

  #initial decomposition on tiles
  ICAp.res.tile0 <- manifoldDecomp_adaptive(gene.exp.tile, L.tile,
                                             k = k.arg, L4 = L4.arg, L4_adaptive = 2,
                                             to_drop = T, save.complete = T, verbose = F)

  #per-tile impact pass: transfer tile-level Z to pixel level
  density.grid  <- unlist(lapply(tile_list, length))
  gene.exp.list <- L.list <- L1.grid <- L2.grid <- shur0.grid <- B.list <- list()

  for (ii in 1:length(tile_list)){
    if (density.grid[ii] >= 3){
      gene.exp.inuse <- as.matrix(gene.exp[, names(tile_list[[ii]])])
      coor.inuse     <- coor[names(tile_list[[ii]]), , drop = F]
      temp           <- L_generate(coor.inuse, opt = "Tri.mesh")

      gene.exp.list[[ii]] <- gene.exp.inuse
      L.list[[ii]]        <- temp$L

      ICAp.res.inuse <- impact_adaptive(ICAp.res.tile0$Z, gene.exp.inuse,
                                         query.L  = temp$L, query.L4 = L4.arg,
                                         to_drop  = F, scale = 1,
                                         max.iter = 200, cor.thr = 0.8, verbose = F)
      L1.grid[[ii]]    <- ICAp.res.inuse$L1
      L2.grid[[ii]]    <- ICAp.res.inuse$L2
      shur0.grid[[ii]] <- ICAp.res.inuse$shur0
      B.list[[ii]]     <- ICAp.res.inuse$B
    }else{
      gene.exp.list[[ii]] <- L.list[[ii]] <- L1.grid[[ii]] <-
        L2.grid[[ii]] <- shur0.grid[[ii]] <- B.list[[ii]] <- NA
    }#else
  }#for ii

  #iterative refinement: repeat tile decomp + pixel transfer until convergence
  normF <- function(x){ sum(x^2) }

  B.agg.pre      <- lapply(B.list, function(x){ apply(x, 1, mean) })
  B.aggregate    <- do.call(cbind, B.agg.pre)
  error.relative <- normF(ICAp.res.tile0$B - B.aggregate) / normF(ICAp.res.tile0$B)
  error.accum    <- c(error.relative)
  count          <- 0
  ICAp.res.tile  <- ICAp.res.tile0

  while ((error.relative >= 1e-3) & (count < max.iter)){
    count <- count + 1
    message("Iteration ", count, " | error: ", error.relative)

    ICAp.res.tile <- manifoldDecomp_adaptive(gene.exp.tile, L.tile,
                                              B     = B.aggregate, k = k.arg, svdres = NULL,
                                              shur0 = ICAp.res.tile0$shur0,
                                              L1    = ICAp.res.tile0$L1, L2 = ICAp.res.tile0$L2,
                                              L4    = L4.arg, L4_adaptive = 2,
                                              to_drop = T, save.complete = T, verbose = F)

    B.list <- list()
    for (ii in 1:length(tile_list)){
      if (density.grid[ii] >= 3){
        ICAp.res.inuse <- impact_adaptive(ICAp.res.tile$Z, gene.exp.list[[ii]],
                                           query.L    = L.list[[ii]],
                                           query.L1   = L1.grid[[ii]], query.L2    = L2.grid[[ii]],
                                           query.L4   = L4.arg,        query.shur0 = shur0.grid[[ii]],
                                           to_drop    = F, scale = 1,
                                           max.iter   = 200, cor.thr = 0.8, verbose = F)
        B.list[[ii]] <- ICAp.res.inuse$B
      }else{
        B.list[[ii]] <- NA
      }#else
    }#for ii

    B.agg.pre      <- lapply(B.list, function(x){ apply(x, 1, mean) })
    B.aggregate    <- do.call(cbind, B.agg.pre)
    error.relative <- normF(ICAp.res.tile$B - B.aggregate) / normF(ICAp.res.tile$B)
    error.accum    <- c(error.accum, error.relative)
  }#while

  B.allspots <- do.call(cbind, B.list)
  Z.allspots <- ICAp.res.tile$Z

  tile.id    <- rep(names(tile_list), times = unlist(lapply(tile_list, length)))
  barcodes   <- unlist(lapply(tile_list, function(x){ names(x) }))
  section.id <- rep(slice.id, times = length(barcodes))

  tile.tab <- data.frame(section         = section.id,
                          tile            = tile.id,
                          barcodes        = barcodes,
                          barcodes.unique = paste0(section.id, "@", barcodes))

  out <- list(gene.exp.tile  = gene.exp.tile,
              L.tile         = L.tile,
              coor.tile      = coor.tile,
              ICAp.res.tile0 = ICAp.res.tile0,
              ICAp.res.tile  = ICAp.res.tile,
              coor           = coor,
              B.list         = B.list,
              B.allspots     = B.allspots,
              Z.allspots     = Z.allspots,
              error.accum    = error.accum,
              tile.tab       = tile.tab)

  saveRDS(out, file = paste0(slice.out, slice.id, "@k=", k.arg, ".RDS"))
  message("Saved: ", slice.out, slice.id, "@k=", k.arg, ".RDS")
  invisible(out)

}#exhaustive




#' Exhaustive decomposition using pre-integrated loadings
#'
#' Applies impact_adaptive tile-by-tile using a fixed Z matrix from a prior
#' cross-section integration step. Produces pixel-level B.allspots for clustering.
#' Use this after running exhaustive() + selfdeco + integration (steps 1-3).
#'
#' @export
exhaustive_integrated <- function(input.folder,              # path to folder containing {slice.id}.RDS files
                                  output.dir,                # path to output folder; a subfolder per slice is created
                                  model.dir,                 # path to integration output folder (contains .RData)
                                  model.name = "all@FineTune@first@default", # name of the .RData file (no extension)
                                  L2_number  = 2000,         # target number of tiles per section
                                  UMI.thr    = -Inf,         # minimum UMI per pixel; -Inf keeps all pixels
                                  slice.id   = NULL){        # section ID (e.g. "43775"); NULL = process all .RDS in input.folder

  #if no slice.id, run all files in the input folder
  if (is.null(slice.id)){
    ids <- sub("\\.RDS$", "", list.files(input.folder, pattern = "\\.RDS$"))
    for (id in ids){
      exhaustive_integrated(input.folder = input.folder,
                            output.dir   = output.dir,
                            model.dir    = model.dir,
                            model.name   = model.name,
                            L2_number    = L2_number,
                            UMI.thr      = UMI.thr,
                            slice.id     = id)
    }#for id
    return(invisible(NULL))
  }#if

  L4.arg <- 50

  slice.out <- paste0(output.dir, slice.id, "/")
  dir.create(slice.out, recursive = T, showWarnings = F)

  #load raw data
  input <- readRDS(paste0(input.folder, slice.id, ".RDS"))

  keep.index     <- which(Matrix::colSums(input$gene.exp) > UMI.thr)
  input$gene.exp <- input$gene.exp[, keep.index]
  input$coor     <- input$coor[keep.index, ]

  gene.exp <- input$gene.exp
  coor     <- input$coor

  #tile the section
  grid_ids        <- tile_the_slice(coor, random.seed = 1, L2_number = L2_number)
  unique_grid_ids <- unique(grid_ids)
  tile_list       <- lapply(unique_grid_ids, function(x){ which(grid_ids == x) })
  names(tile_list) <- paste0("Tile_", unique_grid_ids)

  #load integration model and extract Z
  load(paste0(model.dir, model.name, ".RData"), temp.env <- new.env())
  temp.env <- as.list(temp.env)
  Z.inuse  <- temp.env$integration.res$loading

  #convert Z list to matrix (genes x LVs)
  Z.genes    <- unique(unlist(lapply(Z.inuse, names)))
  Z.inuse.mat <- matrix(0, nrow = length(Z.genes), ncol = length(Z.inuse))
  rownames(Z.inuse.mat) <- Z.genes
  colnames(Z.inuse.mat) <- names(Z.inuse)
  for (ii in 1:length(Z.inuse)){
    Z.inuse.mat[names(Z.inuse[[ii]]), ii] <- Z.inuse[[ii]]
  }#for ii

  #per-tile impact pass with fixed Z
  density.grid  <- unlist(lapply(tile_list, length))
  gene.exp.list <- L.list <- L1.grid <- L2.grid <- shur0.grid <- B.list <- list()

  for (ii in 1:length(tile_list)){
    if (density.grid[ii] >= 3){
      gene.exp.inuse <- as.matrix(gene.exp[, names(tile_list[[ii]])])
      coor.inuse     <- coor[names(tile_list[[ii]]), , drop = F]
      temp           <- L_generate(coor.inuse, opt = "Tri.mesh")

      gene.exp.list[[ii]] <- gene.exp.inuse
      L.list[[ii]]        <- temp$L

      ICAp.res.inuse <- impact_adaptive(Z.inuse.mat, gene.exp.inuse,
                                         query.L  = temp$L, query.L4 = L4.arg,
                                         to_drop  = F, scale = 1,
                                         max.iter = 200, cor.thr = 0.8, verbose = F)
      L1.grid[[ii]]    <- ICAp.res.inuse$L1
      L2.grid[[ii]]    <- ICAp.res.inuse$L2
      shur0.grid[[ii]] <- ICAp.res.inuse$shur0
      B.list[[ii]]     <- ICAp.res.inuse$B
    }else{
      gene.exp.list[[ii]] <- L.list[[ii]] <- L1.grid[[ii]] <-
        L2.grid[[ii]] <- shur0.grid[[ii]] <- B.list[[ii]] <- NA
    }#else
  }#for ii

  B.allspots <- do.call(cbind, B.list)

  tile.id    <- rep(names(tile_list), times = unlist(lapply(tile_list, length)))
  barcodes   <- unlist(lapply(tile_list, function(x){ names(x) }))
  section.id <- rep(slice.id, times = length(barcodes))

  tile.tab <- data.frame(section         = section.id,
                          tile            = tile.id,
                          barcodes        = barcodes,
                          barcodes.unique = paste0(section.id, "@", barcodes))

  out <- list(B.list      = B.list,
              B.allspots  = B.allspots,
              Z.inuse.mat = Z.inuse.mat,
              tile_list   = tile_list,
              grid_ids    = grid_ids,
              tile.tab    = tile.tab)

  saveRDS(out, file = paste0(slice.out, slice.id, ".RDS"))
  message("Saved: ", slice.out, slice.id, ".RDS")
  invisible(out)

}#exhaustive_integrated




#' HD Self-decomposition
#'
#' Loads all per-k exhaustive outputs for a section, runs self_deco() to find
#' consensus latent variables across k values, then fine-tunes with manifoldDecomp_adaptive.
#' Use this after exhaustive() (step 1) and before integration (step 3).
#'
#' @export
selfdeco <- function(input.dir,                # path to exhaustive output folder (contains {slice.id}/ subfolders)
                     output.dir,               # path to output folder; a subfolder per slice is created
                     LVs.filter.thr = 0.9,    # correlation threshold for consensus LV filtering in self_deco
                     freq           = 1,       # minimum frequency (number of k values) an LV must appear in
                     slice.id       = NULL){   # section ID (e.g. "43775"); NULL = process all subfolders in input.dir

  #if no slice.id, run all subfolders in the input dir
  if (is.null(slice.id)){
    ids <- list.dirs(input.dir, full.names = F, recursive = F)
    for (id in ids){
      selfdeco(input.dir      = input.dir,
               output.dir     = output.dir,
               LVs.filter.thr = LVs.filter.thr,
               freq           = freq,
               slice.id       = id)
    }#for id
    return(invisible(NULL))
  }#if

  L4.arg <- 50

  model.folder <- paste0(input.dir, slice.id, "/")
  slice.out    <- paste0(output.dir, slice.id, "/")
  dir.create(slice.out, recursive = T, showWarnings = F)

  #auto-detect all k values from exhaustive output files
  ff         <- list.files(model.folder, pattern = "\\.RDS$")
  k.arg.list <- unlist(lapply(strsplit(ff, "@"), function(x){ x[2] }))
  k.arg.list <- sort(as.numeric(gsub(".RDS", "", gsub("k=", "", k.arg.list))))

  #load all per-k exhaustive results
  proj <- list()
  for (ii in 1:length(k.arg.list)){
    proj[[ii]] <- readRDS(paste0(model.folder, slice.id, "@k=", k.arg.list[ii], ".RDS"))
  }#for ii

  #find consensus LVs across k values
  embed <- self_deco(lapply(proj, function(x){ x$ICAp.res.tile }),
                     LVs.filter.thr = LVs.filter.thr,
                     freq           = freq,
                     opt            = "B")

  #fine-tune with consensus B as initialization
  mat <- proj[[1]]$gene.exp.tile
  L   <- proj[[1]]$L.tile

  ICAp.res <- manifoldDecomp_adaptive(mat, L,
                                       k           = nrow(embed$LVs),
                                       B           = embed$LVs,
                                       L4          = L4.arg,
                                       L4_adaptive = 2,
                                       to_drop     = T,
                                       save.complete = T,
                                       verbose     = F)

  save(embed, ICAp.res, file = paste0(slice.out, slice.id, ".RData"))
  message("Saved: ", slice.out, slice.id, ".RData")
  invisible(list(embed = embed, ICAp.res = ICAp.res))

}#selfdeco




#' HD Horizontal Integration
#'
#' Loads selfdeco outputs for all sections and runs horizontal integration across them.
#' Produces the integrated loading matrix Z used by exhaustive_integrated() (step 5).
#' Use this after selfdeco() (step 2).
#'
#' @export
integration <- function(input.dir,                    # path to selfdeco output folder (contains {slice.id}/ subfolders)
                        output.dir,                   # path to output folder; tile/ subfolder is created inside
                        mod.opt        = "all",       # gene modality mode: "all" (all genes) or "common" (shared genes only)
                        input.opt      = "FineTune",  # which decomp to use: "FineTune" or "OnlyDeco"
                        h.opt          = "first",     # integration method: "first" or "default"
                        L2.in          = "default",   # L2 normalization option: "default", "L2norm", or "L2norm.joint"
                        LVs.filter.thr = 0.8){        # correlation threshold for LV filtering during integration

  tile.out <- paste0(output.dir, "tile/")
  dir.create(tile.out, recursive = T, showWarnings = F)

  #detect all sections from selfdeco output subfolders
  dataset.opts <- list.dirs(input.dir, full.names = F, recursive = F)

  #load selfdeco results for each section
  Y.list       <- list()
  L.list       <- list()
  dav.res.list <- list()

  for (ii in 1:length(dataset.opts)){
    input.file <- paste0(input.dir, dataset.opts[ii], "/", dataset.opts[ii], ".RData")
    load(input.file, verbose = F)

    Y.list[[ii]] <- ICAp.res$Y
    L.list[[ii]] <- ICAp.res$L

    LV.var <- VarianceExplained(ICAp.res$Y, ICAp.res$Z, ICAp.res$B,
                                 option = "simple", normalize = F)

    if (input.opt == "OnlyDeco"){
      tmp      <- list(Z = t(embed$LVs.pair), B = embed$LVs,
                       L4 = ICAp.res$L4, shur0 = ICAp.res$shur0,
                       L1 = ICAp.res$L1,  L2 = ICAp.res$L2)
    }else if (input.opt == "FineTune"){
      tmp <- ICAp.res
    }#else if

    dav.res.list[[ii]] <- tmp
  }#for ii

  #horizontal integration
  if (h.opt == "default"){
    integration.res <- Horizontal.Integration(Y.list, L.list,
                                               dav.res.list   = dav.res.list,
                                               LVs.filter.thr = LVs.filter.thr,
                                               mod            = mod.opt,
                                               remove.LV1     = F,
                                               L2.option      = L2.in,
                                               batch.correction = "harmony")
  }else if (h.opt == "first"){
    integration.res <- Horizontal.Integration.first(Y.list, L.list,
                                                     dav.res.list   = dav.res.list,
                                                     LVs.filter.thr = LVs.filter.thr,
                                                     mod            = mod.opt,
                                                     remove.LV1     = F,
                                                     L2.option      = L2.in,
                                                     batch.correction = "harmony")
  }#else if

  #align names
  embed        <- integration.res$LVs_embeddings
  sids         <- unlist(lapply(strsplit(rownames(embed), "_"), function(x){ paste0(x[-1], collapse = "_") }))
  rownames(embed) <- sids

  out.name <- paste0(tile.out, mod.opt, "@", input.opt, "@", h.opt, "@", L2.in, ".RData")
  save(embed, integration.res, file = out.name)
  message("Saved: ", out.name)
  invisible(list(embed = embed, integration.res = integration.res))

}#integration




#' HD Tile-level clustering
#'
#' Clusters the integrated tile-level embeddings (output of integration()) using
#' Louvain + SNN with adaptive resolution, then applies spatial smoothing per section.
#' Sweeps over k and num.of.clusters; skips combinations that already exist.
#'
#' @export
cluster_tile <- function(integration.dir,                           # path to integration tile/ output folder (contains .RData)
                         exhaustive.dir,                            # path to exhaustive output (for coor.tile per section)
                         output.dir,                                # path to output folder; louvain/ subfolder is created
                         model.name      = "all@FineTune@first@default", # .RData filename (no extension)
                         mn              = "DaVinci",              # model name prefix for output files
                         modality        = "RNA",                   # modality label for output files
                         k.opt.list      = c(30, 35, 40, 45, 50), # SNN k values to sweep
                         num.of.clusters.opts = c(10, 15, 30, 50, 100), # cluster counts to sweep
                         neighbor.arg    = 8){                      # KNN k for spatial smoothing

  louvain.out <- paste0(output.dir, "louvain/")
  dir.create(louvain.out, recursive = T, showWarnings = F)

  #detect sections from exhaustive output
  dataset.list <- list.dirs(exhaustive.dir, full.names = F, recursive = F)

  #load tile coordinates per section
  coor.list <- list()
  for (ii in 1:length(dataset.list)){
    ff    <- list.files(paste0(exhaustive.dir, dataset.list[ii]))
    input <- readRDS(paste0(exhaustive.dir, dataset.list[ii], "/", ff[1]))
    coor.list[[dataset.list[ii]]] <- input$coor.tile
  }#for ii

  #load integration embeddings
  load(paste0(integration.dir, model.name, ".RData"), temp.env <- new.env())
  temp.env <- as.list(temp.env)

  embed        <- temp.env$integration.res$LVs_embeddings
  sids         <- unlist(lapply(strsplit(rownames(embed), "_"), function(x){ paste0(x[-1], collapse = "_") }))
  mat.slice.id <- unlist(lapply(strsplit(rownames(embed), "_"), function(x){ x[1] }))
  mat.slice.id <- dataset.list[as.numeric(mat.slice.id)]
  rownames(embed) <- sids
  mat <- as.matrix(embed)

  mat.norm <- L2Norm(mat, MARGIN = 2)

  for (num.of.clusters in num.of.clusters.opts){
    for (k.opt in k.opt.list){

      out.file <- paste0(louvain.out, mn, "@", modality, "@k=", k.opt, "@n=", num.of.clusters, ".RDS")
      if (file.exists(out.file)){
        cat("Skipping (exists):", basename(out.file), "\n")
        next
      }#if

      set.seed(1)
      snn.res <- Seurat::FindNeighbors(mat.norm, k.param = k.opt,
                                        return.neighbor = F, compute.SNN = T, verbose = F)

      p.res <- leiden_adaptive(snn.res$snn,
                                num.of.cluster   = num.of.clusters,
                                resolution.start = 0.4,
                                adaptive.size    = 2,
                                method           = "louvain",
                                full             = T,
                                remove.singleton = mat.norm)

      partition        <- as.character(p.res$partition)
      names(partition) <- rownames(mat)

      #spatial smoothing per section
      partition.smooth <- partition
      for (ii in 1:length(dataset.list)){
        subpart.index <- which(mat.slice.id == dataset.list[ii])
        p.part <- partition[subpart.index]
        p.part <- refinement(p.part,
                              as.matrix(coor.list[[dataset.list[ii]]])[names(p.part), ],
                              neighbor.option = "KNN",
                              neighbor.arg    = neighbor.arg,
                              tasks           = "discrete")
        partition.smooth[subpart.index] <- p.part
      }#for ii

      saveRDS(list(dataset.list     = dataset.list,
                    num.of.clusters  = num.of.clusters,
                    k.opt            = k.opt,
                    partition        = partition,
                    partition.smooth = partition.smooth,
                    mat.slice.id     = mat.slice.id),
               file = out.file)
      cat("Saved:", basename(out.file), "\n")

    }#for k.opt
  }#for num.of.clusters

  invisible(NULL)

}#cluster_tile




#' HD Pixel-level clustering with sketch
#'
#' Loads pixel-level B.allspots from exhaustive_integrated(), applies Harmony
#' batch correction, then clusters via leverage-weighted sketch + SNN + Louvain
#' on the sketch only, transfers labels to all pixels via weighted kNN, and applies
#' multi-pass spatial smoothing. Sweeps over k and num.of.clusters.
#'
#' @export
cluster_pixel <- function(input.dir,                                    # path to exhaustive_integrated output (contains {slice.id}/ subfolders)
                          raw.dir,                                       # path to raw RDS files (for pixel coordinates)
                          output.dir,                                    # path to output folder
                          dataset.list,                                  # character vector of section IDs to cluster jointly
                          mn                   = "DaVinci",             # model name prefix for output files
                          modality             = "Lipids",               # modality label for output files
                          k.opt.list           = c(30, 35, 40, 45, 50),  # SNN k values to sweep
                          num.of.clusters.opts = c(10, 15, 30, 50, 100), # cluster counts to sweep
                          harmony              = T,                      # apply Harmony batch correction across sections
                          sketch.per.section   = 10000,                  # target sketch pixels per section
                          sketch.method        = "leverage",             # sampling: "leverage" (oversample rare) or "uniform"
                          lev.winsor           = 0.99,                  # winsorize leverage at this quantile; 1 = no cap
                          k.transfer           = 25,                    # kNN neighbors used to transfer sketch labels to all pixels
                          chunk.size           = 200000,                # pixels per transfer chunk (memory control)
                          weighted.vote        = T,                     # adaptive Gaussian weighting for kNN vote; F = majority
                          smooth.neighbor      = 8,                     # KNN k for spatial smoothing
                          smooth.iters         = 1){                    # number of spatial smoothing passes

  dir.create(output.dir, recursive = T, showWarnings = F)

  #internal: per-section leverage scores (hat-matrix diagonal on Harmony embedding)
  .leverage_section <- function(M){
    Mc  <- scale(M, center = T, scale = F)
    G   <- crossprod(Mc)
    G   <- G + diag(1e-8 * (mean(diag(G)) + 1e-12), ncol(G))
    lev <- tryCatch(rowSums((Mc %*% solve(G)) * Mc),
                    error = function(e) rep(1, nrow(M)))
    lev[!is.finite(lev) | lev < 0] <- 0
    lev
  }#.leverage_section

  #internal: vectorised weighted-kNN vote over one chunk
  .vote_chunk <- function(nn.index, nn.dist, sketch.code, classes, weighted){
    C     <- length(classes)
    codes <- matrix(sketch.code[nn.index], nrow = nrow(nn.index), ncol = ncol(nn.index))
    if (weighted){
      bw <- nn.dist[, ncol(nn.dist)]
      bw[!is.finite(bw) | bw <= 0] <- 1
      w  <- exp(-(nn.dist^2) / (2 * bw^2))
    }else{
      w <- matrix(1, nrow = nrow(nn.index), ncol = ncol(nn.index))
    }#else
    wsum   <- matrix(0, nrow = nrow(codes), ncol = C)
    for (cc in 1:C) wsum[, cc] <- rowSums(w * (codes == cc))
    winner <- max.col(wsum, ties.method = "first")
    tot    <- rowSums(wsum); tot[tot == 0] <- 1
    conf   <- wsum[cbind(1:nrow(wsum), winner)] / tot
    list(label = classes[winner], conf = conf)
  }#.vote_chunk

  #load pixel-level B.allspots for each section
  mat.list     <- list()
  slice.id.vec <- c()

  for (sid in dataset.list){
    f <- paste0(input.dir, sid, "/", sid, ".RDS")
    if (!file.exists(f)){ cat("Missing:", f, "\n"); next }
    res <- readRDS(f)
    B   <- t(res$B.allspots)   #pixels x LVs
    rm(res); gc()
    mat.list[[sid]] <- B
    slice.id.vec    <- c(slice.id.vec, rep(sid, nrow(B)))
  }#for sid

  if (length(mat.list) == 0) stop("No data loaded.")

  mat          <- do.call(rbind, mat.list)
  rm(mat.list); gc()
  mat.slice.id <- slice.id.vec
  names(mat.slice.id) <- rownames(mat)
  px.names <- rownames(mat)

  #load pixel coordinates for spatial smoothing
  coor.list <- list()
  for (sid in dataset.list){
    f <- paste0(raw.dir, sid, ".RDS")
    if (!file.exists(f)) next
    input <- readRDS(f)
    coor.list[[sid]] <- input$coor
  }#for sid

  mat.norm <- L2Norm(mat, MARGIN = 2)

  #Harmony batch correction
  if (harmony){
    message("Running Harmony batch correction...")
    meta        <- data.frame(section = mat.slice.id)
    mat.harmony <- harmony::HarmonyMatrix(mat.norm, meta, vars_use = "section",
                                          do_pca = FALSE, verbose = T)
    mat.harmony <- as.matrix(mat.harmony)
    rownames(mat.harmony) <- px.names
    message("Harmony done.")
  }else{
    mat.harmony <- mat.norm
  }#else

  rm(mat, mat.norm); gc()
  N <- nrow(mat.harmony)

  #sketch: per-section leverage-weighted (or uniform) sampling, fixed across all param sweeps
  set.seed(1)
  sketch.idx <- unlist(lapply(dataset.list, function(sid){
    ix     <- which(mat.slice.id == sid)
    if (length(ix) == 0) return(integer(0))
    n.take <- min(length(ix), sketch.per.section)
    if (length(ix) <= n.take) return(ix)
    prob <- NULL
    if (sketch.method == "leverage"){
      prob <- .leverage_section(mat.harmony[ix, , drop = F])
      if (lev.winsor < 1){
        cap <- stats::quantile(prob[prob > 0], lev.winsor, names = F)
        if (is.finite(cap) && cap > 0) prob <- pmin(prob, cap)
      }#if
      if (!any(prob > 0)) prob <- NULL
    }#if
    sample(ix, size = n.take, prob = prob)
  }), use.names = F)
  sketch.idx <- sort(unique(sketch.idx))
  cat("Sketch:", length(sketch.idx), "of", N, "pixels | method:", sketch.method, "\n")

  mat.sketch <- mat.harmony[sketch.idx, , drop = F]

  for (num.of.clusters in num.of.clusters.opts){
    for (k.opt in k.opt.list){

      out.file <- paste0(output.dir, mn, "@", modality, "@k=", k.opt, "@n=", num.of.clusters, ".RDS")
      if (file.exists(out.file)){
        cat("Skipping (exists):", basename(out.file), "\n")
        next
      }#if

      cat("\n--- k =", k.opt, "| n clusters =", num.of.clusters, "---\n")

      #SNN + Louvain on sketch only
      set.seed(1)
      snn.res <- Seurat::FindNeighbors(mat.sketch, k.param = k.opt,
                                        return.neighbor = F, compute.SNN = T, verbose = F)
      p.res <- leiden_adaptive(snn.res$snn,
                                num.of.cluster   = num.of.clusters,
                                resolution.start = 0.4,
                                adaptive.size    = 2,
                                method           = "louvain",
                                full             = T,
                                remove.singleton = mat.sketch)

      partition.sketch        <- as.character(p.res$partition)
      names(partition.sketch) <- px.names[sketch.idx]
      classes     <- sort(unique(partition.sketch))
      sketch.code <- match(partition.sketch, classes)
      rm(snn.res, p.res); gc()

      #transfer labels to all pixels via chunked weighted kNN
      partition      <- character(N); names(partition)      <- px.names
      transfer.score <- numeric(N);   names(transfer.score) <- px.names
      k.use <- min(k.transfer, nrow(mat.sketch))

      for (st in seq(1, N, by = chunk.size)){
        en <- min(st + chunk.size - 1, N)
        nn <- FNN::get.knnx(data  = mat.sketch,
                             query = mat.harmony[st:en, , drop = F],
                             k     = k.use)
        v  <- .vote_chunk(nn$nn.index, nn$nn.dist, sketch.code, classes, weighted.vote)
        partition[st:en]      <- v$label
        transfer.score[st:en] <- v$conf
        rm(nn, v); gc()
      }#for st

      #multi-pass spatial smoothing per section
      partition.smooth <- partition
      for (sid in dataset.list){
        subpart.index <- which(mat.slice.id == sid)
        if (length(subpart.index) == 0 || is.null(coor.list[[sid]])) next
        p.part   <- partition[subpart.index]
        coor.sid <- as.matrix(coor.list[[sid]])[names(p.part), ]
        for (iter in 1:smooth.iters){
          p.part <- refinement(p.part, coor.sid,
                                neighbor.option = "KNN",
                                neighbor.arg    = smooth.neighbor,
                                tasks           = "discrete")
        }#for iter
        partition.smooth[subpart.index] <- p.part
      }#for sid

      saveRDS(list(dataset.list       = dataset.list,
                    num.of.clusters    = num.of.clusters,
                    k.opt              = k.opt,
                    harmony            = harmony,
                    sketch.method      = sketch.method,
                    sketch.per.section = sketch.per.section,
                    k.transfer         = k.use,
                    weighted.vote      = weighted.vote,
                    smooth.neighbor    = smooth.neighbor,
                    smooth.iters       = smooth.iters,
                    sketch.idx         = sketch.idx,
                    partition.sketch   = partition.sketch,
                    partition          = partition,
                    partition.smooth   = partition.smooth,
                    transfer.score     = transfer.score,
                    mat.slice.id       = mat.slice.id),
               file = out.file)
      cat("Saved:", basename(out.file), "\n")

      rm(partition, partition.smooth, transfer.score,
         partition.sketch, sketch.code, classes); gc()

    }#for k.opt
  }#for num.of.clusters

  invisible(NULL)

}#cluster_pixel
