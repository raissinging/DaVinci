#' Find top latent variables per cluster
#'
#' Ranks latent variables (LVs) that best discriminate each cluster from the rest,
#' using elastic net logistic regression and/or naive mean-difference, run per section
#' and averaged. Works with tile-level clusters (from cluster_tile) or pixel-level
#' clusters (from cluster_pixel).
#'
#' @export
find_top_lvs <- function(exhaustive.dir,                # path to exhaustive output folder
                         integration.RData,             # path to integration .RData file
                         clusters.rds,                  # path to cluster_tile/cluster_pixel output .RDS
                         output.dir,                    # path to output folder
                         embed.dir        = NULL,       # path to exhaustive_integrated output for pixel-level embed; NULL = use tile-level integration embed
                         elastic.net      = T,          # run elastic net LV ranking
                         naive            = F,          # run naive mean-diff LV ranking
                         downsample       = NULL,       # target pos:neg ratio; NULL = no downsampling
                         top.n            = -1,         # top N LVs to report; -1 = all non-zero
                         alpha.sweep      = c(0.5),     # alpha values to sweep; best chosen by CV error or BIC
                         cv.folds         = NULL,       # number of CV folds for lambda selection; NULL = BIC (no CV)
                         seed             = 42,         # random seed
                         verbose          = F,          # print progress messages
                         clusters.subset  = NULL){      # restrict to these cluster labels; NULL = all clusters

  set.seed(seed)

  if (verbose) message("Loading data...")

  dataset.list <- list.files(exhaustive.dir)

  if (!is.null(embed.dir)){
    # pixel-level: load B.allspots from each exhaustive_integrated section RDS
    embed.blocks <- lapply(dataset.list, function(sid){
      ff <- list.files(paste0(embed.dir, sid), pattern = "\\.RDS$", full.names = TRUE)
      if (length(ff) == 0) return(NULL)
      obj <- readRDS(ff[1])
      B   <- obj$B.allspots          # LV x pixel
      if (is.null(B)) return(NULL)
      mat.px <- t(B)                 # pixel x LV
      rownames(mat.px) <- paste0(sid, "@", colnames(B))
      mat.px
    })
    embed.blocks <- Filter(Negate(is.null), embed.blocks)
    mat <- do.call(rbind, embed.blocks)
    if (verbose) cat("Pixel-level embed:", nrow(mat), "pixels x", ncol(mat), "LVs\n")
  } else {
    # tile-level: load from integration .RData
    load(integration.RData, temp.env1 <- new.env())
    temp.env1 <- as.list(temp.env1)
    embed    <- temp.env1$integration.res$LVs_embeddings
    full.ids <- rownames(embed)
    sec.idx  <- as.numeric(unlist(lapply(strsplit(full.ids, "_"), function(x) x[1])))
    tile.ids <- unlist(lapply(strsplit(full.ids, "_"), function(x) paste0(x[-1], collapse = "_")))
    sec.names <- dataset.list[sec.idx]
    rownames(embed) <- paste0(sec.names, "@", tile.ids)
    mat <- as.matrix(embed)
  }

  res              <- readRDS(clusters.rds)
  clusters         <- res$partition.smooth
  slice.ids        <- res$mat.slice.id
  names(clusters)  <- paste0(slice.ids, "@", names(clusters))

  common <- intersect(rownames(mat), names(clusters))
  if (verbose) cat("Common entries:", length(common), "\n")

  mat.sub <- mat[common, , drop = F]
  cl.sub  <- clusters[common]
  colnames(mat.sub) <- paste0("LV", 1:ncol(mat.sub))
  mat.sub <- as.matrix(scale(mat.sub))

  sec.sub <- sub("@.*", "", rownames(mat.sub))

  dir.create(output.dir, recursive = T, showWarnings = F)
  clusters.uniq <- sort(unique(cl.sub))
  if (!is.null(clusters.subset))
    clusters.uniq <- clusters.uniq[as.character(clusters.uniq) %in% as.character(clusters.subset)]

  lv.names <- colnames(mat.sub)
  sections  <- unique(sec.sub)

  #naive: per section mean-diff, aggregated across sections
  if (naive){
    if (verbose) message("Running naive LV ranking (per section)...")

    naive.rows <- lapply(clusters.uniq, function(cl){
      if (verbose) message(sprintf("  Cluster %s... ", cl), appendLF = F)

      sec.diffs <- lapply(sections, function(sec){
        sec.idx.s <- which(sec.sub == sec)
        cl.mask   <- cl.sub[sec.idx.s] == cl
        if (sum(cl.mask) < 2 || sum(!cl.mask) < 2) return(NULL)
        colMeans(mat.sub[sec.idx.s[cl.mask],  , drop = F]) -
        colMeans(mat.sub[sec.idx.s[!cl.mask], , drop = F])
      })
      sec.diffs <- sec.diffs[!sapply(sec.diffs, is.null)]
      if (length(sec.diffs) == 0) return(NULL)

      mean.diff <- rowMeans(do.call(cbind, sec.diffs))
      n.sec     <- length(sec.diffs)

      use.top.n <- if (top.n == -1) length(lv.names) else top.n
      top.idx   <- order(abs(mean.diff), decreasing = T)[1:use.top.n]
      top.lvs   <- lv.names[top.idx]

      if (verbose) message(sprintf("done (%d sections)", n.sec))

      data.frame(cluster    = cl,
                 rank       = 1:use.top.n,
                 LV         = top.lvs,
                 mean_diff  = mean.diff[top.idx],
                 n_sections = n.sec,
                 stringsAsFactors = F)
    })
    naive.rows <- naive.rows[!sapply(naive.rows, is.null)]
    naive.df   <- do.call(rbind, naive.rows)
    naive.df$cluster <- tryCatch(as.integer(naive.df$cluster),
      warning = function(w) tryCatch(as.numeric(naive.df$cluster),
        warning = function(w) naive.df$cluster))
    naive.df <- naive.df[order(naive.df$cluster, naive.df$rank), ]
    write.table(naive.df,
                file.path(output.dir, "naive_topLVs.tsv"),
                sep = "\t", quote = F, row.names = F)
    cat("Wrote naive_topLVs.tsv\n")
  }#if naive

  #elastic net: per section, aggregate mean coef_abs across sections
  if (elastic.net){
    if (verbose) message("Running elastic net LV ranking (per section)...")

    en.rows <- list()

    for (cl in clusters.uniq){
      if (verbose) message(sprintf("  Cluster %s... ", cl), appendLF = F)

      sec.coefs <- list()

      for (sec in sections){
        sec.idx.s <- which(sec.sub == sec)
        y.sec     <- as.integer(cl.sub[sec.idx.s] == cl)
        pos.s     <- which(y.sec == 1)
        neg.s     <- which(y.sec == 0)

        if (length(pos.s) < 5 || length(neg.s) < 5) next

        eq.ratio <- length(pos.s) / length(neg.s)
        if (!is.null(downsample)){
          n.neg.s <- min(length(neg.s), round(length(pos.s) / downsample))
          neg.s   <- sample(neg.s, n.neg.s)
        }#if

        all.idx <- c(pos.s, neg.s)
        X.all   <- mat.sub[sec.idx.s[all.idx], ]; y.all <- y.sec[all.idx]

        if (is.null(cv.folds)){
          fit.list <- lapply(alpha.sweep, function(a)
            glmnet::glmnet(X.all, y.all, family = "binomial", alpha = a))
          bic.min  <- sapply(fit.list, function(fit){
            min((1 - fit$dev.ratio) * fit$nulldev + log(nrow(X.all)) * fit$df)
          })
          best.i   <- which.min(bic.min)
          fit.best <- fit.list[[best.i]]
          bic      <- (1 - fit.best$dev.ratio) * fit.best$nulldev + log(nrow(X.all)) * fit.best$df
          best.l   <- fit.best$lambda[which.min(bic)]
          coefs    <- coef(fit.best, s = best.l)[-1]
        }else{
          cv.list  <- lapply(alpha.sweep, function(a)
            glmnet::cv.glmnet(X.all, y.all, family = "binomial", alpha = a, nfolds = cv.folds))
          best.i   <- which.min(sapply(cv.list, function(f) min(f$cvm)))
          cv.fit   <- cv.list[[best.i]]
          coefs    <- coef(cv.fit, s = "lambda.min")[-1]
        }#else

        names(coefs)     <- lv.names
        sec.coefs[[sec]] <- abs(coefs)
      }#for sec

      n.sec <- length(sec.coefs)
      if (n.sec == 0){
        if (verbose) message("skipped (no sections with enough spots)")
        next
      }#if

      coef.mat      <- do.call(cbind, sec.coefs)
      mean.coef.abs <- rowMeans(coef.mat)
      n.nonzero     <- rowSums(coef.mat > 0)

      if (top.n == -1){
        keep    <- which(mean.coef.abs > 0)
        top.idx <- keep[order(mean.coef.abs[keep], decreasing = T)]
      }else{
        top.idx <- order(mean.coef.abs, decreasing = T)[1:top.n]
      }#else
      top.lvs <- lv.names[top.idx]
      use.n   <- length(top.idx)

      en.rows[[cl]] <- data.frame(
        cluster            = cl,
        rank               = 1:use.n,
        LV                 = top.lvs,
        mean_coef_abs      = mean.coef.abs[top.idx],
        n_sections_nonzero = n.nonzero[top.idx],
        n_sections_total   = n.sec,
        stringsAsFactors   = F
      )

      if (verbose) message(sprintf("done (%d sections)", n.sec))
    }#for cl

    en.df <- do.call(rbind, en.rows)
    if (is.null(en.df) || nrow(en.df) == 0){
      cat("No results to write for elastic net.\n")
    }else{
      en.df$cluster <- tryCatch(as.integer(en.df$cluster),
        warning = function(w) tryCatch(as.numeric(en.df$cluster),
          warning = function(w) en.df$cluster))
      en.df <- en.df[order(en.df$cluster, en.df$rank), ]
      write.table(en.df,
                  file.path(output.dir, "elasticnet_topLVs.tsv"),
                  sep = "\t", quote = F, row.names = F)
      cat("Wrote elasticnet_topLVs.tsv\n")
    }#if not empty
  }#if elastic.net

  invisible(NULL)

}#find_top_lvs


#' Evaluate cluster quality across n-cluster solutions
#'
#' Loads pre-computed clusterings and computes FM index and Silhouette
#' score across cluster counts, then plots the curves and reports the optimal n.
#' FM is computed on the full cluster labels; Silhouette may use a subsampled LV
#' embedding (controlled by subsample).
#'
#' @export
n_cluster_metrics <- function(exhaustive.dir,                       # path to exhaustive output
                           integration.RData,                    # path to integration .RData file
                           cluster.dir,                          # exact path to folder containing the .RDS files
                           embed.dir   = NULL,                   # path to exhaustive_integrated output; required when pixel.level = TRUE
                           pixel.level = FALSE,                  # TRUE = pixel-level clusters; FALSE = tile-level
                           k.opt       = 40,                     # which k to select from cluster files in the folder
                           metric      = "both",                 # "silhouette", "fm", or "both"
                           subsample   = 10000,                  # max tiles/pixels for silhouette; NULL = use all
                           seed        = 42,                     # random seed for subsampling
                           title       = NULL,                   # plot title; NULL = auto
                           verbose     = TRUE){                  # print per-n progress

  metric <- match.arg(metric, c("silhouette", "fm", "both"))

  set.seed(seed)

  files.tmp <- list.files(cluster.dir, pattern = paste0("k=", k.opt), full.names = TRUE)

  if (length(files.tmp) == 0)
    stop("No files found in cluster.dir matching k=", k.opt)

  if (pixel.level) {
    if (is.null(embed.dir))
      stop("embed.dir must be provided when pixel.level = TRUE")
    dataset.list <- readRDS(files.tmp[1])$dataset.list
    embed.blocks <- lapply(dataset.list, function(sid) {
      ff <- list.files(file.path(embed.dir, sid), pattern = "\\.RDS$", full.names = TRUE)
      if (!length(ff)) return(NULL)
      B <- readRDS(ff[1])$B.allspots
      if (is.null(B)) return(NULL)
      mat.px <- t(B)
      if (!is.null(colnames(B)))
        rownames(mat.px) <- paste0(sid, "@", colnames(B))
      mat.px
    })
    embed.blocks <- Filter(Negate(is.null), embed.blocks)
    mat <- L2Norm(do.call(rbind, embed.blocks), MARGIN = 1)
  } else {
    dataset.list <- readRDS(files.tmp[1])$dataset.list

    load(integration.RData, e <- new.env())
    e <- as.list(e)

    embed <- e$integration.res$LVs_embeddings

    full.ids    <- rownames(embed)
    sec.idx.num <- as.numeric(unlist(lapply(strsplit(full.ids, "_"), function(x) x[1])))
    tile.ids    <- unlist(lapply(strsplit(full.ids, "_"), function(x) paste0(x[-1], collapse = "_")))

    sec.names <- dataset.list[sec.idx.num]

    rownames(embed) <- paste0(sec.names, "@", tile.ids)

    mat <- L2Norm(as.matrix(embed), MARGIN = 1)
  }#if pixel.level

  if (!is.null(subsample) && nrow(mat) > subsample){
    mat.sub <- mat[sample(nrow(mat), subsample), , drop = FALSE]
  } else {
    mat.sub <- mat
  }

  .fm_index <- function(a, b){
    tab <- table(a, b)

    TP <- sum(choose(tab, 2))
    FP <- sum(choose(rowSums(tab), 2)) - TP
    FN <- sum(choose(colSums(tab), 2)) - TP

    denom <- sqrt((TP + FP) * (TP + FN))

    if (denom == 0)
      return(NA_real_)

    TP / denom
  }#.fm_index

  files <- list.files(cluster.dir, pattern = paste0("k=", k.opt), full.names = TRUE)

  if (length(files) == 0)
    stop("No files found in cluster.dir matching k=", k.opt)

  n.vals <- as.numeric(sub(".*n=([0-9]+)\\.RDS", "\\1", basename(files)))

  files  <- files[order(n.vals)]
  n.vals <- sort(n.vals)

  if (verbose)
    cat("Loading", length(files), "clusterings (k=", k.opt, ")...\n")

  clusterings <- lapply(files, function(f) {
    r  <- readRDS(f)
    cl <- r$partition.smooth

    names(cl) <- paste0(r$mat.slice.id, "@", names(cl))

    cl
  })

  names(clusterings) <- as.character(n.vals)

  results <- do.call(rbind, lapply(seq_along(n.vals), function(i){

    cl.full <- clusterings[[i]]
    fm      <- NA_real_

    if (metric %in% c("fm", "both")){

      fmi.vals <- numeric(0)

      if (i > 1){
        common2 <- intersect(names(cl.full), names(clusterings[[i - 1]]))

        if (length(common2) > 0){
          fmi.vals <- c(
            fmi.vals,
            .fm_index(cl.full[common2], clusterings[[i - 1]][common2])
          )
        }
      }

      if (i < length(n.vals)){
        common2 <- intersect(names(cl.full), names(clusterings[[i + 1]]))

        if (length(common2) > 0){
          fmi.vals <- c(
            fmi.vals,
            .fm_index(cl.full[common2], clusterings[[i + 1]][common2])
          )
        }
      }

      fm <- if (length(fmi.vals) > 0) mean(fmi.vals, na.rm = TRUE) else NA_real_
    }#if fm

    sil <- NA_real_

    if (metric %in% c("silhouette", "both")){

      common <- intersect(rownames(mat.sub), names(cl.full))

      if (length(common) >= 2){

        cl <- cl.full[common]
        sm <- mat.sub[common, , drop = FALSE]

        cl.int <- as.integer(as.factor(cl))

        if (length(unique(cl.int)) >= 2 && length(unique(cl.int)) < length(cl.int)){
          sil <- mean(cluster::silhouette(cl.int, stats::dist(sm))[, 3])
        }
      }
    }#if silhouette

    if (verbose){

      msg <- sprintf("n=%d", n.vals[i])

      if (metric %in% c("fm", "both"))
        msg <- paste0(msg, sprintf("  FM=%.3f", fm))

      if (metric %in% c("silhouette", "both"))
        msg <- paste0(msg, sprintf("  sil=%.3f", sil))

      cat(msg, "\n")
    }#if verbose

    data.frame(
      n          = n.vals[i],
      FM         = fm,
      Silhouette = sil
    )
  }))

  best.sil <- if (
    metric %in% c("silhouette", "both") &&
      any(!is.na(results$Silhouette))
  ) {
    results$n[which.max(ifelse(is.na(results$Silhouette), -Inf, results$Silhouette))]
  } else {
    NULL
  }

  best.fm <- if (
    metric %in% c("fm", "both") &&
      any(!is.na(results$FM))
  ) {
    results$n[which.max(ifelse(is.na(results$FM), -Inf, results$FM))]
  } else {
    NULL
  }

  if (verbose){

    parts <- c()

    if (!is.null(best.sil))
      parts <- c(parts, paste0("Silhouette: ", best.sil))

    if (!is.null(best.fm))
      parts <- c(parts, paste0("FM: ", best.fm))

    if (length(parts) > 0)
      cat(sprintf("\nCandidate n — %s\n", paste(parts, collapse = " | ")))
  }#if verbose

  plot.title <- if (!is.null(title)) {
    title
  } else {
    suf <- if (metric %in% c("silhouette", "both") && !is.null(subsample)) {
      paste0(" (sil subsample=", subsample, ")")
    } else {
      ""
    }

    paste0("Cluster quality — k=", k.opt, suf)
  }

  active.colors <- c(
    "FM"         = "purple",
    "Silhouette" = "forestgreen"
  )

  if (metric == "both"){

    fm.rng  <- range(results$FM, na.rm = TRUE)
    sil.rng <- range(results$Silhouette, na.rm = TRUE)

    if (!all(is.finite(fm.rng)) || !all(is.finite(sil.rng))){
      warning("Cannot make dual-axis plot because FM or Silhouette has no finite values.")
      print(results)
      return(invisible(results))
    }

    if (diff(fm.rng) == 0)
      fm.rng <- fm.rng + c(-0.01, 0.01)

    if (diff(sil.rng) == 0)
      sil.rng <- sil.rng + c(-0.01, 0.01)

    .sil_to_fm <- function(x){
      (x - sil.rng[1]) / diff(sil.rng) * diff(fm.rng) + fm.rng[1]
    }

    .fm_to_sil <- function(x){
      (x - fm.rng[1]) / diff(fm.rng) * diff(sil.rng) + sil.rng[1]
    }

    plot.df <- rbind(
      data.frame(
        n      = results$n,
        metric = "FM",
        y      = results$FM
      ),
      data.frame(
        n      = results$n,
        metric = "Silhouette",
        y      = .sil_to_fm(results$Silhouette)
      )
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = plot.df,
        ggplot2::aes(x = n, y = y, color = metric),
        linewidth = 0.8,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = plot.df,
        ggplot2::aes(x = n, y = y, color = metric),
        size = 2,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = active.colors) +
      ggplot2::scale_x_log10(breaks = n.vals, labels = n.vals) +
      ggplot2::scale_y_continuous(
        name = "FM",
        sec.axis = ggplot2::sec_axis(
          trans = ~ .fm_to_sil(.),
          name  = "Silhouette"
        )
      ) +
      ggplot2::labs(
        title = plot.title,
        x     = "Number of clusters (n)",
        color = NULL
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title.y.left  = ggplot2::element_text(color = active.colors["FM"]),
        axis.text.y.left   = ggplot2::element_text(color = active.colors["FM"]),
        axis.title.y.right = ggplot2::element_text(color = active.colors["Silhouette"]),
        axis.text.y.right  = ggplot2::element_text(color = active.colors["Silhouette"])
      )

    if (!is.null(best.fm)){
      p <- p +
        ggplot2::geom_point(
          data = data.frame(
            n = best.fm,
            y = results$FM[match(best.fm, results$n)]
          ),
          ggplot2::aes(x = n, y = y),
          color = active.colors["FM"],
          size = 6,
          shape = 21,
          stroke = 2,
          fill = NA,
          show.legend = FALSE
        )
    }#if best.fm

    if (!is.null(best.sil)){
      p <- p +
        ggplot2::geom_point(
          data = data.frame(
            n = best.sil,
            y = .sil_to_fm(results$Silhouette[match(best.sil, results$n)])
          ),
          ggplot2::aes(x = n, y = y),
          color = active.colors["Silhouette"],
          size = 6,
          shape = 21,
          stroke = 2,
          fill = NA,
          show.legend = FALSE
        )
    }#if best.sil

  } else {

    active.cols <- c(
      if (metric == "fm") "FM",
      if (metric == "silhouette") "Silhouette"
    )

    long.df <- data.frame(
      n      = results$n,
      metric = active.cols,
      value  = results[[active.cols]]
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = long.df,
        ggplot2::aes(x = n, y = value, color = metric),
        linewidth = 0.8,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = long.df,
        ggplot2::aes(x = n, y = value, color = metric),
        size = 2,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = active.colors[active.cols]) +
      ggplot2::scale_x_log10(breaks = n.vals, labels = n.vals) +
      ggplot2::labs(
        title = plot.title,
        x     = "Number of clusters (n)",
        y     = active.cols,
        color = NULL
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

    if (metric == "fm" && !is.null(best.fm)){
      p <- p +
        ggplot2::geom_point(
          data = data.frame(
            n = best.fm,
            y = results$FM[match(best.fm, results$n)]
          ),
          ggplot2::aes(x = n, y = y),
          color = active.colors["FM"],
          size = 6,
          shape = 21,
          stroke = 2,
          fill = NA,
          show.legend = FALSE
        )
    }#if best.fm

    if (metric == "silhouette" && !is.null(best.sil)){
      p <- p +
        ggplot2::geom_point(
          data = data.frame(
            n = best.sil,
            y = results$Silhouette[match(best.sil, results$n)]
          ),
          ggplot2::aes(x = n, y = y),
          color = active.colors["Silhouette"],
          size = 6,
          shape = 21,
          stroke = 2,
          fill = NA,
          show.legend = FALSE
        )
    }#if best.sil
  }#if metric == "both"

  print(p)

  invisible(results)
}#n_cluster_metrics


#' Extract LV loadings from integration output
#'
#' Builds a loading matrix from the integration model, optionally writing the
#' top-N feature table to a TSV file and optionally plotting a z-scored bubble heatmap.
#'
#' @export
lv_rank_genes <- function(integration.RData,   # path to integration .RData file
                              output.dir  = NULL,  # directory to write TSV; NULL = skip
                              num.of.top  = NULL,  # top N features per LV; NULL = all features
                              plot        = FALSE, # whether to print the bubble heatmap
                              label.size  = 20){   # base font size for the plot

  load(integration.RData, e <- new.env())
  e <- as.list(e)
  loadings <- e$integration.res$loading

  feature.unique <- unique(unlist(lapply(loadings, names)))
  tab <- matrix(0, nrow = length(feature.unique), ncol = length(loadings))
  rownames(tab) <- feature.unique
  colnames(tab) <- paste0("LV ", seq_along(loadings))

  for (ii in seq_along(loadings))
    tab[names(loadings[[ii]]), ii] <- loadings[[ii]]

  if (!is.null(output.dir)){
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    feature.names <- c(); feature.values <- c()
    feature.index <- c(); LV.index <- c()
    for (ii in seq_along(loadings)){
      n.use <- if (is.null(num.of.top)) length(loadings[[ii]]) else num.of.top
      tmp <- sort(loadings[[ii]], decreasing = TRUE)[seq_len(n.use)]
      feature.names  <- c(feature.names,  names(tmp))
      feature.values <- c(feature.values, as.numeric(tmp))
      feature.index  <- c(feature.index,  seq_len(n.use))
      LV.index       <- c(LV.index,       rep(paste0("LV ", ii), n.use))
    }
    loading.df <- data.frame(LV      = LV.index,
                             Feature = feature.names,
                             Value   = feature.values,
                             Index   = feature.index)
    tsv.name <- if (is.null(num.of.top)) "loading@all.tsv" else paste0("loading@top=", num.of.top, ".tsv")
    write.table(loading.df,
                file.path(output.dir, tsv.name),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Wrote loading TSV to", output.dir, "\n")
  }#if output.dir

  tab.z  <- scale(tab)
  hc     <- hclust(dist(tab.z))
  to_plot <- data.frame(
    val = c(tab.z),
    y   = rep(rownames(tab.z), times = ncol(tab.z)),
    x   = rep(colnames(tab.z), each  = nrow(tab.z))
  )
  to_plot$y <- factor(to_plot$y, levels = rownames(tab.z)[rev(hc$order)])

  p <- ggplot2::ggplot(to_plot, ggplot2::aes(x = x, y = y, color = val, size = val)) +
    ggplot2::geom_point() +
    ggpubr::theme_pubr(base_size = label.size) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::scale_size_area(max_size = 9) +
    ggplot2::labs(color = "loading\n(z scored)", size = "loading\n(z scored)") +
    ggplot2::theme(
      legend.position             = "bottom",
      legend.justification.bottom = "center",
      legend.box.just             = "center",
      legend.location             = "plot",
      legend.box                  = "horizontal",
      axis.text.x  = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid.major = ggplot2::element_line(colour = "gray", linewidth = 0.1),
      axis.line    = ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_gradientn(
      colors = colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(255),
      limits = c(min(to_plot$val), max(to_plot$val))
    )

  if (plot) print(p)
  invisible(list(loading.tab = tab, loading.tab.z = tab.z, plot = p))

}#lv_rank_genes
