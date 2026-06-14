# library(glmnet)
# library(pROC)
# library(PRROC)


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

  load(integration.RData, temp.env1 <- new.env())
  temp.env1 <- as.list(temp.env1)

  embed    <- temp.env1$integration.res$LVs_embeddings
  full.ids <- rownames(embed)
  sec.idx  <- as.numeric(unlist(lapply(strsplit(full.ids, "_"), function(x) x[1])))
  tile.ids <- unlist(lapply(strsplit(full.ids, "_"), function(x) paste0(x[-1], collapse = "_")))
  sec.names <- dataset.list[sec.idx]
  rownames(embed) <- paste0(sec.names, "@", tile.ids)
  mat <- as.matrix(embed)

  res              <- readRDS(clusters.rds)
  clusters         <- res$partition.smooth
  slice.ids        <- res$mat.slice.id
  names(clusters)  <- paste0(slice.ids, "@", names(clusters))

  common <- intersect(rownames(mat), names(clusters))
  if (verbose) cat("Common tiles:", length(common), "\n")

  mat.sub <- mat[common, , drop = F]
  cl.sub  <- clusters[common]
  colnames(mat.sub) <- paste0("LV", 1:ncol(mat.sub))
  mat.sub <- as.matrix(scale(mat.sub))

  sec.sub <- sub("@.*", "", rownames(mat.sub))

  dir.create(output.dir, recursive = T, showWarnings = F)
  clusters.uniq <- sort(unique(cl.sub))
  if (!is.null(clusters.subset))
    clusters.uniq <- intersect(clusters.uniq, clusters.subset)

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

      sec.coefs   <- list()
      sec.metrics <- list()

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

        train.pos <- sample(pos.s, floor(0.8 * length(pos.s)))
        train.neg <- sample(neg.s, floor(0.8 * length(neg.s)))
        train.idx <- c(train.pos, train.neg)
        test.idx  <- setdiff(c(pos.s, neg.s), train.idx)

        X.train <- mat.sub[sec.idx.s[train.idx], ]; y.train <- y.sec[train.idx]
        X.test  <- mat.sub[sec.idx.s[test.idx],  ]; y.test  <- y.sec[test.idx]

        if (is.null(cv.folds)){
          #no CV: fit glmnet for each alpha, pick best alpha+lambda by BIC
          fit.list <- lapply(alpha.sweep, function(a)
            glmnet::glmnet(X.train, y.train, family = "binomial", alpha = a))
          bic.min  <- sapply(fit.list, function(fit){
            min((1 - fit$dev.ratio) * fit$nulldev + log(nrow(X.train)) * fit$df)
          })
          best.i   <- which.min(bic.min)
          fit.best <- fit.list[[best.i]]
          bic      <- (1 - fit.best$dev.ratio) * fit.best$nulldev + log(nrow(X.train)) * fit.best$df
          best.l   <- fit.best$lambda[which.min(bic)]
          coefs    <- coef(fit.best, s = best.l)[-1]
          best.a   <- alpha.sweep[best.i]
          pred.prob <- as.numeric(predict(fit.best, X.test, s = best.l, type = "response"))
        }else{
          #CV: sweep alpha values, pick best alpha by min CV error
          cv.list  <- lapply(alpha.sweep, function(a)
            glmnet::cv.glmnet(X.train, y.train, family = "binomial", alpha = a, nfolds = cv.folds))
          best.i   <- which.min(sapply(cv.list, function(f) min(f$cvm)))
          cv.fit   <- cv.list[[best.i]]
          coefs    <- coef(cv.fit, s = "lambda.min")[-1]
          best.a   <- alpha.sweep[best.i]
          pred.prob <- as.numeric(predict(cv.fit, X.test, s = "lambda.min", type = "response"))
        }#else

        names(coefs)      <- lv.names
        sec.coefs[[sec]]  <- abs(coefs)

        pred.class <- ifelse(pred.prob > 0.5, 1, 0)
        sec.metrics[[sec]] <- list(
          accuracy         = mean(pred.class == y.test),
          AUC              = tryCatch(as.numeric(pROC::auc(pROC::roc(y.test, pred.prob, quiet = T))),
                                      error = function(e) NA_real_),
          AUPRC            = tryCatch(PRROC::pr.curve(scores.class0 = pred.prob[y.test == 1],
                                                      scores.class1 = pred.prob[y.test == 0],
                                                      curve = F)$auc.integral,
                                      error = function(e) NA_real_),
          best_alpha       = best.a,
          downsample_ratio = if (!is.null(downsample)) eq.ratio else NA_real_
        )
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

      mean.acc   <- round(mean(sapply(sec.metrics, `[[`, "accuracy"),  na.rm = T), 4)
      mean.auc   <- round(mean(sapply(sec.metrics, `[[`, "AUC"),       na.rm = T), 4)
      mean.auprc <- round(mean(sapply(sec.metrics, `[[`, "AUPRC"),     na.rm = T), 4)
      mean.ds    <- mean(sapply(sec.metrics, `[[`, "downsample_ratio"), na.rm = T)

      en.rows[[cl]] <- data.frame(
        cluster            = cl,
        rank               = 1:use.n,
        LV                 = top.lvs,
        mean_coef_abs      = mean.coef.abs[top.idx],
        n_sections_nonzero = n.nonzero[top.idx],
        n_sections_total   = n.sec,
        accuracy           = mean.acc,
        AUC                = mean.auc,
        AUPRC              = mean.auprc,
        downsample_ratio   = if (!is.null(downsample)) mean.ds else NA_real_,
        stringsAsFactors   = F
      )

      if (verbose) message(sprintf("done (%d sections)", n.sec))
    }#for cl

    en.df <- do.call(rbind, en.rows)
    en.df$cluster <- tryCatch(as.integer(en.df$cluster),
      warning = function(w) tryCatch(as.numeric(en.df$cluster),
        warning = function(w) en.df$cluster))
    en.df <- en.df[order(en.df$cluster, en.df$rank), ]
    write.table(en.df,
                file.path(output.dir, "elasticnet_topLVs.tsv"),
                sep = "\t", quote = F, row.names = F)
    cat("Wrote elasticnet_topLVs.tsv\n")
  }#if elastic.net

  invisible(NULL)

}#find_top_lvs


#' Evaluate cluster quality across n-cluster solutions
#'
#' Loads pre-computed Louvain clusterings and computes FM index and Silhouette
#' score across cluster counts, then plots the curves and reports the optimal n.
#' FM is computed on the full cluster labels; Silhouette uses a subsampled LV
#' embedding (controlled by subsample).
#'
#' @export
eval_clusters <- function(exhaustive.dir,                       # path to exhaustive output (to resolve section names)
                          integration.RData,                    # path to integration .RData file
                          cluster.dir,                          # path to cluster/ folder containing .RDS files
                          model.name  = "all@FineTune@first@default", # .RData filename (no extension)
                          k.opt       = 40,                     # which SNN k to select from cluster files
                          subsample   = 10000,                   # max tiles for silhouette; NULL = all
                          seed        = 42,                     # random seed for subsampling
                          title       = NULL,                   # plot title; NULL = auto
                          verbose     = T){                     # print per-n progress

  set.seed(seed)

  #load LV matrix from integration output
  dataset.list <- list.files(exhaustive.dir)
  load(paste0(integration.RData), e <- new.env())
  e     <- as.list(e)
  embed <- e$integration.res$LVs_embeddings

  sids         <- unlist(lapply(strsplit(rownames(embed), "_"), function(x) paste0(x[-1], collapse = "_")))
  sec.idx      <- unlist(lapply(strsplit(rownames(embed), "_"), function(x) x[1]))
  mat.slice.id <- dataset.list[as.numeric(sec.idx)]
  rownames(embed) <- sids
  mat <- L2Norm(as.matrix(embed), MARGIN = 1)

  #optionally subsample for silhouette (O(n^2) memory)
  if (!is.null(subsample) && nrow(mat) > subsample){
    mat.sub <- mat[sample(nrow(mat), subsample), ]
  }else{
    mat.sub <- mat
  }#else

  #FM index helper: pair-counting similarity between two label vectors
  .fm_index <- function(a, b){
    tab <- table(a, b)
    TP  <- sum(choose(tab, 2))
    FP  <- sum(choose(rowSums(tab), 2)) - TP
    FN  <- sum(choose(colSums(tab), 2)) - TP
    TP / sqrt((TP + FP) * (TP + FN))
  }#.fm_index

  #load pre-computed clusterings filtered by k.opt
  files  <- list.files(cluster.dir, pattern = paste0("k=", k.opt), full.names = T)
  if (length(files) == 0) stop("No files found in cluster.dir matching k=", k.opt)
  n.vals <- as.numeric(sub(".*n=([0-9]+)\\.RDS", "\\1", basename(files)))
  files  <- files[order(n.vals)]
  n.vals <- sort(n.vals)

  if (verbose) cat("Loading", length(files), "clusterings (k=", k.opt, ")...\n")
  clusterings <- lapply(files, function(f) readRDS(f)$partition.smooth)
  names(clusterings) <- as.character(n.vals)

  #compute FM + silhouette per n
  results <- do.call(rbind, lapply(seq_along(n.vals), function(i){
    cl.full <- clusterings[[i]]

    #FM: compare to adjacent n solutions
    fmi.vals <- numeric(0)
    if (i > 1){
      common2  <- intersect(names(cl.full), names(clusterings[[i - 1]]))
      if (length(common2) > 0)
        fmi.vals <- c(fmi.vals, .fm_index(cl.full[common2], clusterings[[i - 1]][common2]))
    }#if
    if (i < length(n.vals)){
      common2  <- intersect(names(cl.full), names(clusterings[[i + 1]]))
      if (length(common2) > 0)
        fmi.vals <- c(fmi.vals, .fm_index(cl.full[common2], clusterings[[i + 1]][common2]))
    }#if
    fm <- if (length(fmi.vals) > 0) mean(fmi.vals) else NA_real_

    #silhouette: restrict to mat.sub pixels that were clustered
    common <- intersect(rownames(mat.sub), names(cl.full))
    sil    <- NA_real_
    if (length(common) >= 2){
      cl     <- cl.full[common]
      sm     <- mat.sub[common, ]
      cl.int <- as.integer(as.factor(cl))
      if (length(unique(cl.int)) >= 2)
        sil <- mean(cluster::silhouette(cl.int, dist(sm))[, 3])
    }#if

    if (verbose) cat(sprintf("n=%d  FM=%.3f  sil=%.3f\n", n.vals[i], fm, sil))
    data.frame(n = n.vals[i], FM = fm, Silhouette = sil)
  }))

  #optimal n
  best.sil <- if (any(!is.na(results$Silhouette))) results$n[which.max(results$Silhouette)] else NULL
  best.fm  <- if (any(!is.na(results$FM)))         results$n[which.max(ifelse(is.na(results$FM), -Inf, results$FM))] else NULL
  if (verbose) cat(sprintf("\nOptimal n — Silhouette: %s | FM: %s\n",
                            if (is.null(best.sil)) "NA" else best.sil,
                            if (is.null(best.fm))  "NA" else best.fm))

  #plot
  long.df <- reshape(results[, c("n", "FM", "Silhouette")],
                     varying = c("FM", "Silhouette"), v.names = "value",
                     timevar = "metric", times = c("FM", "Silhouette"), direction = "long")

  highlight.rows <- list()
  if (!is.null(best.sil))
    highlight.rows[["Silhouette"]] <- data.frame(
      n = best.sil, metric = "Silhouette",
      y = results$Silhouette[match(best.sil, results$n)])
  if (!is.null(best.fm))
    highlight.rows[["FM"]] <- data.frame(
      n = best.fm, metric = "FM",
      y = results$FM[match(best.fm, results$n)])

  plot.title <- if (!is.null(title)) title else
    paste0("Cluster quality — k=", k.opt,
           if (!is.null(subsample)) paste0(" (sil subsample=", subsample, ")") else "")

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data  = long.df,
                       ggplot2::aes(x = n, y = value, color = metric),
                       linewidth = 0.8, na.rm = T) +
    ggplot2::geom_point(data = long.df,
                        ggplot2::aes(x = n, y = value, color = metric),
                        size = 2, na.rm = T) +
    ggplot2::scale_color_manual(values = c("FM" = "purple", "Silhouette" = "forestgreen")) +
    ggplot2::scale_x_log10(breaks = n.vals, labels = n.vals) +
    ggplot2::labs(title = plot.title, x = "Number of clusters (n)", y = "Score", color = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (length(highlight.rows) > 0){
    hl <- do.call(rbind, highlight.rows)
    p  <- p + ggplot2::geom_point(data = hl,
                                   ggplot2::aes(x = n, y = y, color = metric),
                                   size = 6, shape = 21, stroke = 2,
                                   fill = NA, show.legend = F)
  }#if

  print(p)
  invisible(results)

}#eval_clusters
