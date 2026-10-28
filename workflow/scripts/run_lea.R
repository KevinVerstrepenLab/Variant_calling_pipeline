#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(LEA)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file) && nzchar(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "clade_palette_utils.R"))

option_list <- list(
  make_option(c("--ped"), type = "character"),
  make_option(c("--fam"), type = "character"),
  make_option(c("--out"), type = "character", default = "lea"),
  make_option(c("--k-min"), type = "integer", default = 2, dest = "k_min"),
  make_option(c("--k-max"), type = "integer", default = 20, dest = "k_max"),
  make_option(c("--repetitions"), type = "integer", default = 5),
  make_option(c("--seed"), type = "integer", default = 12345),
  make_option(c("--threads"), type = "integer", default = 1),
  make_option(c("--best-k"), type = "character", default = "auto", dest = "best_k"),
  make_option(c("--plot-k-min"), type = "integer", default = NA, dest = "plot_k_min"),
  make_option(c("--plot-k-max"), type = "integer", default = NA, dest = "plot_k_max"),
  make_option(c("--order-by"), type = "character", default = "bestK", dest = "order_by"),
  make_option(c("--width"), type = "double", default = 20),
  make_option(c("--height"), type = "double", default = 6),
  make_option(c("--multi-height"), type = "double", default = 18, dest = "multi_height"),
  make_option(c("--rename"), type = "character", default = NULL),
  make_option(c("--clades"), type = "character", default = NULL),
  make_option(c("--clade_colors"), type = "character", default = NULL),
  make_option(c("--color_sample_labels_by_clade"), action = "store_true", default = FALSE),
  make_option(c("--clade_color_strip"), action = "store_true", default = FALSE),
  make_option(c("--model-only"), action = "store_true", default = FALSE, dest = "model_only"),
  make_option(c("--summary-only"), action = "store_true", default = FALSE, dest = "summary_only")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$ped) || is.null(opt$fam)) {
  stop("--ped and --fam are required", call. = FALSE)
}
if (!file.exists(opt$ped)) stop("PED file not found: ", opt$ped, call. = FALSE)
if (!file.exists(opt$fam)) stop("FAM file not found: ", opt$fam, call. = FALSE)
if (isTRUE(opt$model_only) && isTRUE(opt$summary_only)) {
  stop("Use either --model-only or --summary-only, not both", call. = FALSE)
}

read_table_flexible <- function(path, header = FALSE) {
  x <- try(read.table(path, header = header, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, comment.char = ""), silent = TRUE)
  if (inherits(x, "try-error")) {
    x <- read.table(path, header = header, sep = "", stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "")
  }
  x
}

read_two_col <- function(path, colnames_out) {
  x <- read_table_flexible(path, header = FALSE)
  if (ncol(x) < 2) stop("Expected at least two columns in ", path, call. = FALSE)
  x <- x[, 1:2, drop = FALSE]
  colnames(x) <- colnames_out
  x
}

apply_rename <- function(x, rename_df) {
  if (is.null(rename_df)) return(x)
  idx <- match(x, rename_df$old)
  replace <- !is.na(idx)
  x[replace] <- rename_df$new[idx[replace]]
  x
}

read_clades_file <- function(path, rename_df) {
  if (is.null(path) || !file.exists(path) || file.info(path)$size == 0) return(NULL)
  x <- read_table_flexible(path, header = FALSE)
  if (ncol(x) < 2) stop("Expected at least two columns in clades file ", path, call. = FALSE)
  x <- x[, 1:2, drop = FALSE]
  first <- tolower(trimws(as.character(unlist(x[1, ], use.names = FALSE))))
  if (length(first) >= 2 && first[1] %in% c("accession", "sample", "sample_id", "label", "tip", "id") &&
      first[2] %in% c("clade", "group", "population", "pop", "lineage")) {
    x <- x[-1, , drop = FALSE]
  }
  colnames(x) <- c("label", "clade")
  x$label <- apply_rename(trimws(as.character(x$label)), rename_df)
  x$clade <- trimws(as.character(x$clade))
  x$clade[tolower(x$clade) == "unknown"] <- NA
  x <- x[!is.na(x$label) & nzchar(x$label) & !is.na(x$clade) & nzchar(x$clade), , drop = FALSE]
  if (nrow(x) == 0) return(NULL)
  x[!duplicated(x$label), , drop = FALSE]
}

set.seed(opt$seed)
ks <- seq(opt$k_min, opt$k_max)
if (length(ks) < 1) stop("Invalid K range", call. = FALSE)
plot_k_min <- ifelse(is.na(opt$plot_k_min), min(ks), opt$plot_k_min)
plot_k_max <- ifelse(is.na(opt$plot_k_max), max(ks), opt$plot_k_max)
plot_ks <- seq(plot_k_min, plot_k_max)
plot_ks <- plot_ks[plot_ks %in% ks]

fam <- read.table(opt$fam, header = FALSE, stringsAsFactors = FALSE)
if (ncol(fam) < 2) stop("FAM file must have at least two columns", call. = FALSE)
sample_names <- trimws(as.character(fam[[2]]))

rename_df <- NULL
if (!is.null(opt$rename)) {
  rename_df <- read_two_col(opt$rename, c("old", "new"))
  rename_df$old <- trimws(as.character(rename_df$old))
  rename_df$new <- trimws(as.character(rename_df$new))
}
sample_labels <- apply_rename(sample_names, rename_df)

clade_df <- read_clades_file(opt$clades, rename_df)
clade_palette_override <- read_clade_colors_file(opt$clade_colors)
sample_clades <- rep(NA_character_, length(sample_names))
clade_colors <- NULL
if (!is.null(clade_df)) {
  sample_clades <- clade_df$clade[match(sample_labels, clade_df$label)]
  sample_clades[is.na(sample_clades) | !nzchar(sample_clades)] <- "Unclassified"
  clade_colors <- make_clade_palette(c(clade_df$clade, sample_clades), palette = clade_palette_override)
  sample_clade_colors <- unname(clade_colors[sample_clades])
  sample_clade_colors[is.na(sample_clade_colors)] <- "grey60"
  write_clade_palette(clade_colors, paste0(opt$out, ".clade_colors.tsv"))
  write.table(
    data.frame(sample = sample_names, label = sample_labels, clade = sample_clades, color = sample_clade_colors, stringsAsFactors = FALSE, check.names = FALSE),
    paste0(opt$out, ".sample_clades.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  message("LEA clade matches: ", sum(sample_clades != "Unclassified"), " / ", length(sample_clades))
}

geno <- paste0(opt$out, ".geno")
project_file <- paste0(opt$out, ".snmfProject")
ce_file <- paste0(opt$out, ".cross_entropy.tsv")
ce_summary_file <- paste0(opt$out, ".cross_entropy_summary.tsv")
summary_done_file <- paste0(opt$out, ".summary.done")

load_or_stop_project <- function(path) {
  if (!file.exists(path)) {
    stop("sNMF project file not found: ", path,
         ". Re-run without --summary-only or delete the checkpoint stamp.", call. = FALSE)
  }
  load.snmfProject(path)
}

cross_entropy_values <- function(x) {
  # LEA versions differ slightly in how cross.entropy() returns values.
  # Some return a plain numeric vector; others may return a named vector,
  # matrix, or data frame. Convert everything to a simple numeric vector so
  # rbind never fails because of inconsistent column names.
  if (is.null(x)) return(numeric(0))
  if (is.data.frame(x) || is.matrix(x) || is.array(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  x <- suppressWarnings(as.numeric(x))
  # Keep original run positions: dropping a failed run would renumber Q-matrix selection.
  x[!is.finite(x)] <- NA_real_
  x
}

write_cross_entropy_tables <- function(project, ks, ce_file, ce_summary_file) {
  ce_rows <- list()
  for (k in ks) {
    ce_values <- cross_entropy_values(cross.entropy(project, K = k))
    if (!any(is.finite(ce_values))) {
      stop("No usable cross-entropy values for K=", k, call. = FALSE)
    }
    ce_rows[[as.character(k)]] <- data.frame(
      K = rep(as.integer(k), length(ce_values)),
      run = seq_along(ce_values),
      cross_entropy = ce_values,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (length(ce_rows) == 0) {
    stop("No usable cross-entropy values were returned by LEA", call. = FALSE)
  }
  ce <- do.call(rbind, unname(ce_rows))
  rownames(ce) <- NULL
  write.table(ce, ce_file, sep = "\t", quote = FALSE, row.names = FALSE)

  ce_summary <- aggregate(cross_entropy ~ K, data = ce, FUN = function(x) mean(x, na.rm = TRUE))
  colnames(ce_summary)[2] <- "mean_cross_entropy"
  ce_summary <- ce_summary[order(ce_summary$K), , drop = FALSE]
  write.table(ce_summary, ce_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  list(ce = ce, ce_summary = ce_summary)
}

if (isTRUE(opt$summary_only)) {
  cat("Loading existing sNMF project:\n")
  cat("  ", project_file, "\n", sep = "")
  project <- load_or_stop_project(project_file)
  if (file.exists(ce_file) && file.exists(ce_summary_file)) {
    ce <- read.table(ce_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    ce_summary <- read.table(ce_summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    tmp <- write_cross_entropy_tables(project, ks, ce_file, ce_summary_file)
    ce <- tmp$ce
    ce_summary <- tmp$ce_summary
  }
} else {
  cat("Converting PED to LEA geno:\n")
  cat("  ", opt$ped, " -> ", geno, "\n", sep = "")
  ped2geno(opt$ped, geno)

  cat("Running snmf for K = ", paste(ks, collapse = ","), " with ", opt$repetitions, " repetitions\n", sep = "")
  project <- snmf(geno, K = ks, repetitions = opt$repetitions, entropy = TRUE, project = "new", CPU = opt$threads)

  tmp <- write_cross_entropy_tables(project, ks, ce_file, ce_summary_file)
  ce <- tmp$ce
  ce_summary <- tmp$ce_summary

  # LEA normally writes this file, but the exact naming can vary between versions.
  # If needed, try to find the project file that corresponds to this output prefix.
  if (!file.exists(project_file)) {
    candidates <- unique(c(
      paste0(sub("[.]geno$", "", geno), ".snmfProject"),
      paste0(geno, ".snmfProject"),
      list.files(dirname(opt$out), pattern = paste0(basename(opt$out), ".*snmfProject$"), full.names = TRUE)
    ))
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates) > 0) {
      file.copy(candidates[1], project_file, overwrite = TRUE)
    }
  }
}

if (isTRUE(opt$model_only)) {
  cat("Wrote model files:\n")
  cat("  ", geno, "\n", sep = "")
  cat("  ", project_file, "\n", sep = "")
  cat("  ", ce_file, "\n", sep = "")
  cat("  ", ce_summary_file, "\n", sep = "")
  quit(save = "no", status = 0)
}

if (nrow(ce_summary) == 0) stop("No cross-entropy summary values available", call. = FALSE)
auto_best <- ce_summary$K[which.min(ce_summary$mean_cross_entropy)]
best_k <- if (tolower(opt$best_k) == "auto") auto_best else as.integer(opt$best_k)
if (!(best_k %in% ks)) stop("best_k is not within tested K range", call. = FALSE)
best_run <- which.min(cross_entropy_values(cross.entropy(project, K = best_k)))
writeLines(c(paste0("bestK=", best_k), paste0("bestRun=", best_run)), paste0(opt$out, ".bestK.txt"))

pdf(paste0(opt$out, ".cross_entropy.pdf"), width = 7, height = 5)
plot(ce_summary$K, ce_summary$mean_cross_entropy, type = "b", pch = 19,
     xlab = "Number of ancestral populations (K)", ylab = "Mean cross-entropy criterion")
abline(v = best_k, lty = 2)
dev.off()

q_best <- Q(project, K = best_k, run = best_run)
if (nrow(q_best) != length(sample_names)) {
  warning("Q-matrix row count does not match sample count in FAM file")
}
colnames(q_best) <- paste0("Cluster", seq_len(ncol(q_best)))
q_best_df <- data.frame(sample = sample_names, label = sample_labels, clade = sample_clades, q_best, check.names = FALSE)
write.table(q_best_df, paste0(opt$out, ".K", best_k, ".Q.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

dominant_cluster <- apply(q_best, 1, which.max)
order_by <- tolower(opt$order_by)
if (order_by %in% c("bestk", "dominant", "dominant_cluster")) {
  order_idx <- order(dominant_cluster)
} else if (order_by %in% c("clade", "clades")) {
  order_idx <- order(sample_clades, dominant_cluster, sample_labels)
} else if (order_by %in% c("clade_bestk", "clade_then_bestk", "clade_dominant", "clade_then_dominant")) {
  order_idx <- order(sample_clades, dominant_cluster, sample_labels)
} else {
  order_idx <- seq_along(sample_names)
}

cluster_cols <- function(k) rainbow(k)

add_clade_annotations <- function(bp, ordered_samples, ordered_clades, show_labels = TRUE) {
  if (is.null(clade_colors) || all(is.na(ordered_clades))) {
    if (show_labels) {
      usr <- par("usr")
      text(bp, usr[3] - 0.015 * diff(usr[3:4]), labels = ordered_samples, srt = 90,
           adj = 1, xpd = NA, cex = 0.6, col = "black")
    }
    return(invisible(NULL))
  }

  ordered_cols <- unname(clade_colors[ordered_clades])
  ordered_cols[is.na(ordered_cols)] <- "grey60"
  usr <- par("usr")
  y_range <- diff(usr[3:4])

  if (isTRUE(opt$clade_color_strip)) {
    half_width <- if (length(bp) > 1) min(diff(bp)) / 2 else 0.5
    rect(bp - half_width, -0.055, bp + half_width, -0.010, col = ordered_cols, border = NA, xpd = NA)
  }

  if (show_labels) {
    label_cols <- if (isTRUE(opt$color_sample_labels_by_clade)) ordered_cols else rep("black", length(ordered_samples))
    text(bp, usr[3] - 0.015 * y_range, labels = ordered_samples, srt = 90,
         adj = 1, xpd = NA, cex = 0.6, col = label_cols)
  }
  invisible(NULL)
}

plot_single_admixture <- function(qmatrix, k, order_idx, samples, file) {
  pdf(file, width = opt$width, height = opt$height)
  par(mar = c(9, 4, 2, 1))
  ylim <- if (!is.null(clade_colors) && isTRUE(opt$clade_color_strip)) c(-0.07, 1) else c(0, 1)
  bp <- barplot(t(qmatrix[order_idx, , drop = FALSE]), col = cluster_cols(k), border = NA,
                names.arg = rep("", length(samples)), las = 2, cex.names = 0.6, space = 0,
                ylim = ylim, ylab = "Ancestry proportion", main = paste("LEA sNMF K =", k))
  add_clade_annotations(bp, samples[order_idx], sample_clades[order_idx], show_labels = TRUE)
  dev.off()
}

plot_single_admixture(q_best, best_k, order_idx, sample_labels, paste0(opt$out, ".K", best_k, ".admixture.pdf"))

if (length(plot_ks) > 0) {
  pdf(paste0(opt$out, ".K", min(plot_ks), "_to_K", max(plot_ks), ".admixture.pdf"),
      width = opt$width, height = opt$multi_height)
  par(mfrow = c(length(plot_ks), 1), mar = c(0.3, 4, 0.5, 6), oma = c(8, 0, 1, 0))
  ylim <- if (!is.null(clade_colors) && isTRUE(opt$clade_color_strip)) c(-0.07, 1) else c(0, 1)
  for (k in plot_ks) {
    run_k <- which.min(cross_entropy_values(cross.entropy(project, K = k)))
    q_k <- Q(project, K = k, run = run_k)
    show_labels <- k == max(plot_ks)
    bp <- barplot(t(q_k[order_idx, , drop = FALSE]), col = cluster_cols(k), border = NA,
                  names.arg = rep("", length(sample_names)), las = 2, cex.names = 0.7, space = 0,
                  ylim = ylim, ylab = "Q")
    add_clade_annotations(bp, sample_labels[order_idx], sample_clades[order_idx], show_labels = show_labels)
    mtext(paste("K =", k), side = 4, line = 1, cex = 1.0, las = 1)
  }
  mtext("Samples", side = 1, line = 6, outer = TRUE)
  dev.off()
}

cat("Wrote:\n")
cat("  ", geno, "\n", sep = "")
cat("  ", ce_file, "\n", sep = "")
cat("  ", ce_summary_file, "\n", sep = "")
cat("  ", paste0(opt$out, ".bestK.txt"), "\n", sep = "")
cat("  ", paste0(opt$out, ".K", best_k, ".Q.tsv"), "\n", sep = "")
cat("  ", paste0(opt$out, ".K", best_k, ".admixture.pdf"), "\n", sep = "")
if (!is.null(clade_colors)) {
  cat("  ", paste0(opt$out, ".clade_colors.tsv"), "\n", sep = "")
  cat("  ", paste0(opt$out, ".sample_clades.tsv"), "\n", sep = "")
}
stopifnot(file.copy(paste0(opt$out, ".K", best_k, ".Q.tsv"),
                   paste0(opt$out, ".best.Q.tsv"), overwrite = TRUE))
stopifnot(file.copy(paste0(opt$out, ".K", best_k, ".admixture.pdf"),
                   paste0(opt$out, ".best.admixture.pdf"), overwrite = TRUE))
writeLines(c(
  paste0("completed=", format(Sys.time(), usetz = TRUE)),
  paste0("bestK=", best_k),
  paste0("bestRun=", best_run),
  paste0("clade_colors=", ifelse(!is.null(clade_colors), paste0(opt$out, ".clade_colors.tsv"), ""))
), summary_done_file)
cat("  ", summary_done_file, "\n", sep = "")
