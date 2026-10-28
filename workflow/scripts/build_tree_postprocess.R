#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(ggplot2)
  library(dplyr)
})

has_ggtree <- requireNamespace("ggtree", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggtree) suppressPackageStartupMessages(library(ggtree))

option_list <- list(
  make_option(c("--dist"), type = "character"),
  make_option(c("--ids"), type = "character"),
  make_option(c("--method"), type = "character", default = "nj"),
  make_option(c("--rename"), type = "character", default = NULL),
  make_option(c("--pca_sample_list"), type = "character", default = NULL),
  make_option(c("--tree_sample_list"), type = "character", default = NULL),
  make_option(c("--outgroup"), type = "character", default = NULL),
  make_option(c("--clades"), type = "character", default = NULL),
  make_option(c("--clade_colors"), type = "character", default = NULL),
  make_option(c("--out"), type = "character", default = "tree"),
  make_option(c("--layout"), type = "character", default = "fan"),
  make_option(c("--branch_length"), type = "character", default = "yes"),
  make_option(c("--width"), type = "double", default = 18),
  make_option(c("--height"), type = "double", default = 18),
  make_option(c("--dpi"), type = "integer", default = 300),
  make_option(c("--label_tips"), action = "store_true", default = FALSE),
  make_option(c("--label_size"), type = "double", default = 2.0),
  make_option(c("--tip_size"), type = "double", default = 1.2),
  make_option(c("--color_branches"), action = "store_true", default = FALSE),
  make_option(c("--color_tip_labels"), action = "store_true", default = FALSE),
  make_option(c("--pca_eigenvec"), type = "character", default = NULL),
  make_option(c("--pca_eigenval"), type = "character", default = NULL),
  make_option(c("--pca_pc_x"), type = "integer", default = 1),
  make_option(c("--pca_pc_y"), type = "integer", default = 2),
  make_option(c("--pca_point_size"), type = "double", default = 1.8),
  make_option(c("--pca_label_points"), type = "character", default = NULL),
  make_option(c("--pca_label_size"), type = "double", default = NULL),
  make_option(c("--pca_alpha"), type = "double", default = 0.9)
)

opt <- parse_args(OptionParser(option_list = option_list))
has_tree <- !is.null(opt$dist) && !is.null(opt$ids)
if (!has_tree && is.null(opt$pca_eigenvec)) {
  stop("Provide --dist and --ids, or --pca_eigenvec", call. = FALSE)
}

as_bool_option <- function(x, name) {
  if (is.null(x)) return(TRUE)
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("1", "true", "yes", "y", "on")) return(TRUE)
  if (s %in% c("0", "false", "no", "n", "off", "none")) return(FALSE)
  stop(name, " must be yes or no", call. = FALSE)
}

branch_length_enabled <- as_bool_option(opt$branch_length, "--branch_length")
branch_length_mode <- if (branch_length_enabled) "branch.length" else "none"
tree_plot_prefix <- paste0(opt$out, "-", opt$layout)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file) && nzchar(script_file)) dirname(normalizePath(script_file)) else "."
source(file.path(script_dir, "clade_palette_utils.R"))

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

read_clades_file <- function(path, rename_df) {
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

make_tree_node_annotations <- function(tr, tip_df) {
  n_tip <- length(tr$tip.label)
  n_total <- n_tip + tr$Nnode
  node_ids <- seq_len(n_total)

  if (!"clade_label" %in% names(tip_df)) {
    tip_df$clade_label <- NA_character_
  }
  if (!"display_label" %in% names(tip_df)) {
    tip_df$display_label <- tip_df$label
  }

  tip_match <- match(tr$tip.label, tip_df$label)
  tip_clades <- tip_df$clade_label[tip_match]
  tip_clades[is.na(tip_clades) | !nzchar(tip_clades)] <- "Unclassified"
  display_labels <- tip_df$display_label[tip_match]
  display_labels[is.na(display_labels) | !nzchar(display_labels)] <- tr$tip.label[is.na(display_labels) | !nzchar(display_labels)]

  child_map <- split(tr$edge[, 2], tr$edge[, 1])
  desc_clades <- function(node) {
    if (node <= n_tip) return(tip_clades[node])
    kids <- child_map[[as.character(node)]]
    if (is.null(kids) || length(kids) == 0) return(NA_character_)
    unlist(lapply(kids, desc_clades), use.names = FALSE)
  }

  branch_clade <- rep(NA_character_, n_total)
  tip_clade <- rep(NA_character_, n_total)
  display_label <- rep(NA_character_, n_total)
  tip_clade[seq_len(n_tip)] <- tip_clades
  display_label[seq_len(n_tip)] <- display_labels

  for (node in node_ids) {
    cls <- unique(desc_clades(node))
    cls <- cls[!is.na(cls) & nzchar(cls)]
    if (length(cls) == 1) {
      branch_clade[node] <- cls
    }
  }

  data.frame(
    node = node_ids,
    branch_clade = branch_clade,
    tip_clade = tip_clade,
    display_label = display_label,
    stringsAsFactors = FALSE
  )
}

apply_rename <- function(x, rename_df) {
  if (is.null(rename_df)) return(x)
  idx <- match(x, rename_df$old)
  replace <- !is.na(idx)
  x[replace] <- rename_df$new[idx[replace]]
  x
}

# make_clade_palette(), read_clade_colors_file(), and write_clade_palette() are
# sourced from scripts/clade_palette_utils.R so tree, PCA, and LEA share one palette.

is_numeric_like <- function(x) {
  suppressWarnings(all(!is.na(as.numeric(as.character(x)))))
}

read_plink_eigenvec <- function(path) {
  raw <- read_table_flexible(path, header = FALSE)
  if (nrow(raw) == 0 || ncol(raw) < 3) {
    stop("Expected PLINK eigenvec file with sample IDs and PCs in ", path, call. = FALSE)
  }

  first <- trimws(as.character(unlist(raw[1, ], use.names = FALSE)))
  has_header <- first[1] %in% c("#FID", "FID", "#IID", "IID") || any(grepl("^PC[0-9]+$", first))
  if (has_header) {
    header <- sub("^#", "", first)
    raw <- raw[-1, , drop = FALSE]
    colnames(raw) <- header
  }

  if (has_header && all(c("FID", "IID") %in% colnames(raw))) {
    eig <- raw
    pc_cols <- grep("^PC[0-9]+$", colnames(eig), value = TRUE)
  } else if (has_header && "IID" %in% colnames(raw) && !("FID" %in% colnames(raw))) {
    eig <- data.frame(FID = raw[["IID"]], IID = raw[["IID"]], raw[, setdiff(colnames(raw), "IID"), drop = FALSE], check.names = FALSE)
    pc_cols <- grep("^PC[0-9]+$", colnames(eig), value = TRUE)
  } else {
    # Headerless PLINK files can be either FID IID PC1... or IID PC1...
    # In the latter case column 2 is numeric, which was the source of the
    # unclassified PCA: it was being treated as IID instead of PC1.
    one_id_format <- is_numeric_like(raw[[2]])
    if (one_id_format) {
      pc_n <- ncol(raw) - 1
      colnames(raw) <- c("FID", paste0("PC", seq_len(pc_n)))
      eig <- data.frame(FID = raw$FID, IID = raw$FID, raw[, paste0("PC", seq_len(pc_n)), drop = FALSE], check.names = FALSE)
    } else {
      pc_n <- ncol(raw) - 2
      colnames(raw) <- c("FID", "IID", paste0("PC", seq_len(pc_n)))
      eig <- raw
    }
    pc_cols <- grep("^PC[0-9]+$", colnames(eig), value = TRUE)
  }

  if (length(pc_cols) < 2) {
    stop("Expected at least PC1 and PC2 in ", path, call. = FALSE)
  }
  for (pc in pc_cols) eig[[pc]] <- as.numeric(as.character(eig[[pc]]))
  eig$FID <- trimws(as.character(eig$FID))
  eig$IID <- trimws(as.character(eig$IID))
  eig
}

choose_pca_clade_key <- function(eig, clade_df, rename_df) {
  fid_name <- apply_rename(trimws(eig$FID), rename_df)
  iid_name <- apply_rename(trimws(eig$IID), rename_df)

  count_matches <- function(x) sum(x %in% clade_df$label)
  fid_matches <- count_matches(fid_name)
  iid_matches <- count_matches(iid_name)

  # Prefer the original working behaviour, FID. Fall back to IID only when it
  # gives more clade matches. This handles both one-ID and two-ID PLINK output.
  if (iid_matches > fid_matches) {
    list(name = iid_name, key = "IID", matches = iid_matches, fid_matches = fid_matches, iid_matches = iid_matches)
  } else {
    list(name = fid_name, key = "FID", matches = fid_matches, fid_matches = fid_matches, iid_matches = iid_matches)
  }
}

rename_df <- NULL
if (!is.null(opt$rename)) {
  rename_df <- read_two_col(opt$rename, c("old", "new"))
  rename_df$old <- trimws(rename_df$old)
  rename_df$new <- trimws(rename_df$new)
}

read_sample_list <- function(path) {
  if (is.null(path)) return(NULL)
  samples <- trimws(readLines(path, warn = FALSE))
  samples <- unique(samples[nzchar(samples) & !startsWith(samples, "#")])
  if (length(samples) == 0) stop("No sample IDs in sample_list: ", path, call. = FALSE)
  unique(apply_rename(samples, rename_df))
}

sample_list_matches <- function(requested, available, analysis) {
  missing <- setdiff(requested, available)
  if (length(missing) > 0) {
    warning(analysis, " sample_list IDs not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  matched <- intersect(requested, available)
  if (length(matched) == 0) stop(analysis, " sample_list matched no samples", call. = FALSE)
  matched
}

clade_df <- NULL
if (!is.null(opt$clades) && file.exists(opt$clades) && file.info(opt$clades)$size > 0) {
  clade_df <- read_clades_file(opt$clades, rename_df)
}
clade_colors <- read_clade_colors_file(opt$clade_colors)
if (has_tree) {
  ids <- read_table_flexible(opt$ids, header = FALSE)
  if (ncol(ids) < 2) stop("PLINK .mdist.id file must have at least two columns", call. = FALSE)
  labels <- trimws(ids[[2]])
  labels <- apply_rename(labels, rename_df)

  M <- as.matrix(read_table_flexible(opt$dist, header = FALSE))
  storage.mode(M) <- "numeric"
  if (nrow(M) != length(labels) || ncol(M) != length(labels)) {
    stop("Distance matrix dimensions do not match number of IDs", call. = FALSE)
  }
  rownames(M) <- colnames(M) <- labels
  D <- as.dist(M)

  method <- tolower(opt$method)
  tr <- if (method == "upgma") {
    as.phylo(hclust(D, method = "average"))
  } else if (method == "nj") {
    nj(D)
  } else {
    stop("Unsupported tree method: ", opt$method, ". Use nj or upgma.", call. = FALSE)
  }

  if (!is.null(opt$outgroup) && opt$outgroup %in% tr$tip.label) {
    tr <- root(tr, outgroup = opt$outgroup, resolve.root = TRUE)
  } else if (!is.null(opt$outgroup)) {
    warning("Outgroup not found among tip labels: ", opt$outgroup)
  }

  write.tree(tr, paste0(opt$out, ".nwk"))
  inferred_tip_count <- length(tr$tip.label)
  requested <- read_sample_list(opt$tree_sample_list)
  if (!is.null(requested)) {
    selected <- sample_list_matches(requested, tr$tip.label, "TREE")
    if (length(selected) < 2) stop("TREE sample_list requires at least two matching samples", call. = FALSE)
    tr <- keep.tip(tr, selected)
    message("TREE display samples: ", length(selected), " / ", inferred_tip_count)
  }

  tip_df <- data.frame(label = tr$tip.label, display_label = tr$tip.label, stringsAsFactors = FALSE)
  clade_df <- NULL
  if (!is.null(opt$clades) && file.exists(opt$clades) && file.info(opt$clades)$size > 0) {
    clade_df <- read_clades_file(opt$clades, rename_df)
  }
  if (!is.null(clade_df)) {
    tip_df <- left_join(tip_df, clade_df, by = "label")
    tip_df$clade_label <- ifelse(is.na(tip_df$clade) | !nzchar(tip_df$clade), "Unclassified", tip_df$clade)
  }
  clade_palette_override <- read_clade_colors_file(opt$clade_colors)
  clade_colors <- NULL
  if (!is.null(clade_df) && "clade_label" %in% names(tip_df)) {
    # Build one named palette for all clade labels. This same object is used for
    # tree branches, tree tips, tree labels, PCA points, and LEA clade annotations.
    clade_colors <- make_clade_palette(c(tip_df$clade_label, clade_df$clade), palette = clade_palette_override)
    write_clade_palette(clade_colors, paste0(opt$out, ".clade_colors.tsv"))
  }
  node_ann <- make_tree_node_annotations(tr, tip_df)

  plot_tree <- function(filename) {
    if (has_ggtree) {
      p <- ggtree::ggtree(tr, layout = opt$layout, branch.length = branch_length_mode, size = 0.25)
      p$data <- left_join(p$data, node_ann, by = "node")

      if (!is.null(clade_df) && isTRUE(opt$color_branches)) {
        p <- p + aes(color = branch_clade)
      }

      if (!is.null(clade_df)) {
        p <- p +
          ggtree::geom_tippoint(aes(color = tip_clade), size = opt$tip_size, alpha = 0.9) +
          scale_color_manual(values = clade_colors, na.value = "grey70", drop = FALSE) +
          theme(legend.position = "right") +
          labs(color = "Clade")
      } else {
        p <- p +
          ggtree::geom_tippoint(size = opt$tip_size, alpha = 0.9) +
          theme(legend.position = "none")
      }

      if (isTRUE(opt$label_tips)) {
        if (!is.null(clade_df) && isTRUE(opt$color_tip_labels)) {
          p <- p + ggtree::geom_tiplab(aes(label = display_label, color = tip_clade), size = opt$label_size)
        } else {
          p <- p + ggtree::geom_tiplab(aes(label = display_label), size = opt$label_size, color = "black")
        }
      }
      ggsave(filename, p, width = opt$width, height = opt$height, dpi = opt$dpi, limitsize = FALSE)
    } else {
      if (grepl("[.]pdf$", filename)) pdf(filename, width = opt$width, height = opt$height)
      else png(filename, width = opt$width, height = opt$height, units = "in", res = opt$dpi)

      edge_cols <- rep("black", nrow(tr$edge))
      tip_cols <- rep("black", length(tr$tip.label))
      if (!is.null(clade_df) && !is.null(clade_colors)) {
        tip_clades <- node_ann$tip_clade[match(seq_along(tr$tip.label), node_ann$node)]
        if (isTRUE(opt$color_tip_labels)) {
          tip_cols <- unname(clade_colors[tip_clades])
          tip_cols[is.na(tip_cols)] <- "grey60"
        }
        if (isTRUE(opt$color_branches)) {
          edge_clades <- node_ann$branch_clade[match(tr$edge[, 2], node_ann$node)]
          edge_cols <- unname(clade_colors[edge_clades])
          edge_cols[is.na(edge_cols)] <- "grey70"
        }
      }

      plot(tr, type = ifelse(opt$layout %in% c("fan", "unrooted", "radial"), opt$layout, "phylogram"),
           use.edge.length = branch_length_enabled,
           cex = ifelse(isTRUE(opt$label_tips), opt$label_size / 2, 0.001),
           tip.color = tip_cols, edge.color = edge_cols, no.margin = TRUE)
      if (!is.null(clade_df) && !is.null(clade_colors)) {
        tip_clades <- node_ann$tip_clade[match(seq_along(tr$tip.label), node_ann$node)]
        tip_cols_points <- unname(clade_colors[tip_clades])
        tip_cols_points[is.na(tip_cols_points)] <- "grey60"
        ape::tiplabels(pch = 16, col = tip_cols_points, cex = opt$tip_size / 3)
      }
      dev.off()
    }
  }

  writeLines(c(
    paste0("tree_method	", method),
    paste0("tree_layout	", opt$layout),
    paste0("branch_length	", ifelse(branch_length_enabled, "yes", "no")),
    paste0("branch_length_mode	", branch_length_mode),
    paste0("label_tips	", ifelse(isTRUE(opt$label_tips), "yes", "no")),
    paste0("label_size	", opt$label_size),
    paste0("tip_size	", opt$tip_size),
    paste0("color_branches	", ifelse(isTRUE(opt$color_branches), "yes", "no")),
    paste0("color_tip_labels	", ifelse(isTRUE(opt$color_tip_labels), "yes", "no")),
    paste0("clade_colors_file	", ifelse(!is.null(clade_colors), paste0(opt$out, ".clade_colors.tsv"), "")),
    paste0("tree_plot_pdf	", paste0(tree_plot_prefix, ".pdf")),
    paste0("tree_plot_png	", paste0(tree_plot_prefix, ".png")),
    paste0("outgroup	", ifelse(is.null(opt$outgroup), "", opt$outgroup)),
    paste0("tips	", length(tr$tip.label)),
    paste0("inferred_tips\t", inferred_tip_count),
    paste0("sample_list\t", ifelse(is.null(opt$tree_sample_list), "", opt$tree_sample_list)),
    paste0("displayed_samples\t", paste(tr$tip.label, collapse = ",")),
    paste0("internal_nodes	", tr$Nnode)
  ), con = paste0(opt$out, ".tree_info.tsv"))

  plot_tree(paste0(tree_plot_prefix, ".pdf"))
  plot_tree(paste0(tree_plot_prefix, ".png"))
}

if (!is.null(opt$pca_eigenvec)) {
  eig <- read_plink_eigenvec(opt$pca_eigenvec)
  requested <- read_sample_list(opt$pca_sample_list)
  if (!is.null(requested)) {
    fid <- apply_rename(eig$FID, rename_df)
    iid <- apply_rename(eig$IID, rename_df)
    selected <- sample_list_matches(requested, unique(c(fid, iid)), "PCA")
    keep <- fid %in% selected | iid %in% selected
    message("PCA display samples: ", sum(keep), " / ", nrow(eig))
    eig <- eig[keep, , drop = FALSE]
  }

  if (!is.null(opt$pca_eigenval)) {
    evals <- scan(opt$pca_eigenval, quiet = TRUE)
    varprop <- evals / sum(evals) * 100
  } else {
    varprop <- rep(NA_real_, length(grep("^PC[0-9]+$", colnames(eig))))
  }

  if (!is.null(clade_df)) {
    key_info <- choose_pca_clade_key(eig, clade_df, rename_df)
    eig$name <- key_info$name
    eig <- left_join(eig, clade_df, by = c("name" = "label"))
    eig$clade_label <- ifelse(is.na(eig$clade) | !nzchar(eig$clade), "Unclassified", eig$clade)
    message(
      "PCA clade matches using ", key_info$key, ": ", key_info$matches, " / ", nrow(eig),
      " (FID matches: ", key_info$fid_matches, "; IID matches: ", key_info$iid_matches, ")"
    )
    clade_colors <- make_clade_palette(c(names(clade_colors), clade_df$clade, eig$clade_label), palette = clade_colors)
    write_clade_palette(clade_colors, paste0(opt$out, ".clade_colors.tsv"))
  } else {
    # No clades file: keep PCA uncoloured and label by the original FID-like field.
    eig$name <- apply_rename(trimws(eig$FID), rename_df)
  }

  pcx <- paste0("PC", opt$pca_pc_x)
  pcy <- paste0("PC", opt$pca_pc_y)
  if (!(pcx %in% colnames(eig)) || !(pcy %in% colnames(eig))) {
    stop("Requested PCA axes are not present in eigenvec file: ", pcx, " and/or ", pcy, call. = FALSE)
  }
  xlab <- if (!is.na(varprop[opt$pca_pc_x])) sprintf("%s (%.1f%%)", pcx, varprop[opt$pca_pc_x]) else pcx
  ylab <- if (!is.na(varprop[opt$pca_pc_y])) sprintf("%s (%.1f%%)", pcy, varprop[opt$pca_pc_y]) else pcy

  write.table(eig, paste0(opt$out, ".pca_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (!is.null(clade_df)) {
    p <- ggplot(eig, aes(x = .data[[pcx]], y = .data[[pcy]], color = clade_label)) +
      geom_point(size = opt$pca_point_size, alpha = opt$pca_alpha) +
      scale_color_manual(values = clade_colors, drop = FALSE) +
      theme_bw() +
      labs(x = xlab, y = ylab, color = "Clade")
  } else {
    p <- ggplot(eig, aes(x = .data[[pcx]], y = .data[[pcy]])) +
      geom_point(size = opt$pca_point_size, alpha = opt$pca_alpha) +
      theme_bw() +
      labs(x = xlab, y = ylab)
  }

  pca_labels <- if (is.null(opt$pca_label_points)) isTRUE(opt$label_tips) else as_bool_option(opt$pca_label_points, "--pca_label_points")
  pca_label_size <- if (is.null(opt$pca_label_size)) opt$label_size else opt$pca_label_size
  if (pca_labels && has_ggrepel) {
    p <- p + ggrepel::geom_text_repel(aes(label = name), size = pca_label_size)
  }
  ggsave(paste0(opt$out, ".pca.pdf"), p, width = 8, height = 6, dpi = opt$dpi, limitsize = FALSE)
  ggsave(paste0(opt$out, ".pca.png"), p, width = 8, height = 6, dpi = opt$dpi, limitsize = FALSE)
}

cat("Wrote:\n")
if (has_tree) {
  cat(paste0("  ", opt$out, ".nwk\n"))
  cat(paste0("  ", tree_plot_prefix, ".pdf\n"))
  cat(paste0("  ", tree_plot_prefix, ".png\n"))
  cat(paste0("  ", opt$out, ".tree_info.tsv\n"))
  if (!is.null(clade_colors)) cat(paste0("  ", opt$out, ".clade_colors.tsv\n"))
  cat(paste0("Tree method: ", method, "\n"))
  cat(paste0("Tree layout/type: ", opt$layout, "\n"))
  cat(paste0("Branch lengths used: ", ifelse(branch_length_enabled, "yes", "no"), "\n"))
  cat(paste0("Tip labels shown: ", ifelse(isTRUE(opt$label_tips), "yes", "no"), "\n"))
  cat(paste0("Label size: ", opt$label_size, "\n"))
  cat(paste0("Tip size: ", opt$tip_size, "\n"))
  cat(paste0("Branches coloured by clade: ", ifelse(isTRUE(opt$color_branches), "yes", "no"), "\n"))
  cat(paste0("Tip/accession labels coloured by clade: ", ifelse(isTRUE(opt$color_tip_labels), "yes", "no"), "\n"))
}
if (!is.null(opt$pca_eigenvec)) {
  cat(paste0("  ", opt$out, ".pca_scores.tsv\n"))
  cat(paste0("  ", opt$out, ".pca.pdf\n"))
  cat(paste0("  ", opt$out, ".pca.png\n"))
}
