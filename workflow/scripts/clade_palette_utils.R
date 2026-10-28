# Shared clade colour utilities for tree/PCA/LEA plots.
# This file is sourced by build_tree_postprocess.R and run_lea.R so that the
# same clade names always resolve to the same colours across all figures.

read_clade_colors_file <- function(path) {
  if (is.null(path) || !file.exists(path) || file.info(path)$size == 0) return(NULL)
  x <- try(read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, comment.char = ""), silent = TRUE)
  if (inherits(x, "try-error")) {
    x <- read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE,
                    check.names = FALSE, comment.char = "")
  }
  if (ncol(x) < 2 || nrow(x) == 0) return(NULL)
  x <- x[, 1:2, drop = FALSE]
  first <- tolower(trimws(as.character(unlist(x[1, ], use.names = FALSE))))
  if (length(first) >= 2 && first[1] %in% c("clade", "group", "population", "pop", "lineage") &&
      first[2] %in% c("color", "colour", "hex", "hex_color", "hex_colour")) {
    x <- x[-1, , drop = FALSE]
  }
  colnames(x) <- c("clade", "color")
  x$clade <- trimws(as.character(x$clade))
  x$color <- trimws(as.character(x$color))
  x <- x[!is.na(x$clade) & nzchar(x$clade) & !is.na(x$color) & nzchar(x$color), , drop = FALSE]
  if (nrow(x) == 0) return(NULL)
  x <- x[!duplicated(x$clade), , drop = FALSE]
  stats::setNames(x$color, x$clade)
}

make_clade_palette <- function(labels, palette = NULL) {
  labs <- unique(as.character(labels))
  labs <- labs[!is.na(labs) & nzchar(labs)]
  if (length(labs) == 0) return(NULL)

  ordered_labs <- c(sort(setdiff(labs, "Unclassified")), intersect("Unclassified", labs))
  pal <- stats::setNames(rep(NA_character_, length(ordered_labs)), ordered_labs)

  if (!is.null(palette)) {
    if (is.data.frame(palette)) {
      palette <- stats::setNames(as.character(palette[[2]]), as.character(palette[[1]]))
    }
    palette <- palette[!is.na(names(palette)) & nzchar(names(palette)) & !is.na(palette) & nzchar(palette)]
    shared <- intersect(names(pal), names(palette))
    if (length(shared) > 0) pal[shared] <- palette[shared]
  }

  missing_non_unclassified <- names(pal)[is.na(pal) & names(pal) != "Unclassified"]
  if (length(missing_non_unclassified) > 0) {
    pal[missing_non_unclassified] <- grDevices::hcl.colors(length(missing_non_unclassified), palette = "Dark 3")
  }
  if ("Unclassified" %in% names(pal) && is.na(pal["Unclassified"])) {
    pal["Unclassified"] <- "grey60"
  }

  pal[labs]
}

write_clade_palette <- function(palette, path) {
  if (is.null(palette) || length(palette) == 0 || is.null(path)) return(invisible(FALSE))
  x <- data.frame(
    clade = names(palette),
    color = unname(palette),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  x <- x[!is.na(x$clade) & nzchar(x$clade) & !duplicated(x$clade), , drop = FALSE]
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(TRUE)
}
