
#' @export
pheno_to_name_suffix <- function(pheno) {
  name <- gsub("_", " ",pheno)
  name <- toupper(name)
  name
}

#' @export
to_file_suffix <- function(name ) {
  name <- gsub("\\.", "_", name)
  return(name)
}

#' @export
build_image_name <- function(..., extension=".png") {
  name <- c(...)
  name <- paste(name, collapse='')
  name <- to_file_suffix(name)
  name <- paste0(name, ".png")
  name <- gsub(" ", "_", name)
  name <- gsub(",", "", name)
  name <- gsub("=", "_", name)
  return(name)
}

#' @export
wrap_sentence <- function(x, ...) {
  paste(strwrap(x, ...), collapse = "\n")
}

#' @export
parse_multi_predixcan_label <- function(name, number_name) {
  #number_regexp <- "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"
  scientific_regexp <- "[-+]?[\\d.]+(?:e[-+]?\\d+)?"
  scientific_regexp_end <- paste0("\\_", scientific_regexp, "\\.txt$")
  result <- list()
  if (grepl(scientific_regexp_end, name, perl=TRUE)) {
    part <- stringr::str_extract(name, scientific_regexp_end)
    part <- gsub(".txt", "", part)
    part <- gsub("_", "", part)
    number <- as.double(part)
    result[[number_name]] <- number
    result[["pheno"]] <- stringr::str_replace(name, scientific_regexp_end, "")
  }

  result
}

#' @export
load_multi_summary_predixcan_estimate <- function(path, number_name) {
  name <- split_path(path)[1]
  parts <- parse_multi_predixcan_label(name, number_name)
  k <- read.delim(path)
  k$pheno <- parts$pheno
  k[[number_name]] <- parts[[number_name]]

  k
}

#' @export
load_multi_summary_predixcan_estimate_folder <- function(folder, number_name, label) {
  results <- data.frame()
  names <- list.files(folder)
  for (name in names) {
    path <- file.path(folder,name)
    k <- load_multi_summary_predixcan_estimate (path, number_name)
    k$label <- label
    results <- rbind(results, k)
  }

  results
}

#' @export
load_multi_predixcan_estimate <- function(file_path, column_for_facet=NULL) {
  d <- read.delim(file_path)
  d
}

convert_predixcan_style_to_metaxcan_style <- function(predixcan) {
  p <- data.frame(predixcan)
  p$gene <- remove_id_from_ensemble(p$gene)
  if ("model"%in% names(p)) {
    p$tissue <- gsub("TW_","",p$model)
  }

  p
}

supplement_cross <- function(cross, factor_prefix=NULL) {
  p <- data.frame(cross)
  names <- names(p)
  if (! "facet_w" %in% names) {
    if ("cer" %in% names) {
      prefix <- if(is.null(factor_prefix)) "cer: " else factor_prefix
      p$facet_w <- build_levels(p$cer, prefix)
    } else if ("cthr" %in% names) {
      prefix <- if(is.null(factor_prefix)) "cthr: " else factor_prefix
      p$facet_w <- build_levels(p$cthr, "cthr: ")
    } else if ("ctr" %in% names) {
      prefix <- if(is.null(factor_prefix)) "ctr: " else factor_prefix
      p$facet_w <- build_levels(p$ctr, "ctr: ")
    } else {
      stop("could not build facet")
    }
  }

  p
}

#' @export
multi_predixcan_metaxcan_merge <- function(cross_predixcan, cross_metaxcan, full_outer_join=FALSE, factor_prefix=NULL){
  c_ <- convert_predixcan_style_to_metaxcan_style(cross_predixcan)
  cm_ <- supplement_cross(cross_metaxcan)
  if (full_outer_join) {
    r <- c_ %>% dplyr::full_join(cm_, by="gene")
  } else {
    r <- c_ %>% dplyr::inner_join(cm_, by="gene")
  }

  r
}

#' @export
multi_predixcan_vs_estimate_plot_config <- function(columns=NULL, scales="free", truncate = NULL, truncate_to_x=FALSE,
                                                    metaxcan_merge=FALSE, no_color=FALSE, facet_prefix=NULL, point_size=2) {
  config <- list()
  config$columns <- columns
  config$scales <- scales
  config$truncate <- truncate
  config$truncate_to_x <- truncate_to_x
  config$metaxcan_merge <- metaxcan_merge
  config$no_color <- no_color
  config$facet_prefix <- facet_prefix
  config$point_size <- point_size
  config
}

#' @export
multi_predixcan_vs_estimate_plot <- function(cross_predixcan, cross_estimate, config=multi_predixcan_vs_estimate_plot_config()) {

  m <- if (config$metaxcan_merge) {
    multi_predixcan_metaxcan_merge(cross_predixcan, cross_estimate, factor_prefix=config$facet_prefix)
  } else {
    cross_estimate <- supplement_cross(cross_estimate, factor_prefix=config$facet_prefix)
    cross_predixcan %>% dplyr::inner_join(cross_estimate, by="gene")
  }

  m$x <- -log10(m$pvalue.x)
  m$y <- -log10(m$pvalue.y)

  if (!is.null(config$truncate)) {
    m$x <- pmin(config$truncate, m$x)
    m$y <- pmin(config$truncate, m$y)
  }
  if (config$truncate_to_x) {
    k <- max(m$x[is.finite(m$x)])
    m$y <- pmin(m$y, k)
  }

  map <-ggplot2::aes(x=x, y=y)
  p <- ggplot2::ggplot(data = m) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(mapping=map, alpha=1, size=config$point_size) +
    ggplot2::geom_abline(intercept=0, slope=1, colour='red')

  if("facet_w" %in% colnames(m) && !is.null(config$columns)) {
    p <- p + ggplot2::facet_wrap(~facet_w, scales=config$scales,ncol=config$columns)
  }
  
  p
}
