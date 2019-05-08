
#' Save as a tab-delimited file
#' @export
ra_save_delim <- function(data, path){
  write.table(data, file=path, row.names = FALSE, sep="\t", quote=FALSE)
}

#' Removes everything after a "."character, meant to remove id number of gencode names.
#' (An ugly convention)
#' 
#' @param v Strings that will be transformed.
#' @return Transformed strings
#' @keywords gencode
#' @export
remove_id_from_ensemble <- function(v) {
  k <- gsub("\\.(.*)", "", v)
  k
}

#' Checks if text matches anything in a list of patterns
#' (An ugly convention)
#'
#' @param v Strings that will be transformed.
#' @return Transformed strings
#' @keywords gencode
#' @export
is_listed <- function(text, list=NULL) {
  is <- FALSE;
  for (w in list) {
    if (grepl(w, text)) {
      is <- TRUE
      break;
    }
  }
  is
}

#' @export
load_metaxcan_result_file <- function(path) {
  d <- read.csv(path)
  d <- d %>% dplyr::filter(is.finite(zscore))
  d
}

#' @export
load_metaxcan_folder <- function(folder, remove_postfix=NULL, white_list=NULL, black_list=NULL,
                            s=c("gene","zscore","pheno", "tissue"), verbose=FALSE) {
  names <- list.files(folder)
  names <- sort(names)
  r <- data.frame()
  for(i in 1:length(names)) {
    name <- names[i]
    if (!is.null(white_list)) {
      skip <- !is_listed(name, white_list)
      if (skip) { next; }
    }
    if (!is.null(black_list)) {
      skip <- is_listed(name, black_list)
      if (skip) { next; }
    }
    if (verbose) {
      print(paste0("Loading ", name))
    }
    path <- file.path(folder,name)
    d <- load_metaxcan_result_file(path)
    parts <- pheno_tissue_from_file_name(name, remove_postfix)
    d$pheno <- parts$pheno
    d$tissue <- parts$tissue
    d <- d %>% dplyr::select(s)
    r <- rbind(r, d)
  }

  r
}

#' @export
pheno_tissue_from_file_name <- function(text, remove_postfix=NULL) {
  if (!is.null(remove_postfix)) {
    text <- gsub(remove_postfix, "", text)
  }
  pheno <- NULL
  tissue <- NULL
  if (grepl("_DGN", text)) {
    pheno <- strsplit(text, "_DGN")[[1]][1]
    tissue <- "DGN_WB"
  } else if (grepl("-TW_", text)) {
    comps = strsplit(text, "-TW_")
    pheno <- comps[[1]][1]
    tissue <- comps[[1]][2]
  } else if (grepl("_TW_", text)) {
    comps = strsplit(text, "_TW_")
    pheno <- comps[[1]][1]
    tissue <- comps[[1]][2]
  } else if (grepl("_eQTL_", text)) {
    comps = strsplit(text, "_eQTL_")
    pheno <- comps[[1]][1]
    tissue <- comps[[1]][2]
  }else if (grepl("_gEUVADIS", text)) {
    comps <- strsplit(text, "_gEUVADIS")
    pheno <- comps[[1]][1]
    tissue <- paste0("Geuvadis",comps[[1]][2])
  }

  list(pheno=pheno,tissue=tissue)
}

#' @export
build_levels <- function(figures, prefix) {
  l <- unique(figures)
  l <- sort(l[!is.na(l) & l!=0], decreasing = TRUE)
  l <- c(0, l)
  l <- paste0(prefix, l)
  level <- factor(paste0(prefix, figures), levels=l)

  level
}

#' @export
split_path <- function(x) {
  if (dirname(x)==x) 
    return(x)
  else 
    return(c(basename(x),split_path(dirname(x))))
}

#' @export
save_plot <- function(plot, path, height, width, res=NA) {
  png(path, height=height, width=width, res=res)
  print(plot)
  dev.off()
}

#' @export
save_plot_svg <- function(plot, path, height, width) {
  svg(path, height=height, width=width)
  print(plot)
  dev.off()
}
