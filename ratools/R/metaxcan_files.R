

#' Loads metaxcan result file
#'
#' @param path filesystem path to the file
#' @param finite keep only finite results
#' @return a dataframe with metaxcan results
#' @keywords metaxcan
#' @export
load_metaxcan_result_file <- function(path, finite=TRUE) {
  d <- read.csv(path)
  if (finite){
    d <- d %>% dplyr::filter(is.finite(zscore))
  }

  d
}

#' Heuristic to get phenotype name and model name from a metaxcan results file name
#'
#' @param text String to be parsed
#' @param remove_postfix String to be removed from the input text before parsing
#' @return list with phenotype and model
#' @keywords metaxcan
#' @export
pheno_model_from_file_name <- function(text, remove_postfix=NULL) {
  if (!is.null(remove_postfix)) {
    text <- gsub(remove_postfix, "", text)
  }
  pheno <- NULL
  model <- NULL
  if (grepl("_DGN", text)) {
    pheno <- strsplit(text, "_DGN")[[1]][1]
    model <- "DGN_WB"
  } else if (grepl("-TW_", text)) {
    comps = strsplit(text, "-TW_")
    pheno <- comps[[1]][1]
    model <- comps[[1]][2]
  } else if (grepl("_TW_", text)) {
    comps = strsplit(text, "_TW_")
    pheno <- comps[[1]][1]
    model <- comps[[1]][2]
  } else if (grepl("_eQTL_", text)) {
    comps = strsplit(text, "_eQTL_")
    pheno <- comps[[1]][1]
    model <- comps[[1]][2]
  }else if (grepl("_gEUVADIS", text)) {
    comps <- strsplit(text, "_gEUVADIS")
    pheno <- comps[[1]][1]
    model <- paste0("Geuvadis",comps[[1]][2])
  }

  list(pheno=pheno,model=model)
}

#' Loads all metaxcan results in a folder
#'
#' @param folder path where the result files are
#' @param remove_postfix pattern to be removed from file names before their parsing
#' @param white_list List of patterns for selecting files to be loaded
#' @param black_list List of patterns for discarding files not to be loaded
#' @param s List of columns from the results file to keep
#' @return Dataframe with metaxcan results
#' @keywords metaxcan
#' @export
load_metaxcan_result_folder_3 <- function(folder, remove_postfix=NULL, white_list=NULL, black_list=NULL,
                            s=c("gene","zscore","pheno", "model")) {
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
    path <- file.path(folder,name)
    d <- load_metaxcan_result_file(path)
    parts <- pheno_model_from_file_name(name, remove_postfix)
    d$pheno <- parts$pheno
    d$model <- parts$model
    d <- d %>% dplyr::select_(.dots=s)
    r <- rbind(r, d)
  }

  r
}

