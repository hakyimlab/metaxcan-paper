
#' @export
run_load_expression_ <- function(expression_folder, verbose=FALSE) {
  if (verbose) {
    print("loading expression")
    start.time <- Sys.time()
  }
  
  expression <- load_expression(expression_folder)
  
  if (verbose) {
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  }
  
  expression
}

#' @export
run_get_multi_predixcan_association_ <- function(expression, fam, binomial, verbose=FALSE) {
  if (verbose) {
    log <- "Calculating multi PrediXcan"
    if (binomial) log <- paste0(log, " (binomial)")
    print(log)
    start.time <- Sys.time()
  }
  
  multi_predixcan <- get_multivariate_predixcan_association(expression, fam, binomial)
  
  if (verbose) {
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  }
  
  multi_predixcan
}

#' @export
run_get_predixcan_association_ <- function(expression, fam, binomial, verbose=FALSE) {
  if (verbose) {
    log <- "Calculating Marginal PrediXcan"
    if (binomial) log <- paste0(log, " (binomial)")
    print(log)
    start.time <- Sys.time()
  }
  
  predixcan <- get_predixcan_association(expression, fam, binomial)
  
  if (verbose) {
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  }

  predixcan
}

#' @export
run_get_multi_predixcan_estimate_progression_ <- 
  function(expression, predixcan, progression_tol, verbose=FALSE) {
  if (verbose) {
    print("Calculating multi PrediXcan Estimate")
    start.time <- Sys.time() 
  }
  
  multi_predixcan <- get_multi_predixcan_estimate_progression(expression, predixcan, tol=progression_tol)
  
  if (verbose) {
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  }
  
  multi_predixcan
}


#' @export
run_get_multi_predixcan_estimate_progression_m_ <- 
  function(expression, metaxcan, progression_tol, verbose=FALSE) {
  if (verbose) {
    print("Calculating multi PrediXcan Estimate from metaxcan")
    start.time <- Sys.time()
  }
  
  multi_predixcan_s <- get_multi_predixcan_estimate_progression(
    expression, metaxcan, estimate=get_multi_metaxcan_estimate, tol=progression_tol)
  
  if (verbose) {
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  }
  
  multi_predixcan_s
}

#' @export
run_cpm_multi_predixcan <- function(expression_folder, fam, results, binomial=FALSE, verbose=TRUE) {
  expression <- run_load_expression_(expression_folder, verbose)
  fam <- load_fam(fam)
  multi_predixcan <- run_get_multi_predixcan_association_(expression, fam, binomial, verbose)
  ra_save_delim(multi_predixcan, results)
}

#' @export
run_cpm_predixcan <- 
  function(expression_folder, fam, predixcan_results, multi_predixcan_estimate_results, 
           progression_tol, binomial=FALSE, verbose=TRUE) {
  expression <- run_load_expression_(expression_folder, verbose)
  fam <- load_fam(fam)
  predixcan <- run_get_predixcan_association_(expression, fam, binomial, verbose)
  multi_predixcan_estimate <- run_get_multi_predixcan_estimate_progression_(expression, predixcan, progression_tol, verbose=verbose)
  ra_save_delim(predixcan, predixcan_results)
  ra_save_delim(multi_predixcan_estimate, multi_predixcan_estimate_results)
}

#' @export
run_cpm_multi_predixcan_estimate_summary <- 
  function(expression_folder, metaxcan_folder, estimate_results, 
           progression_tol, dichotomy=FALSE, verbose=TRUE) {
  expression <- run_load_expression_(expression_folder, verbose)
  expression <- convert_expression_to_metaxcan(expression)
  metaxcan <-  load_metaxcan_folder(metaxcan_folder, ".csv")
  multi_metaxcan_estimate <- run_get_multi_predixcan_estimate_progression_m_(expression, metaxcan, progression_tol, verbose=verbose)
  ra_save_delim(multi_metaxcan_estimate, estimate_results)
}

#' @export
run_cpm <- function(expression_folder, fam_path, metaxcan_folder, results_dir, binomial, progression_tol, verbose=TRUE) {
  if(!dir.exists(results_dir)) dir.create(results_dir)
  
  expression <- run_load_expression_(expression_folder, verbose)
  fam <- load_fam(fam_path)
  
  #
  multi_predixcan <- run_get_multi_predixcan_association_(expression, fam, binomial, verbose)
  path <- file.path(results_dir, "multi_predixcan.txt")
  ra_save_delim(multi_predixcan, path)
  
  #
  predixcan <- run_get_predixcan_association_(expression, fam, binomial, verbose)
  multi_predixcan_estimate <- run_get_multi_predixcan_estimate_progression_(expression, predixcan, progression_tol, verbose=verbose)
  path <- file.path(results_dir, "predixcan.txt")
  ra_save_delim(predixcan, path)
  path <- file.path(results_dir, "multi_predixcan_estimate.txt")
  ra_save_delim(multi_predixcan_estimate, path)
  
  #
  expression <- convert_expression_to_metaxcan(expression)
  metaxcan <-  load_metaxcan_folder(metaxcan_folder, ".csv")
  multi_metaxcan_estimate <- run_get_multi_predixcan_estimate_progression_m_(expression, metaxcan, progression_tol, verbose=verbose)
  path <- file.path(results_dir, "multi_summary_predixcan_estimate_r.txt")
  ra_save_delim(multi_metaxcan_estimate, path)
}