
#' @export
expr_name_ <- function(name) {
  name <- gsub(".txt", "", name)
  name <- gsub("_0.5.expr", "", name)
  name <- gsub("-", "_", name)
  
  name
}

normalize_v_ <- function(x) {
  m <- mean(x) # column mean
  std <- sd(x) # column (sample) sd
  x <- (x - m) / std
  
  x
}

#' Load predicted expression files fro ma folder
#' @export
load_expression <- function(folder, white_list=NULL, metaxcan_style=FALSE) {
  results <- list()
  names <- list.files(folder)
  
  for (name in names) {
    if (!is.null(white_list) && !grepl(white_list, name)) {
      next;
    }
    
    path <- file.path(folder,name)
    k <- read.delim(path)
    k <- apply(k, 2, normalize_v_)
    k <- data.frame(k)
    
    k <- k[, !apply(is.na(k), 2, all)]
    
    key <- expr_name_(name)
    if (metaxcan_style) {
      names(k) <- remove_id_from_ensemble(names(k))
      key <- gsub("TW_","",key)
    }
    results[[key]] <- k
  }
  
  results
}

#' @export
convert_expression_to_metaxcan <- function(expression) {
  result <- list()
  models <- names(expression)
  for (model in models) {
    a <- expression[[model]]
    names(a) <- remove_id_from_ensemble(names(a))
    
    key <- gsub("TW_", "", model)
    result[[key]] <- a
    
    expression[[model]] <- NULL
  }
  
  result
}

#'keep only those expressions you want
#' @export
select_expressions <- function(expression, models) {
  r <- list()
  for (model in models) {
    r[[model]] <- expression[[model]]
  }
  
  r
}

#'load plink fam file
#' @export
load_fam <- function(path, standardize=FALSE) {
  fam <- read.delim(path, sep="", header=FALSE)
  colnames(fam) <- c("fam", "id", "a", "b", "c", "pheno")
  if (standardize) {
    fam$pheno <- scale(fam$pheno)
  }
  
  fam
}

#' all genes in expression
#' @export
get_genes <- function(expression) {
  genes <- c()
  expression_names <- names(expression)
  for (name in expression_names) {
    k <- expression[[name]]
    genes <- c(genes, colnames(k))
  }
  genes <-unique(genes)
  
  genes
}

#' @export
get_cleared_tissues <- function(expression, gene) {
  tissues <- c()
  expression_names <- names(expression)
  for (name in expression_names) {
    k <- expression[[name]]
    if (gene %in% names(k)) {
      tissues <- c(tissues, name) 
    }
  }
  tissues <-unique(tissues)
  
  tissues
}

#' assorted data to do predixcan fo ra given gene
#' @export
get_data_for_gene <- function(expression, fam, gene) {
  results <- fam %>% dplyr::select(pheno)
  expression_names <- sort(names(expression))
  for (name in expression_names) {
    k <- expression[[name]]
    if (gene %in% colnames(k)) {
      e <- k %>% dplyr::select(dplyr::contains(gene))
      e <- e %>% dplyr::rename_(.dots=setNames(gene, name))
      results <- cbind(results, e)
    }
  }
  
  results
}

#' get covariance fro maforementioend data
#' @export
get_cov <- function(data) {
  k <- data %>% dplyr::select(-pheno)
  models <- colnames(k)
  M <- c()
  for (model in models) {
    m <- data[[model]]
    M <- cbind(M,m)
  }
  colnames(M) <- models
  covariance <- cov(M)
  
  covariance
}

#' alternative method to get covariance from expression, for a certain gene, wit ha certain orde rof models
#' @export
get_cov_2 <- function(expression, gene, models){
  M <- c()
  for (model in models) {
    d <- expression[[model]]
    M <- cbind(M, d[[gene]])
  }
  colnames(M) <- models
  covariance <- cov(M)
  
  covariance
}

###############################################################################
#' @export
get_predixcan_association_from_data <- function(data, binomial = FALSE) {
  k <- data %>% dplyr::select(-pheno)
  models <- colnames(k)
  y <- data$pheno
  results <- data.frame()
  for (m in models) {
    x <- data[[m]]
    coefs <- if (binomial) {
      fit <- glm(y ~x, family="binomial", data=data, maxit=10)
      coef(summary(fit))[-1,]
    } else {
      fit <- lm(y ~ x -1)
      coef(summary(fit))
    }
    r <- data.frame(model=m, zscore=coefs[3], pvalue=coefs[4])
    results <- rbind(results, r)
  }
  
  results
}

#' @export
get_predixcan_association <- function(expression, fam, binomial = FALSE) {
  models <- names(expression)
  results <- data.frame()
  genes <- get_genes(expression)
  for (gene in genes) {
    data <- get_data_for_gene(expression, fam, gene)
    r <- get_predixcan_association_from_data(data, binomial)
    r$gene <- gene
    results <- rbind(results, r)
  }
  
  results
}

###############################################################################
# multivariate chi-square wald test

get_wald_p_ <- function(covariance, zscores, tol = sqrt(.Machine$double.eps), epsilon=NULL) {
  if (nrow(covariance) > 1 && !is.null(epsilon)) {
    e <- diag(rep(epsilon, nrow(covariance)))
    covariance <- covariance + e
  }
  inv <- MASS::ginv(covariance, tol)
  z<- zscores$zscore
  w_2 <- t(z) %*% (inv %*% z)
  w_2 <- w_2[1][1]
  p_w_2 <- pchisq(w_2, length(z), lower.tail = FALSE)
  
  p_w_2
}

###############################################################################
# estimate cross model from marginal analysis

get_extended_wald_p <- function(data, tol = sqrt(.Machine$double.eps)) {
  k <- data %>% dplyr::select(-pheno)
  models <- colnames(k)
  covariance <- get_cov(data, models)
  predixcan <- get_predixcan_association_from_data(data)
  p_w_2 <- get_wald_p_(covariance, predixcan, tol)
  
  p_w_2
}

get_wald_stats_for_gene_ <- function(expression, predixcan, the_gene, tol = sqrt(.Machine$double.eps), epsilon=NULL) {
  cleared_tissues <- get_cleared_tissues(expression, the_gene)
  if (length(cleared_tissues) < 2){
    return(NULL)
  }
  p <- predixcan %>% dplyr::filter(gene == the_gene)
  cov <- get_cov_2(expression, the_gene, p$model)
  p_w <- get_wald_p_(cov, p, tol=tol, epsilon=epsilon)
  p_i_best <- min(p$pvalue)
  p_i_worst <- max(p$pvalue)
  
  list("w"=p_w, "n"=nrow(cov), "p_i_best"=p_i_best, "p_i_worst"=p_i_worst)
}

get_wald_stats_for_gene_m_ <- function(expression, metaxcan, the_gene, tol = sqrt(.Machine$double.eps), epsilon=NULL) {
  cleared_tissues <- get_cleared_tissues(expression, the_gene)
  if (length(cleared_tissues) < 2){
    return(NULL)
  }
  
  m <- metaxcan %>% dplyr::filter(gene == the_gene, tissue %in% cleared_tissues)
  if (nrow(m) == 0) {
    return(NULL)
  }
  
  m$pvalue <- 2*pnorm(-abs(m$zscore))
  m$model <- m$tissue
  cov <- get_cov_2(expression, the_gene, m$model)
  p_w <- get_wald_p_(cov, m, tol=tol, epsilon=epsilon)
  p_i_best <- min(m$pvalue)
  p_i_worst <- max(m$pvalue)
  
  list("w"=p_w, "n"=nrow(cov), "p_i_best"=p_i_best, "p_i_worst"=p_i_worst)
}

#' @export
get_multi_predixcan_estimate <- function(expression, predixcan, tol = sqrt(.Machine$double.eps), epsilon=NULL ) {
  genes <- unique(predixcan$gene)
  results <- data.frame()
  for (gene in genes) {
    l <- get_wald_stats_for_gene_(expression, predixcan, gene, tol, epsilon)
    if (is.null(l)) {
      next;
    }
    r <- data.frame(gene=gene, pvalue = l$w, n = l$n, p_i_best=l$p_i_best, p_i_worst=l$p_i_worst)
    results <- rbind(results, r)
  }
  
  results
}

#' @export
get_multi_metaxcan_estimate <- function(expression, metaxcan, tol = sqrt(.Machine$double.eps), epsilon=NULL ) {
  genes <- unique(metaxcan$gene)
  results <- data.frame()
  for (gene in genes) {
    l <- get_wald_stats_for_gene_m_(expression, metaxcan, gene, tol, epsilon)
    if (is.null(l)) {
      next;
    }
    r <- data.frame(gene=gene, pvalue = l$w, n = l$n, p_i_best=l$p_i_best, p_i_worst=l$p_i_worst)
    results <- rbind(results, r)
  }
  
  results
}

###############################################################################

get_multi_fit <- function(data, binomial=FALSE) {
  f <- colnames(data %>% dplyr::select(-pheno))
  
  formula <- "pheno ~ "
  formula <- paste0(formula, paste0(f, collapse=" + "))
  if (!binomial) { formula <- paste0(formula, " - 1") }
  
  fit <- if(binomial) {
    glm(formula,family="binomial",data=data,maxit=10)
  } else {
    lm(formula, data) 
  }
  
  fit
}

get_stats_from_multi_fit <- function(fit, gene) {
  value <- summary(fit)$fstatistic[["value"]]
  df1 <- summary(fit)$fstatistic[["numdf"]]
  df2<- summary(fit)$fstatistic[["dendf"]]
  p <- pf(value, df1, df2, lower.tail=FALSE)
  
  pvalues <- unname(coef(summary(fit))[,4])
  p_best <- min(pvalues)
  p_worst <- max(pvalues)
  
  r <- data.frame(
    gene = gene,
    pvalue = p,
    p_i_best = p_best,
    p_i_worst = p_worst
  )
  
  r
}

get_stats_from_multi_fit_binomial <- function(fit, gene) {
  #the first coefficient is the intercept for binomial, so we must lose it (-1 in the indexing)
  p <- pchisq(fit$null.deviance-fit$deviance, fit$df.null - fit$df.residual, lower.tail = FALSE )
  
  pvalues <- unname(summary(fit)$coefficients[-1,4])
  p_best <- min(pvalues)
  p_worst <- max(pvalues)
  
  r <- data.frame(
    gene = gene,
    pvalue = p,
    p_i_best = p_best,
    p_i_worst = p_worst
  )
  
  r
}

#' @export
get_multivariate_predixcan_association <- function(expression, fam, binomial=FALSE) {
  results <- data.frame()
  genes <- get_genes(expression)
  for (gene in genes) {
    data <- get_data_for_gene(expression, fam, gene)
    r <- if (binomial) {
      fit <- get_multi_fit(data, binomial=TRUE)
      get_stats_from_multi_fit_binomial(fit, gene)
    } else {
      fit <- get_multi_fit(data)
      get_stats_from_multi_fit(fit, gene)
    }
    results <- rbind(results, r)
  }
  results <- results %>% dplyr::arrange(pvalue)
  
  results
}

###############################################################################
# High level driver

#' @export
get_multi_predixcan_estimate_progression <- function(expression, predixcan, estimate=get_multi_predixcan_estimate, epsilon=c(0), tol=c(0)){
  c_ <- data.frame()
  for (e_ in epsilon) {
    for (t_ in tol) {
      c_a <- estimate(expression, predixcan, epsilon=e_, tol=t_)
      c_a$epsilon <- e_
      c_a$tolerance <- t_
      c_ <- rbind(c_, c_a)
    }
  }
  
  c_
}
