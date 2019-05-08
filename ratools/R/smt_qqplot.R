
###############################################################################
#' @export
get_bonferroni_significant_genes_ <- function(data, x=2.8, y=1.5) {
  config <- get_sa_config_("")
  qq_data <- prepare_qq_data_(data, config$threshold)
  result <- qq_data$result
  decoration <- qq_data$decoration
  nb <- nrow(result %>% dplyr::filter(y > decoration$b))
  d <- data.frame(x=x, y=y, label=paste0("# bonferroni genes:", nb), colour=config$colours[1], stringsAsFactors = FALSE)
}

#' @export
get_bonferroni_significant_genes_sa_ <- function(data, selected_tissue, x=2.8, y=0) {

  if(!("pvalue" %in% colnames(data))) {
    data$pvalue <- 2*pnorm(-abs(data$zscore))
  }
  config <- get_sa_config_(selected_tissue)
  m <- get_sa_(data, selected_tissue)
  qq_data <- prepare_qq_faceted_data_(m, config$facet_column, threshold=config$threshold, facet_order=config$facet_order)
  decoration <- qq_data$decoration


  get_n <- function(data, decoration, label) {
    bt <- (decoration %>% dplyr::filter(facet_w == label))$b
    f <- data %>% dplyr::filter(label_w == label, -log10(pvalue)>bt)
    u <- unique(f$gene)
    length(u)
  }
  bs <- get_n(m, decoration, selected_tissue)
  ba <- get_n(m, decoration, "All tissues")
  facet <- factor(c(selected_tissue, "All tissues"), c(selected_tissue, "All tissues"))
  data.frame(x=x, y=y, label=paste0("# bonferroni genes:", c(bs, ba)), facet_w=facet, colour=config$colours[2], stringsAsFactors = FALSE)
}

#' @export
get_sa_ <- function(marginal, selected_tissue) {
  m1 <- marginal %>% dplyr::filter(tissue == selected_tissue) %>% dplyr::mutate(label_w = selected_tissue)
  m2 <- marginal %>% dplyr::mutate(label_w = "All tissues")
  rbind(m1,m2)
}

#' @export
get_sa_config_ <- function(selected_tissue, labels=c("Multi-Tissue PrediXcan", "Univariate PrediXcan")) {
  qq_plot_config(point_size=2, columns=2, scales="fixed",
                 facet_column="label_w",
                 facet_order=c(selected_tissue, "All tissues"),
                 labels=labels)
}

#' @export
multi_predixcan_vs_marginalsa_qqplot <- function(multi_estimate,
                                                  marginal,
                                                  selected_tissue,
                                                  labels=c("Multi-Tissue PrediXcan", "Univariate PrediXcan")) {
  m <- get_sa_(marginal, selected_tissue)
  config <- get_sa_config_(selected_tissue, labels)
  qq_data_ <- prepare_qq_data_(multi_estimate, config$threshold)
  qq_additional_data_ <- prepare_qq_faceted_data_(m, config$facet_column, threshold=config$threshold, facet_order=config$facet_order)
  build_qq_plot_(qq_data_, qq_additional_data_, config)
}


###############################################################################
