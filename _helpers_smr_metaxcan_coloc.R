

get_coloc_proportions_2_ <- function(data) {
  data %>%
    summarize(n=n(), p3=sum(p_h3>0.5, na.rm=T), p4=sum(p_h4>0.5, na.rm=T)) %>%
    mutate(prop_p4 = round(p4/n,2), prop_p3 = round(p3/n,2), prop_p4_ = round(p4/(p4+p3),2), prop_p3_ = round(p3/(p3+p4), 2) )
}

get_smr_available_ <- function(data, exclude_genes) {
  bf <- 0.05/nrow(data)
  data %>% dplyr::filter(!(gene_name %in% exclude_genes), p_smr < bf)
}

get_smr_available_heidi_ <- function(data, exclude_genes) {
  bf <- 0.05/nrow(data)
  data %>% dplyr::filter(!(gene_name %in% exclude_genes), p_smr < bf, p_heidi > 0.05)
}

heidi_filter_ <- function(data) {
  data %>% dplyr::filter(p_heidi > 0.05)
}

smr_predixcan_coloc_fraction_data_process_ <- function(data, exclude_genes, additional_filter=NULL) {
  ###############################################################################
  smr_available <- data %>% get_smr_available_(exclude_genes)
  if (!is.null(additional_filter)) { smr_available <- additional_filter(smr_available) }
  if (nrow(smr_available) == 0) {
    return(NULL)
  }
  
  pheno <- data$phenotype %>% unique()
  
  max_smr_p <- max(smr_available$p_smr)
  
  result_smr_ <- smr_available %>% 
    get_coloc_proportions_2_() %>% mutate(phenotype=pheno, method='smr', genes="no_hla")
  
  ###############################################################################
  spredixcan_ <- data %>% dplyr::filter(!(gene_name %in% exclude_genes), pred_perf_pval<0.05, pval<max_smr_p)
  #spredixcan_ <- spredixcan_ %>% dplyr::filter(p_h3+p_h4>0.5)
  #if (!is.null(additional_filter)) { spredixcan_ <- additional_filter(spredixcan_) }
  result_predixcan_ <- spredixcan_ %>% 
    get_coloc_proportions_2_() %>% mutate(phenotype=pheno, method='predixcan', genes="no_hla")
  
  rbind(result_smr_, result_predixcan_)
}


smr_predixcan_coloc_fraction_data <- function(connection, phenos, gencode, verbose=FALSE, additional_filter=NULL) {
  #exclude hla genes
  HLA <- select_around_genes_2(gencode, chromosome="chr6", start=28866528, end=33775446)
  exclude_genes <- unique(HLA$gene_name)
  
  results <- data.frame()
  for (pheno in sort(phenos)) {
    if (verbose) { message("Processing " %&% pheno) }
    data <- build_data_2(connection, c(pheno))
    results_ <- smr_predixcan_coloc_fraction_data_process_(data, exclude_genes, additional_filter)
    if (is.null(results_)) { next }
    results <- rbind(results, results_)
  }
  
  results
}

smr_predixcan_coloc_fraction_plot_ <- function(data, ylabel=NULL) {
  p <- data %>%
    ggplot2::ggplot(ggplot2::aes(x=predixcan, y=smr, label=phenotype)) +
    ggplot2::geom_point(size = 5, color = 'grey35') + 
    ggplot2::geom_abline(slope=1, intercept=0) + 
    ggrepel::geom_label_repel(fontface = 'bold', size=6, segment.color = 'grey20', box.padding = grid::unit(0.5, "lines")) +
    ggplot2::coord_cartesian(xlim=c(0.0, 1.0), ylim=c(0.0, 1.0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust=0, size=28, face="bold"),
                   plot.subtitle= ggplot2::element_text(size=20, hjust=0.5, face="italic", color="black"),
                   axis.text = ggplot2::element_text(vjust = 0.5, hjust = 0.5, lineheight = 0,  size = 20,  face="bold"),
                   axis.title.x = ggplot2::element_text(size = 22, vjust=1,  face="bold"),
                   axis.title.y = ggplot2::element_text(size = 22, vjust=1,  face="bold"))
  if (!is.null(ylabel)) { p <- p + ggplot2::ylab(ylabel) }
  p
}

smr_predixcan_p3_fraction_plot <- function(data, ylabel=NULL) {
  data %>%
    dplyr::select(phenotype, prop_p3, method) %>% tidyr::spread(key=method, value=prop_p3) %>% 
    smr_predixcan_coloc_fraction_plot_(ylabel) +
    ggplot2::ggtitle("Proportion of Non-Colocalization:\n SMR vs PrediXcan",subtitle="\n#independent/#total: independent:P3>0.5; HLA region excluded") 
}

smr_predixcan_p3_b_fraction_plot <- function(data, ylabel=NULL) {
  data %>%
    dplyr::select(phenotype, prop_p3_, method) %>% tidyr::spread(key=method, value=prop_p3_) %>% 
    smr_predixcan_coloc_fraction_plot_(ylabel) +
    ggplot2::ggtitle("Proportion of Non-Colocalization:\n SMR vs PrediXcan",subtitle="\n#independent/(#independent + #colocalized); HLA region excluded") 
}

smr_predixcan_p4_fraction_plot <- function(data, ylabel=NULL) {
  data %>%
    dplyr::select(phenotype, prop_p4, method) %>% tidyr::spread(key=method, value=prop_p4) %>% 
    smr_predixcan_coloc_fraction_plot_(ylabel) +
    ggplot2::ggtitle("Proportion of Colocalization:\n SMR vs PrediXcan",subtitle="\n#colocalized/#total: colocalized:P4>0.5; HLA region excluded") 
}

smr_predixcan_p4_b_fraction_plot <- function(data, ylabel=NULL) {
  data %>%
    dplyr::select(phenotype, prop_p4_, method) %>% tidyr::spread(key=method, value=prop_p4_) %>% 
    smr_predixcan_coloc_fraction_plot_(ylabel) +
    ggplot2::ggtitle("Proportion of Colocalization:\n SMR vs PrediXcan",subtitle="\n#colocalized/(#independent + #colocalized) ; HLA region excluded") 
}