

get_coloc_proportions_ <- function(data) {
  data %>% 
    summarize(n=n(), p3=sum(p_h3>0.5, na.rm=T), p4=sum(p_h4>0.5, na.rm=T)) %>%
    mutate(prop_p4 = round(p4/n,2), prop_p3 = round(p3/n,2) )
}

twas_predixcan_coloc_fraction_data <- function(metaxcan_twas, gencode) {
  HLA <- select_around_genes_2(gencode, chromosome="chr6", start=28866528, end=33775446)
  hla_genes <- unique(HLA$gene_name)
  
  # So, analyze anything we have from s-twas
  twas_available <- metaxcan_twas %>% 
    dplyr::filter(!is.na(twas_pvalue), !(gene_name %in% hla_genes))
  
  results <- data.frame()
  for (pheno in sort(unique(twas_available$phenotype))) {
    # We'll restrict S-PrediXcan to the same tissues available in twas, and predixcan_pvalue < largest twas (and good prediction performance)
    twas_ <- twas_available %>% 
      dplyr::filter(phenotype == pheno, !is.na(twas_pvalue)) 
    twas_tissues <- unique(twas_$tissue)
    max_twas_p <- max(twas_$twas_pvalue,na.rm=TRUE)
    
    result_twas_ <- twas_ %>% 
      get_coloc_proportions_() %>% mutate(phenotype=pheno, method='twas', genes="all")

    result_twas_no_hla_ <- twas_ %>% 
      dplyr::filter(!(gene_name %in% hla_genes)) %>% 
      get_coloc_proportions_() %>% mutate(phenotype=pheno, method='twas', genes="no_hla")
    
    #
    spredixcan_ <- metaxcan_twas %>% 
      dplyr::filter(phenotype == pheno, tissue %in% twas_tissues, pred_perf_pval<0.05, pval<max_twas_p)

    result_predixcan_ <- spredixcan_ %>% 
      get_coloc_proportions_() %>% mutate(phenotype=pheno, method='predixcan', genes="all")
    
    #
    result_predixcan_no_hla_ <- spredixcan_ %>%
      dplyr::filter(!(gene_name %in% hla_genes)) %>% 
      get_coloc_proportions_() %>% mutate(phenotype=pheno, method='predixcan', genes="no_hla")
    
    results <- rbind(results, result_twas_, result_twas_no_hla_, result_predixcan_, result_predixcan_no_hla_)
  }
  
  results
}

twas_predixcan_plot__ <- function(data) {
  data %>%
    ggplot2::ggplot(ggplot2::aes(x=predixcan, y=twas, label=phenotype)) +
    ggplot2::geom_point(size = 5, color = 'grey35') + 
    ggplot2::geom_abline(slope=1, intercept=0) +
    ggrepel::geom_label_repel(fontface = 'bold', size=6, segment.color = 'grey20', box.padding = grid::unit(0.5, "lines")) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust=0, size=28, face="bold"),
                   plot.subtitle= ggplot2::element_text(size=20, hjust=0.5, face="italic", color="black"),
                   axis.text = ggplot2::element_text(vjust = 0.5, hjust = 0.5, lineheight = 0,  size = 20,  face="bold"),
                   axis.title.x = ggplot2::element_text(size = 22, vjust=1,  face="bold"),
                   axis.title.y = ggplot2::element_text(size = 22, vjust=1,  face="bold"))
}

twas_predixcan_coloc_fraction_plot_ <- function(data) {
  twas_predixcan_plot__(data) + ggplot2::coord_cartesian(xlim=c(0.0, 1.0), ylim=c(0.0, 1.0))
}

twas_predixcan_n_plot <- function(data) {
  data %>%
    dplyr::select(phenotype, n, method) %>% tidyr::spread(key=method, value=n) %>% 
    twas_predixcan_plot__() +
    ggplot2::coord_cartesian(xlim=c(0.0, 1000), ylim=c(0.0, 1000)) +
    ggplot2::ggtitle("Number of significant genes:\n TWAS (BSLMM) vs PrediXcan (Elastic Net)")
}

twas_predixcan_p3_fraction_plot <- function(data) {
  data %>%
    dplyr::select(phenotype, prop_p3, method) %>% tidyr::spread(key=method, value=prop_p3) %>% 
    twas_predixcan_coloc_fraction_plot_() +
    ggplot2::ggtitle("Proportion of Non-Colocalization:\n TWAS (BSLMM) vs PrediXcan (Elastic Net)",subtitle="\n#independent/#total: independent:P3>0.5; HLA region excluded") 
}

twas_predixcan_p4_fraction_plot <- function(data) {
  data %>%
    dplyr::select(phenotype, prop_p4, method) %>% tidyr::spread(key=method, value=prop_p4) %>% 
    twas_predixcan_coloc_fraction_plot_() +
    ggplot2::ggtitle("Proportion of Colocalization:\n TWAS (BSLMM) vs PrediXcan (Elastic Net)",subtitle="\n#colocalized/#total: colocalized:P4>0.5; HLA region excluded") 
}