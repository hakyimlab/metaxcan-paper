
do_ecaviar_plots <- function(caviar_files, caviar_folder, twas, metaxcan_twas, gencode, output_folder) {
  HLA <- select_around_genes_2(gencode, chromosome="chr6", start=28866528, end=33775446)
  hla_genes <- unique(HLA$gene_name)
  
  phenos <- twas$phenotype %>% unique()
  for (i in 1:length(phenos)) {
    pheno_ <- phenos[i]
    t_ <- twas %>% dplyr::filter(phenotype == pheno_)
    m_ <- metaxcan_twas %>% dplyr::filter(phenotype == pheno_)
    
    caviar_file <- caviar_files %>% dplyr::filter(related_gwas == pheno_)
    ecaviar_path <- file.path(caviar_folder, caviar_file$file)
    ecaviar <- read.delim(ecaviar_path, sep= " ") %>% 
      rename(tissue=Tissue, rsid=rsID, gene=geneID, clpp=CLPP) %>%
      mutate(gene= gsub("\\.(.*)", "", gene))
    
    
    output_path <- file.path(output_folder, paste0("ecaviar_",pheno_,".png"))
    ecaviar_plot_(ecaviar, m_,  hla_genes) %>% ratools::save_plot(output_path, height=600, width=600)
  }
}

ecaviar_plot_ <- function(ecaviar, metaxcan_twas_p, exclude_genes=NULL, bare_bones=FALSE) {
  e_ <- ecaviar %>% dplyr::group_by(gene, tissue) %>% dplyr::top_n(1,clpp)

  twas_available <- metaxcan_twas_p %>% dplyr::filter(!is.na(twas_pvalue))
  if (!is.null(exclude_genes)){ twas_available <- twas_available %>%  dplyr::filter(!(gene_name %in% exclude_genes)) }

  twas_tissues <- unique(twas_available$tissue)
  max_twas_p <- max(twas_available$twas_pvalue,na.rm=TRUE)

  #
  spredixcan_ <- metaxcan_twas_p %>% 
    dplyr::filter(phenotype == pheno, tissue %in% twas_tissues, pred_perf_pval<0.05, pval<max_twas_p)

  twas_ <- twas_available %>% dplyr::inner_join(e_, by=c("gene", "tissue")) %>% dplyr::select(gene, tissue, phenotype, clpp) %>% dplyr::mutate(method="TWAS")
  spredixcan_ <- spredixcan_ %>% dplyr::inner_join(e_, by=c("gene", "tissue")) %>% dplyr::select(gene, tissue, phenotype, clpp) %>% dplyr::mutate(method="S-Predixcan")

  d_ <- rbind(twas_, spredixcan_)

  p <- ggplot2::ggplot(d_)
  p <- if (bare_bones) {
    p + ggplot2::geom_boxplot(ggplot2::aes(x=phenotype, y=clpp, fill=method), outlier.shape=NA)
  } else {
    p + ggplot2::geom_boxplot(ggplot2::aes(x=phenotype, y=clpp, fill=method))
  }
  p
}
