

find_genes_more_significant_than_gwas <- function(gwas, gencode, data, pheno_tag, window_semi=1e6) {
  m <- data %>% 
    dplyr::filter(phenotype == pheno_tag, tissue != "DGN_WB") %>% 
    dplyr::select(phenotype, gene, gene_name, tissue, pval)
  
  m_p <- m %>% dplyr::inner_join(gencode, by = c("gene" = "gene_key", "gene_name"="gene_name"))
  #explicit logic as dplyr doesnt handel inequality join
  # and sqldf is too slow
  genes <- unique(m_p$gene)
  n <- length(genes)
  
  bonferroni <- 0.05/nrow(m_p)
  m_p <- m_p %>% dplyr::filter(pval < bonferroni) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::top_n(n = 1, wt = -pval) %>%  
    dplyr::mutate(window_start = ifelse(start_location-window_semi<0, 0, start_location-window_semi),
           window_end = end_location+window_semi)
  top_snp <- data.frame()
  for (i in 1:length(genes)) {
    the_gene <- genes[i]
    row <- m_p %>% dplyr::filter(gene == the_gene)
    g_chr <- row$chr[1]
    w_s <- row$window_start[1]
    w_e <- row$window_end[1]
    snps <- gwas %>% dplyr::filter(chrx==g_chr, w_s<=start_position, start_position<=w_e)
    n_ <- nrow(snps)
    snps <- snps[which.min(snps$pvalue),]
    if (length(snps$rsid) <1) next
    if (snps$pvalue[1] < row$pval[1]) next
    top_snp <- rbind(top_snp, data.frame(gene=the_gene, top_rsid=snps$rsid[1], top_snp_pvalue=snps$pvalue[1]))
    print(paste0(i," ", m$gene_name[i], " n_snps:", n_))
  }
  m_p <- m_p %>% 
    dplyr::inner_join(top_snp, by="gene") %>%
    dplyr::select(phenotype, gene, tissue, gene_name, pval, top_rsid, top_snp_pvalue)
  m_p
}