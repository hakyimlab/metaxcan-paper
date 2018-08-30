
metaZn =  function(z_discovery, z_replication, sample_size_discovery, sample_size_replication)
{
  wdisc = sqrt(sample_size_discovery)
  wrepl = sqrt(sample_size_replication)
  ( z_discovery * wdisc  + z_replication * wrepl )/sqrt( wdisc^2 + wrepl^2 )
}

pi1f <- function(pvalue) {
  if(length(pvalue) == 0) return(NA)
  qobj <- qvalue::qvalue(p = pvalue)
  1 - qobj$pi0
}

metaxcan_replication_stats <- function(metaxcan_results, metaxcan_results_replication, replication_relationship, pheno_info) {
  results <- data.frame()
  for (pheno_discovery in replication_relationship$pheno_disc) {
    pheno_replication <- replication_relationship %>% dplyr::filter(pheno_disc == pheno_discovery) %>% .$pheno_repl
    
    #
    discovery <- metaxcan_results %>% dplyr::filter(phenotype == pheno_discovery, tissue!='DGN_WB')
    replication <- metaxcan_results_replication %>% dplyr::filter(phenotype == pheno_replication, tissue!='DGN_WB')
    pBF_thres <- 0.05/nrow(discovery)
    
    discovery_significant_genes <- discovery %>% dplyr::filter(pval< pBF_thres) %>%  dplyr::count(gene_name) %>% nrow()
    
    #
    m <- replication  %>% 
      dplyr::left_join(discovery, by=c("gene"="gene", "gene_name"="gene_name", "tissue"="tissue"), suffix=c("_r", "_d")) %>%
      dplyr::select(gene, gene_name, tissue, zscore_d, zscore_r, pval_d, pval_r, p_h3_d, p_h3_r, p_h4_d, p_h4_r, effect_size_d, effect_size_r)
      
    
    available_genes_replication <- m %>% dplyr::filter(pval_d < pBF_thres, !is.na(pval_r)) %>%  count(gene_name) %>% nrow()
    discovery_n <- pheno_info$sample_size[pheno_info$tag== pheno_discovery]
    replication_n <- pheno_info$sample_size[pheno_info$tag== pheno_replication]
    
    #
    m <- m %>% 
      dplyr::mutate(meta_pval = 2*pnorm(-abs(metaZn(zscore_d, zscore_r, discovery_n, replication_n))),
                    replication = pval_d<pBF_thres & pval_r<0.05 & sign(zscore_r)==sign(zscore_d) & meta_pval < pBF_thres)
    
    n_genes_repl = m %>% dplyr::filter(replication) %>% dplyr::count(gene_name) %>% nrow()
    
    r_ <- data.frame(phenotype_discovery = pheno_discovery,
                     phenotype_replication = pheno_replication,
                     discovery_significant_genes = discovery_significant_genes,
                     available_genes_replication = available_genes_replication,
                     pi1_all_inrepl = pi1f(m$pval_r),
                     pi1_dsig_inrep = pi1f(m$pval_r[m$pval_d < pBF_thres]),
                     prop_genes_repl = n_genes_repl / available_genes_replication,
                     n_genes_repl = n_genes_repl,
                     n_metasig = m %>% dplyr::filter(meta_pval < pBF_thres) %>% dplyr::count(gene_name) %>% nrow(),
                     n_genes_repl_exclude_p3 = m %>% dplyr::filter(replication, (p_h3_r < 0.5 | is.na(p_h3_r)), (p_h3_d < 0.5 | is.na(p_h3_d))) %>% dplyr::count(gene_name) %>% nrow(),
                     n_metasig_exclude_p3 = m %>% dplyr::filter(meta_pval < pBF_thres, (p_h3_r < 0.5 | is.na(p_h3_r)), (p_h3_d < 0.5 | is.na(p_h3_d))) %>% dplyr::count(gene_name) %>% nrow(),
                     n_genes_repl_keep_p4 = m %>% dplyr::filter(replication , p_h4_r > 0.5, p_h4_d > 0.5) %>% count(gene_name) %>% nrow(),
                     n_metasig_keep_p4 = m %>% dplyr::filter(meta_pval < pBF_thres, p_h4_r > 0.5, p_h4_d > 0.5) %>% count(gene_name) %>% nrow())
    results <- rbind(results, r_)
  }
  results
}

prop_annotation_d_ <- function(m, zthres) {
  results <- data.frame()
  l_ <- zthres-1
  for (pheno in unique(m$phenotype_discovery)) {
    prop <- m %>% dplyr::filter(phenotype_discovery == pheno) %>%
      dplyr::mutate(pos_dis=sign(Z_disc)>0,pos_repl=sign(Z_repl)>0)
    
    tag <- paste0(prop$pheno_disc_label[1], "\nReplication: ", prop$pheno_repl_label[1])
    prop <- prop  %>% dplyr::select(pos_repl,pos_dis) %>% table() %>% prop.table() %>% round(2) %>% data.frame()
    
    a <- c(prop[prop$pos_repl == FALSE & prop$pos_dis == FALSE,"Freq"],
           prop[prop$pos_repl == FALSE & prop$pos_dis == TRUE,"Freq"],
           prop[prop$pos_repl == TRUE & prop$pos_dis == FALSE,"Freq"],
           prop[prop$pos_repl == TRUE & prop$pos_dis == TRUE,"Freq"])
    
    r_ <- data.frame(x=c(-5, -5, 5, 5),
               y=c(-l_, l_, -l_, l_),
               t=paste0("Proportion: ",a),
               tag=tag)
    results <- rbind(results, r_)
  }

  results
}

plot_replication_ <- function(m, point_size=1) {
  zthres <- 10
  zmin <- 2
  m <- m %>% 
    dplyr::filter(abs(zscore_d)>zmin & abs(zscore_r)>zmin )  %>% 
    dplyr::mutate(Z_disc = pmin(abs(zscore_d), zthres) * sign(zscore_d), 
                  Z_repl=pmin(abs(zscore_r),zthres)*sign(zscore_r),
                  tag = paste0(pheno_disc_label, "\nReplication: ", pheno_repl_label))
  
  annotation <- prop_annotation_d_(m, zthres)
  
  m %>% 
    ggplot2::ggplot(mapping=ggplot2::aes(x=Z_disc, y=Z_repl)) + 
    ggplot2::geom_jitter(size=point_size) +
    ggplot2::geom_label(data=annotation, mapping=ggplot2::aes(x=x, y=y, label=t), fontface="bold", size=8, colour = "#335533") +
    ggplot2::geom_hline(color="gray", yintercept = 0) + ggplot2::geom_vline(color="gray", xintercept = 0) +
    ggplot2::xlim(-zthres,zthres) + ggplot2::ylim(-zthres,zthres) +
    ggplot2::theme_bw()
}

plot_replication <- function(metaxcan_results, metaxcan_results_replication, replication_relationship, point_size=1) {
  m <- metaxcan_results %>%
    dplyr::inner_join(replication_relationship, by=c("phenotype"="pheno_disc")) %>%
    dplyr::mutate(phenotype_discovery = phenotype, phenotype_replication=pheno_repl) %>%
    dplyr::inner_join(metaxcan_results_replication, by=c("phenotype_replication"="phenotype", "tissue"="tissue", "gene"="gene", "gene_name"="gene_name"), suffix=c("_d", "_r")) %>%
    dplyr::select(phenotype_discovery, phenotype_replication, zscore_d, zscore_r, pheno_disc_label, pheno_repl_label)
  plot_replication_(m, point_size=point_size)
}

plot_replication_simple <- function(metaxcan_results, metaxcan_results_replication, replication_relationship, point_size=1) {
  plot_replication(metaxcan_results, metaxcan_results_replication, replication_relationship) +
    ggplot2::theme(strip.text = ggplot2::element_text(face="bold", size=27),
                   plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=37),
                   plot.subtitle= ggplot2::element_text(size=30, hjust=0.5, face="italic"),
                   axis.title = ggplot2::element_text(face="bold", size=30),
                   axis.text = ggplot2::element_text(face="bold", size=25)) +
    ggplot2::xlab(paste0("Z-score: ", replication_relationship$pheno_disc_label)) + ggplot2::ylab(paste0("Z-score: ",replication_relationship$pheno_repl_label)) +
    ggplot2::ggtitle("S-PrediXcan Results Replication", subtitle = paste0(replication_relationship$pheno_disc_label, " and ", replication_relationship$pheno_repl_label))
}

plot_replication_faceted <- function(metaxcan_results, metaxcan_results_replication, replication_relationship, point_size=1) {
  ncols <- replication_relationship$pheno_disc %>% unique() %>% length() %>% sqrt() %>% ceiling()
  plot_replication(metaxcan_results, metaxcan_results_replication, replication_relationship) +
    ggplot2::facet_wrap(~tag, ncols) +
    ggplot2::theme(strip.text = ggplot2::element_text(face="bold", size=27),
                   strip.background = ggplot2::element_rect(fill="#DDDDDD"),
                   plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=37),
                   plot.subtitle= ggplot2::element_text(size=30, hjust=0.5, face="italic"),
                   axis.title = ggplot2::element_text(face="bold", size=30),
                   axis.text = ggplot2::element_text(face="bold", size=25)) +
    ggplot2::xlab("Z-score: Discovery") + ggplot2::ylab("Z-score: Replication") +
    ggplot2::ggtitle("S-PrediXcan Results Replication")
}

