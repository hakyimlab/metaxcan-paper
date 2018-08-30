###############################################################################
build_heatmap_data <- function(d, pheno_selected) {
  result <- data.frame()
  for(phenoname in pheno_selected$pheno) {
    t <- d %>% dplyr::filter(phenotype == phenoname) 
    n <- nrow(t)
    t <- t %>% dplyr::group_by(tissue) %>% 
      summarise(meanZ2 = mean(zscore^2), n.signif.BF = sum(pval< 0.05/n), n.signif.e4 = sum(pval<1e-4),
                meanZ2.thres.1 = mean(zscore^2*(abs(zscore)>1) ), meanZ2.thres.2 = mean(zscore^2*(abs(zscore)>2) ), n.genes = n())
    if (length(t$meanZ2) == 0)
      next;
    t$phenotype <- phenoname
    result <- rbind(result, t)
  }
  result <-  result %>% dplyr::left_join(pheno_selected %>% 
                              select(pheno,pheno.short),by=c("phenotype"="pheno"))
  
  result
}

heatmap_plot <- function(data, tissue_colors) {
  colors <- tissue_colors %>% 
    dplyr::rename(tissue.short = tissue_site_detail_abbr) %>%
    dplyr::select(tissue_site_detail_id,tissue.short)
  t <- data %>%
    dplyr::left_join(colors, by=c("tissue"="tissue_site_detail_id"))

  
  t <- t %>% dplyr::filter(!(pheno.short %in% c('HEIGHT','RA')))
  p <- t %>% 
    ggplot2::ggplot(ggplot2::aes(pheno.short, tissue.short,fill=meanZ2)) +
    ggplot2::geom_tile() + 
    ggplot2::scale_fill_continuous(low="white", high="blue", name=expression(Mean~Z^2)) + 
    ggplot2::theme(
          plot.title = ggplot2::element_text(hjust=0.5, face="bold", colour="black", size=27),
          axis.text.y = ggplot2::element_text(face="bold", colour="black"),
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, lineheight = 10,hjust = 1,  face="bold", colour="black"),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank()) +
    ggplot2::ggtitle(expression(Mean~Z^2))
  p
}

plot_metaxcan_heatmap_data <- function(metaxcan, pheno_selected, tissue_colors, filename, h=800, w=800) {
  data <- build_heatmap_data(metaxcan, pheno_selected)
  p <- heatmap_plot(data, tissue_colors)
  
  png(filename, height=h, width=w)
  print(p)
  dev.off()
}