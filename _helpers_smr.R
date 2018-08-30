load_smr <- function(folder, filter) {
  files <- list.files(folder)
  files <- files[grepl(filter, files, perl=TRUE)]
  results <- data.frame()
  for (file in files) {
    path <- file.path(folder, file)
    key <- stringr::str_match(file,filter)[,2]
    smr <- read.delim(path)
    smr$key <- key
    results <- rbind(results, smr)
  }
  return(results)
}

smr_breakdown_2 <- function(smr) {
  result <- rbind(data.frame(), smr)
  result <- result %>% select(gene=probeID, p_SMR, p_top_gwas=p_GWAS, p_top_eqtl=p_eQTL)
  return(result)
}

smr_vs_topeqtl_plot <- function(smr, subtitle=NULL) {
  ggplot2::ggplot(data=smr, mapping=ggplot2::aes(x=eqtl, y=smr)) +
    ggplot2::geom_point() + ggplot2::xlim(0, 50) + ggplot2::ylim(0, 50) +
    ggplot2::geom_abline(intercept=0, color='darkgray', size = 1) +
    ggplot2::labs(x=expression(bold('-log'[10]*'(p'['Top-eQTL']*')')), y=expression(bold('-log'[10]*'(p'['SMR']*')'))) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=37),
                   plot.subtitle= ggplot2::element_text(size=25, hjust=0.5, face="italic", color="black"),
                   axis.title = ggplot2::element_text(size=30, face="bold"),
                   axis.text = ggplot2::element_text(size=27, face="bold")) +
    ggplot2::ggtitle("SMR vs top eQTL", subtitle=subtitle)
}

smr_vs_topeqtlgwas_plot <- function(smr, subtitle=NULL) {
  ggplot2::ggplot(data=smr, mapping=ggplot2::aes(x=gwas, y=smr)) + 
    ggplot2::geom_point() + ggplot2::xlim(0, 50) + ggplot2::ylim(0, 50) +
    ggplot2::geom_abline(intercept=0, color='darkgray', size = 1) +
    ggplot2::labs(x=expression(bold('-log'[10]*'(p'['GWAS for top eQTL']*')')), y=expression(bold('-log'[10]*'(p'['SMR']*')'))) + 
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=37),
                   plot.subtitle= ggplot2::element_text(size=25, hjust=0.5, face="italic", color="black"),
                   axis.title = ggplot2::element_text(size=30, face="bold"),
                   axis.text = ggplot2::element_text(size=27, face="bold")) +
    ggplot2::ggtitle("SMR vs top eQTL's GWAS", subtitle=subtitle)
}