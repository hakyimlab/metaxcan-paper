###############################################################################

is_listed <- function(text, list=NULL) {
  is <- FALSE;
  for (w in list) {
    if (grepl(w, text)) {
      is <- TRUE
      break;
    }
  }
  is
}

load_coloc <- function(path, pheno_selected, verbose=TRUE) {
  files <- list.files(path)
  filter_ <- vapply(files, is_listed,TRUE, list=pheno_selected$pheno)
  files <- files[filter_]
  
  results <- data.frame()
  for (pheno_ in pheno_selected$pheno) {
    filter_ <- vapply(files, grepl, TRUE, pattern=pheno_)
    files_ <- files[filter_]
    for (file_ in files_) {
      if (verbose) {print(file_)}
      name <- gsub(".txt", "", file_)
      comps <- strsplit(name, "_eQTL_")[[1]]
      
      p_ <- file.path(path, file_)
      d <- read.delim(p_)
      d$gene <- ratools::remove_id_from_ensemble(d$gene_id)
      d$phenotype <- gsub("_1KGf", "", comps[1]) #Remember that some phenotypes have allele frequencies from 1KG
      d$tissue <- comps[2]
      
      results <- rbind(results, d)
    }
  }
  results
}

overwrite_coloc_results <- function(metaxcan_results, coloc) {
  m <- metaxcan_results %>% 
    dplyr::left_join(coloc, by=c("phenotype", "tissue", "gene")) %>%
    dplyr::mutate(
      p_h0=ifelse(is.na(p_h0), P_H0, p_h0),
      p_h1=ifelse(is.na(p_h1), P_H1, p_h1),
      p_h2=ifelse(is.na(p_h2), P_H2, p_h2),
      p_h3=ifelse(is.na(p_h3), P_H3, p_h3),
      p_h4=ifelse(is.na(p_h4), P_H4, p_h4)
    ) %>%
    dplyr::select(-P_H0, -P_H1, -P_H2, -P_H3, -P_H4)
}

plot_coloc_tern <- function(d, title=NULL, vjust=-5, point_alpha=1.0, density_alpha=1.0) {
  p <- ggtern::ggtern(d, ggtern::aes(p_h3,p_h4, p_h0+p_h1+p_h2)) +
    ggplot2::scale_fill_gradient(low = "#EEEEFF", high = "#000020") +
    ggtern::stat_density_tern( mapping=ggtern::aes(fill=..level..), geom="polygon", base = "identity", alpha=density_alpha) +
    ggplot2::geom_point(shape = 16, size = 2, alpha=point_alpha, color='black') +
    ggtern::scale_T_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25)) +
    ggtern::scale_L_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25)) +
    ggtern::scale_R_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25)) +
    ggplot2::xlab("") + 
    ggplot2::ylab("") +
    ggtern::zlab("") +
    ggtern::Larrowlab("P3") +
    ggtern::Tarrowlab("P4") +
    ggtern::Rarrowlab("P0+P1+P2") +
    ggtern::theme_bw() +
    ggtern::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust=vjust, size=30, face="bold"),
                   tern.axis.arrow.text = ggplot2::element_text(hjust = 0.5, vjust=0.5, size=25, face="bold"),
                   tern.axis.text = ggplot2::element_text(hjust = 0.5, vjust=0.5, size=20, face="bold"),
                   panel.border = ggplot2::element_rect(fill = NA, colour = "black")) +
    ggtern::theme_showarrows() + 
    ggplot2::guides(fill = "none", alpha = "none")
  
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  
  p
}
