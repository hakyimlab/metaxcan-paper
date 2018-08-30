
###############################################################################
# It is very important that they all match chromosome
select_around_genes <- function(gencode, genes, window) {
  s <- gencode %>% dplyr::filter(gene_name %in% genes)
  pos_min <- min(s$start_location) - window
  pos_max <- max(s$end_location) + window
  chrom <- unique(s$chr)
  if (length(chrom) != 1) {
    stop("Bad list of genes")
  }
  
  gencode %>% dplyr::filter(chr == chrom, start_location >= pos_min, end_location <= pos_max)
}

select_around_genes_2 <- function(gencode, chromosome, start, end) {
  gencode %>% dplyr::filter(chr == chromosome, start_location >= start, end_location <= end)
}

###############################################################################
wrap_sentence <- function(x, ...) {
  paste(strwrap(x, ...), collapse = "\n")
}

#unusued
printable_string <- function(data, column) {
  k <- data[[column]]
  k <- gsub("_", " ", k)
  k <- gsub("\\.", " ", k)
  k
}

to_wrapped <- function(strings, length) {
  s <- stringi::stri_replace_all_fixed(strings, "_", " ")
  s <- stringi::stri_replace_all_fixed(s, ".", " ")
  s <- vapply(s, function(x){
    w <- stringi::stri_wrap(x, length, 0.0)
    return(paste(w, collapse="\n"))
  },FUN.VALUE="")
  s <- stringi::stri_replace_all_fixed(s, " ", "_")
  return(s)
}

printable_facet <- function(data, column, width=20) {
  s <- data[[column]]
  l <- sort(unique(s))
  k <- to_wrapped(l, length=width)
  d_ <- data.frame(x=l, y=k)
  d__ <- data.frame(x=s)
  d__ <- d__ %>% dplyr::inner_join(d_,by="x")
  data[[column]] <- factor(d__$y, k)
  data
}

###############################################################################
master_zscore_figure <- function(output, data, title, facet, col_x, col_y, ncols, w, h) {
  d <- data %>% printable_facet("phenotype") %>% printable_facet("tissue", width=16)
  facet_font_size <- if(facet == "tissue") { 12 } else { 14 }
  axis_text_font_size <- if (facet == "tissue") { 
    if (col_x == "rank(pred_perf_pval)"){ 4 } else { 8 }
  } else {12}
  
  pp <- 
    ggplot2::ggplot(d, ggplot2::aes_string(x=col_x, y=col_y)) + 
    ggplot2::geom_smooth() + 
    ggplot2::ggtitle(title) + 
    ggplot2::facet_wrap(facet,scales="free_y",ncol=ncols) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="bottom",
                   legend.text=ggplot2::element_text(size=25),
                   legend.title=ggplot2::element_text(size=25),
                   strip.text = ggplot2::element_text(size=facet_font_size, face="bold"),
                   strip.background = ggplot2::element_rect(fill="white"),
                   plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=27),
                   axis.title = ggplot2::element_text(size=25),
                   axis.text = ggplot2::element_text(size=axis_text_font_size)) + 
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

  png(output,width=w,height=h)
  print(pp)
  dev.off()
}

figure_zscore2_pred_perf_r2 <- function(output_name, data, facet_string, ncols, w, h) {
  facet <- as.formula(paste("~", facet_string))
  title <- expression(paste(zscore^{2}, " vs ", R^2))
  master_zscore_figure(output_name, data, title, facet, "pred_perf_r2", "zscore^2", ncols, w, h)
}

figure_zscore2_pred_perf_pval <- function(output_name, data, facet_string, ncols, w, h) {
  facet <- as.formula(paste("~", facet_string))
  title <- expression(paste(zscore^{2}, " vs ", p-value))
  master_zscore_figure(output_name, data,title, facet, "rank(pred_perf_pval)", "zscore^2", ncols, w, h)
}

###############################################################################

plot_density <- function(data, column, title=NULL, label=NULL) {
  ggplot2::ggplot(data=data) +
    ggplot2::geom_density(ggplot2::aes_string(x=column), size=1, adjust=2) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(label) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=30),
                   axis.title = ggplot2::element_text(size=25),
                   axis.text = ggplot2::element_text(size=20))
}

plot_histogram <- function(data, column, title=NULL, label=NULL, nbins=22, breaks=NULL) {

  p <- ggplot2::ggplot(data=data) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(label) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=37),
                   axis.title = ggplot2::element_text(size=30, face="bold"),
                   axis.text = ggplot2::element_text(size=25))
  if (!is.null(breaks)) {
    p + ggplot2::geom_histogram(ggplot2::aes_string(x=column), breaks=breaks, colour="black")
  } else {
    p + ggplot2::geom_histogram(ggplot2::aes_string(x=column), bins=nbins, colour="black")
  }
}


###############################################################################

c_theme_ <- function(plot, title_size=30, axis_title_size=25, axis_text_size=20) {
  scale <- 50.0
  plot <- plot + ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(limits = c(0,scale)) +
    ggplot2::scale_y_continuous(limits = c(0,scale)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust=0, size=title_size, face="bold"),
          plot.subtitle= ggplot2::element_text(size=25, hjust=0.5, face="italic", color="black"),
          axis.text = ggplot2::element_text(vjust = 0.5, hjust = 0.5, lineheight = 0,  size = axis_text_size,  face="bold"),
          axis.title = ggplot2::element_text(size = axis_title_size, vjust=1,  face="bold"),
          aspect.ratio = 1)
  return(plot)
}

smr_vs_metaxcan_plot <- function(d, subtitle=NULL){
  p <- ggplot2::ggplot(data = d , mapping = ggplot2::aes(x = log_p_metax, y = log_p_smr)) +
    ggplot2::geom_point(size = 2, color = "black") +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#404040") +
    ggplot2::labs(x = expression(bold('-log'[10]*'(p'[S-PrediXcan]*')')),
         y = expression(bold('-log'[10]*'(p'[SMR]*')')))
  
  p <- c_theme_(p, title_size = 37, axis_title_size = 30, axis_text_size = 27)
  p <- p + ggplot2::ggtitle("SMR vs. S-PrediXcan", subtitle=subtitle)
  return(p)
}

twas_vs_metaxcan_plot <- function(d, subtitle=NULL){
  p <- ggplot2::ggplot(data = d, mapping = ggplot2::aes(x = log_p_metax, y = log_p_twas)) + 
    ggplot2::geom_point(size = 2, color = "black") +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#404040") +
    ggplot2::labs(x = expression(bold('-log'[10]*'(p'[S-PrediXcan]*')')),
                  y = expression(bold('-log'[10]*'(p'[TWAS]*')')))
  
  p <- c_theme_(p)
  p <- p + ggplot2::ggtitle("S-TWAS vs. S-PrediXcan", subtitle = subtitle)
  return(p)
}

###############################################################################
#Depends  on "_postgre_utilities.R"

get_pheno_summary_ <- function(tag, connection) {
    pheno <- query_2_(connection, tag)
    b <- 0.05/nrow(pheno)
    pheno <- pheno %>% dplyr::mutate(p_h0_h1_h2=p_h0+p_h1+p_h2)
    
    significant <- pheno %>% dplyr::filter(pval < b)
    has_coloc <- sum(!is.na(significant$p_h4)) > 0
        
    b2 <- 0.05/nrow(significant)
    significant_pf <- significant %>% dplyr::filter(pred_perf_pval < b2)
    has_coloc_pf <- sum(!is.na(significant_pf$p_h4)) > 0
    
    data.frame(
      identifier = tag,
      significant_gene_tissue_pairs =  significant %>% nrow(),
      #significant_gene_tissue_pairs_pred_filtered = significant_pf %>% nrow(),
      #significant_gene_tissue_pairs_pred_filtered_exclude_p3 = if (has_coloc_pf) { significant_pf %>% dplyr::filter(p_h3<0.5) %>% nrow() } else { NA },
      #significant_gene_tissue_pairs_pred_filtered_keep_p4 = if (has_coloc_pf) { significant_pf %>% dplyr::filter(p_h4>0.5) %>% nrow() } else { NA },
      significant_unique_genes = significant %>% .$gene_name %>% unique() %>% length(),
      significant_unique_genes_pf = significant_pf %>% .$gene_name %>% unique() %>% length(),
      significant_unique_genes_pf_exclude_p3 = if (has_coloc_pf) { significant_pf %>% dplyr::filter(p_h3 < 0.5) %>% .$gene_name %>% unique() %>% length() } else { NA },
      significant_unique_genes_pf_keep_p4= if (has_coloc_pf) { significant_pf %>% dplyr::filter(p_h4 > 0.5) %>% .$gene_name %>% unique() %>% length() } else { NA }
    )
}

get_pheno_summary <- function(pheno_info, connection) {
  p_ <- pheno_info %>% data.frame() %>%
    dplyr::filter(hidden == FALSE) %>%
    dplyr::select(consortium=consortium,
                  study=name,
                  identifier=tag,
                  pubmed_paper_link=pubmed_paper_link,
                  sample_size=sample_size,
                  population=population) %>%
    dplyr::mutate(consortium=gsub("IM_GERA", "GERA", consortium))
  
  s_ <- data.frame()
  for (identifier in unique(p_$identifier)) {
    print(identifier)
    s_ <- rbind(s_, get_pheno_summary_(identifier, connection))
  }
  
  p_ %>% dplyr::inner_join(s_, by="identifier") %>% dplyr::arrange(consortium, study)
}

get_pheno_summary_2_ <- function(tag_, metaxcan_results) {
  pheno <- metaxcan_results %>% filter(phenotype == tag_)
  b <- 0.05/nrow(pheno)
  pheno <- pheno %>% dplyr::mutate(p_h0_h1_h2=p_h0+p_h1+p_h2)
  
  significant <- pheno %>% dplyr::filter(pval < b)
  has_coloc <- sum(!is.na(significant$p_h4)) > 0
  
  b2 <- 0.05/nrow(significant)
  significant_pf <- significant %>% dplyr::filter(pred_perf_pval < b2)
  has_coloc_pf <- sum(!is.na(significant_pf$p_h4)) > 0
  
  data.frame(
    identifier = tag_,
    significant_gene_tissue_pairs =  significant %>% nrow(),
    #significant_gene_tissue_pairs_pred_filtered = significant_pf %>% nrow(),
    #significant_gene_tissue_pairs_pred_filtered_exclude_p3 = if (has_coloc_pf) { significant_pf %>% dplyr::filter(p_h3<0.5) %>% nrow() } else { NA },
    #significant_gene_tissue_pairs_pred_filtered_keep_p4 = if (has_coloc_pf) { significant_pf %>% dplyr::filter(p_h4>0.5) %>% nrow() } else { NA },
    significant_unique_genes = significant %>% .$gene_name %>% unique() %>% length(),
    significant_unique_genes_pf = significant_pf %>% .$gene_name %>% unique() %>% length(),
    significant_unique_genes_pf_exclude_p3 = if (has_coloc_pf) { significant_pf %>% dplyr::filter(p_h3 < 0.5) %>% .$gene_name %>% unique() %>% length() } else { NA },
    significant_unique_genes_pf_keep_p4= if (has_coloc_pf) { significant_pf %>% dplyr::filter(p_h4 > 0.5) %>% .$gene_name %>% unique() %>% length() } else { NA }
  )
}

get_pheno_summary_2 <- function(pheno_info, metaxcan_results) {
  p_ <- pheno_info %>% data.frame() %>%
    dplyr::filter(hidden == FALSE) %>%
    dplyr::select(consortium=consortium,
                  study=name,
                  identifier=tag,
                  pubmed_paper_link=pubmed_paper_link,
                  sample_size=sample_size,
                  population=population) %>%
    dplyr::mutate(consortium=gsub("IM_GERA", "GERA", consortium))
  
  s_ <- data.frame()
  for (identifier in unique(p_$identifier)) {
    print(identifier)
    s_ <- rbind(s_, get_pheno_summary_2_(identifier, metaxcan_results))
  }
  
  p_ %>% dplyr::inner_join(s_, by="identifier") %>% dplyr::arrange(consortium, study)
}


ra_save_delim <- function(data, path){
  write.table(data, file=path, row.names = FALSE, sep="\t", quote=FALSE)
}

save_plot <- function(plot, path, height, width, res=NA) {
  png(path, height=height, width=width, res=res)
  print(plot)
  dev.off()
}

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
