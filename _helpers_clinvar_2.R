source("_helpers_clinvar.R")

###############################################################################
prepare_qq_data_c <- function(data, threshold=30, transform_data_callback=NULL, transform_decoration_callback=NULL){
  d <- data.frame(data)
  if (!is.null(transform_data_callback)) {
    d <- transform_data_callback(d)    
  }
  
  d$y <- -log10(d$pvalue)
  d$y <- pmin(d$y, threshold) #upper threshold value
  d <- d %>% dplyr::arrange(-y)
  
  nn <- nrow(d)
  d$x <- -log10((1:nn)/(nn+1))
  
  result <- d %>% dplyr::select(x, y, label)
  
  b <- -log10(0.05/nn) #bonferroni
  decoration <- data.frame(b=b, stringsAsFactors = FALSE)
  if (!is.null(transform_decoration_callback)) {
    decoration <- transform_decoration_callback(decoration)
  }
  list("result"=result, "decoration"=decoration)
}

#' @export
prepare_qq_faceted_data_c <- function(data, facet_column=NULL, threshold=30, facet_order=NULL, transform_data_callback=NULL, transform_decoration_callback=NULL) {
  facets <- unique(data[[facet_column]])
  result <- data.frame()
  decoration <- data.frame()
  for (facet in facets) {
    d <- data %>% dplyr::filter_(paste0(facet_column,"=='", facet, "'"))
    qq_ <- prepare_qq_data_c(d, threshold=threshold, transform_data_callback=transform_data_callback, transform_decoration_callback=transform_decoration_callback)
    
    r <- qq_$result %>% dplyr::mutate(facet_w = facet)
    d <- qq_$decoration %>% dplyr::mutate(facet_w = facet)

    result <- rbind(result, r)
    decoration <- rbind(decoration, d)
  }
  
  if (!is.null(facet_order)) {
    result$facet_w <- factor(result$facet_w , facet_order)
    decoration$facet_w  <- factor(decoration$facet_w , facet_order)
  }
  
  list("result"=result, "decoration"=decoration)
}

#' @export
prepare_qq_faceted_data_2_c <- function(data, config, transform_data_callback=NULL, transform_decoration_callback) {
  prepare_qq_faceted_data_c(data, facet_column=config$facet_column, threshold = config$threshold, facet_order = config$facet_order, transform_data_callback=transform_data_callback, transform_decoration_callback=transform_decoration_callback)
}

###############################################################################
qq_plot_config_c <- function(columns=2,
                           point_size=1,
                           scales="free",
                           facet_column="tissue",
                           facet_order=NULL,
                           threhsold=30,
                           labels=c("Main", "Secondary"),
                           colours = c("black", "dodgerblue3")) {
  config <- list()
  config$point_size <- point_size
  config$scales <- scales
  config$columns <- columns
  config$facet_column <- facet_column
  config$facet_order <- facet_order
  config$threshold <- 30
  
  config$colours <- if(is.null(colours) || length(colours) < 2) {
    c("black", "dodgerblue3")
  } else {
    colours
  }
  
  config$labels <- if(is.null(labels) || length(labels) < 2) {
    c("Main", "Secondary")
  } else {
    labels
  }
  
  config
}

build_qq_plot_c <- function(qq_data, qq_additional_data=NULL, config=qq_plot_config_c()) {
  data <- qq_data$result
  decoration <- qq_data$decoration
  additional_data <- if(is.null(qq_additional_data)) { NULL } else {qq_additional_data$result}
  additional_decoration <- if(is.null(qq_additional_data)) { NULL } else {qq_additional_data$decoration}
  
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::xlab(expression(Expected~~-log[10](italic(p)))) +
    ggplot2::ylab(expression(Observed~~-log[10](italic(p))))
  
  if("facet_w" %in% names(data) ||
     (!is.null(additional_data) && "facet_w" %in% names(additional_data))) {
    p <- p + ggplot2::facet_wrap(~facet_w, scales=config$scales,ncol=config$columns)
  }
  
  p <- p +
    ggplot2::geom_abline(data=decoration, mapping=ggplot2::aes(intercept=b, slope=0, colour=label)) +
    ggplot2::geom_abline(data = data, mapping=ggplot2::aes(intercept=0, slope=1), colour='black') +
    ggplot2::geom_point(data = data, mapping=ggplot2::aes(x=x, y=y, colour=label), size=config$point_size)
  
  if (!is.null(additional_data)) {
    additional_map <- ggplot2::aes(x=x, y=y, colour=label)
    p <- p +
      ggplot2::geom_abline(data=additional_decoration, mapping=ggplot2::aes(intercept=b, slope=0, colour=label)) +
      ggplot2::geom_point(data=additional_data, mapping=additional_map, size=config$point_size)
    
    p <- p +
      ggplot2::theme(legend.position="bottom",legend.direction="horizontal")
    #p <- p + ggplot2::guides(color=TRUE)
  } else {
    p <- p + ggplot2::guides(color=FALSE)
  }

  p <- p +
    ggplot2::scale_colour_manual(name = 'Colours:',
                                values = config$colours, 
                                labels = config$labels)
  
  p
}

###############################################################################
clinvar_qq_plot <- function(metaxcan_results, pheno_selected, clinvar,  chosen_clinvar_tags, additional_filter=NULL, title_addition=NULL) {
  clinvar_pheno_selected <- clinvar_pheno_selected_(pheno_selected, chosen_clinvar_tags)
  
  m <- metaxcan_results %>%
    dplyr::inner_join(clinvar_pheno_selected %>% dplyr::select(phenotype=pheno, label=label), by="phenotype") %>%
    dplyr::rename("pvalue"="pval")
  
  
  clinvar_filter <- build_clinvar_filter(clinvar, clinvar_pheno_selected)
  mc <- m %>% 
    dplyr::inner_join(clinvar_filter, by=c("phenotype", "gene_name"))
  if (!is.null(additional_filter)) {
    mc <- additional_filter(mc)
  }
  
  labels <- c("All Results", "Clinvar")
  colours <- c("black", "dodgerblue3")
  qq_config <- qq_plot_config_c(columns=3,
                                point_size=2,
                                scales="free",
                                facet_column="label",
                                facet_order=levels(clinvar_pheno_selected$label),
                                threhsold=30,
                                colours=colours,
                                labels=labels)
  
  title <- "QQ plot:\n Summary PrediXcan results for all genes,\n and results for genes in Clinvar"
  if (!is.null(title_addition)) {
    title <- paste0(title, "\n", title_addition)
  }
  
  all_results_colour_ <- function(d) {
    d %>% dplyr::mutate(label = "All results")
  }
  qq_ <- prepare_qq_faceted_data_2_c(m, qq_config, transform_data_callback=all_results_colour_, transform_decoration_callback=all_results_colour_)
  
  p_colour_ <- function(d) {
    d %>% dplyr::mutate(label="Clinvar")
  }
  qqa_ <- prepare_qq_faceted_data_2_c(mc, qq_config, transform_data_callback = p_colour_, transform_decoration_callback=p_colour_)
  
  build_qq_plot_c(qq_, qqa_, qq_config) + 
    ggplot2::theme(legend.position="bottom",
                   legend.text=ggplot2::element_text(size=25),
                   legend.title=ggplot2::element_text(size=25),
                   strip.text = ggplot2::element_text(size=18, face="bold"),
                   strip.background = ggplot2::element_rect(fill="white"),
                   plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=27),
                   axis.title = ggplot2::element_text(size=25),
                   axis.text = ggplot2::element_text(size=20)) +
    ggplot2::guides(color=ggplot2::guide_legend("Colors :", override.aes = list(size=1, shape=c(19)))) + 
    ggplot2::ggtitle(title)
}

###############################################################################
clinvar_qq_plot_c <- function(metaxcan_results, pheno_selected, clinvar,  chosen_clinvar_tags, additional_filter=NULL, title_addition=NULL) {
  clinvar_pheno_selected <- clinvar_pheno_selected_(pheno_selected, chosen_clinvar_tags)
  
  m <- metaxcan_results %>%
    dplyr::inner_join(clinvar_pheno_selected %>% dplyr::select(phenotype=pheno, label=label), by="phenotype") %>%
    dplyr::rename("pvalue"="pval")
    
  
  clinvar_filter <- build_clinvar_filter(clinvar, clinvar_pheno_selected)
  mc <- m %>% 
    dplyr::inner_join(clinvar_filter, by=c("phenotype", "gene_name"))
  if (!is.null(additional_filter)) {
    mc <- additional_filter(mc)
  }
  
  labels <- c("All Results", "Clinvar", "Clinvar, p4>0.5", "Clinvar, p3>0.5", "Clinvar, Undetermined Colocalization")
  colours <- c("black", "dodgerblue3", "green", "yellow", "orange", "gray")
  qq_config <- qq_plot_config_c(columns=3,
                              point_size=2,
                              scales="free",
                              facet_column="label",
                              facet_order=levels(clinvar_pheno_selected$label),
                              threhsold=30,
                              colours=colours,
                              labels=labels)
  
  title <- "QQ plot:\n Summary PrediXcan results for all genes,\n and results for genes in Clinvar"
  if (!is.null(title_addition)) {
    title <- paste0(title, "\n", title_addition)
  }
  
  all_results_colour_ <- function(d) {
    d <- data.frame(d)
    d$label <- "All results"
    d
  }
  qq_ <- prepare_qq_faceted_data_2_c(m, qq_config, transform_data_callback=all_results_colour_, transform_decoration_callback=all_results_colour_)
  
  p_colour_ <- function(d) {
    d %>% dplyr::mutate(label = ifelse(is.na(p_h0), "Clinvar",
                                ifelse(d$p_h4>0.5, "Clinvar, p4>0.5",
                                ifelse(d$p_h3>0.5, "Clinvar, p3>0.5", "Clinvar, Undetermined Colocalization"))))
  }
  c_colour_ <- function(d) {
    d <- data.frame(d)
    d$label <- "Clinvar"
    d
  }
  qqa_ <- prepare_qq_faceted_data_2_c(mc, qq_config, transform_data_callback = p_colour_, transform_decoration_callback=c_colour_)
  
  build_qq_plot_c(qq_, qqa_, qq_config) + 
    ggplot2::theme(legend.position="bottom",
                   legend.text=ggplot2::element_text(size=25),
                   legend.title=ggplot2::element_text(size=25),
                   strip.text = ggplot2::element_text(size=18, face="bold"),
                   strip.background = ggplot2::element_rect(fill="white"),
                   plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=27),
                   axis.title = ggplot2::element_text(size=25),
                   axis.text = ggplot2::element_text(size=20)) +
    ggplot2::guides(color=ggplot2::guide_legend("Colors :", override.aes = list(size=1, shape=c(19)))) + 
    ggplot2::ggtitle(title)
}