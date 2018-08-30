

prepare_qq_data_ <- function(data, transform=NULL, threshold=30){
  d <- data %>% data.frame()
  if (!is.null(transform)) {
    d <- d %>% transform()
  }
  if(!("pvalue" %in% colnames(d))) {
    d$pvalue <-  2*pnorm(-abs(d$zscore))
  }
  
  nn <- nrow(d)
  d <- d %>% 
    dplyr::mutate(y=pmin(-log10(pvalue), threshold)) %>%
    dplyr::arrange(-y) %>%
    dplyr::mutate(x = -log10((1:nn)/(nn+1)))
  
  b <- -log10(0.05/nn) #bonferroni
  
  result <- d
  decoration <- data.frame(b=b,  stringsAsFactors = FALSE)
  
  list("result"=result, "decoration"=decoration)
}

prepare_qq_faceted_data_ <- function(data, facet_column=NULL, transform=NULL, threshold=30, facet_order=NULL) {
  facets <- unique(data[[facet_column]])
  result <- data.frame()
  decoration <- data.frame()
  for (facet in facets) {
    d <- data %>% dplyr::filter_(paste0(facet_column,"=='", facet, "'"))
    qq_ <- prepare_qq_data_(d, transform, threshold)
    
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


prepare_qq_ <- function(data, facet_column=NULL, transform=NULL, threshold=30, facet_order=NULL) {
  r <- if (!is.null(facet_column) && facet_column %in% names(data)) {
    prepare_qq_faceted_data_(data, facet_column, transform, threshold, facet_order)
  } else {
    prepare_qq_data_(data, transform, threshold)
  }
  r
}


qq_plot_config <- function(columns=2,
                           point_size=1,
                           legend_label="",
                           colour_labels=c("Data"="Data", "Decoration"="Decoration", "AdditionalData"="AdditionalData", "AdditionalDecoration"="AdditionalDecoration", "line"="line"),
                           colour_limits=c("Data", "AdditionalData"),
                           colours = c("Data"="black", "Decoration"="black", "AdditionalData"="dodgerblue3", "AdditionalDecoration"="dodgerblue3", "line"="black"),
                           scales="free",
                           facet_column="tissue",
                           facet_order=NULL,
                           threshold=30,
                           colour_column=NULL) {
  config <- list()
  config$point_size <- point_size
  config$scales <- scales
  config$columns <- columns
  config$facet_column <- facet_column
  config$facet_order <- facet_order
  config$threshold <- threshold
  
  config$legend_label <- legend_label
  config$colour_labels <- colour_labels
  config$colours <- colours
  config$colour_column <- colour_column
  config$colour_limits <- colour_limits
  
  config
}


build_qq_plot_ <- function(qq_data, config=qq_plot_config()) {
  data <- qq_data$result
  decoration <- qq_data$decoration
  
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::xlab(expression(Expected~~-log[10](italic(p)))) +
    ggplot2::ylab(expression(Observed~~-log[10](italic(p))))
  
  if ("Decoration" %in% config$colour_labels) {
    p <- p + ggplot2::geom_abline(data=decoration, mapping=ggplot2::aes(intercept=b, slope=0, colour='Decoration'), show.legend = F)
  } else {
    p <- p +
      ggplot2::geom_abline(data=decoration, mapping=ggplot2::aes(intercept=b, slope=0), colour="black", show.legend = F)
  }
  
  p <- p + ggplot2::geom_abline(data = data, mapping=ggplot2::aes(intercept=0, slope=1, colour='line'), show.legend = F)
  
  if (!is.null(config$colour_column)) {
    p <- p + ggplot2::geom_point(data = data, mapping=ggplot2::aes_string(x="x", y="y", colour=config$colour_column), size=config$point_size)    
  } else if ("Data" %in% names(config$colour_labels)) {
    p <- p + ggplot2::geom_point(data = data, mapping=ggplot2::aes(x=x, y=y, colour='Data'), size=config$point_size)
  } else {
    p <- p + ggplot2::geom_point(data = data, mapping=ggplot2::aes(x=x, y=y), colour="black",  size=config$point_size)
  }
  
  p <- p + ggplot2::scale_colour_manual(values = config$colours, labels = config$colour_labels, limits = config$colour_limits, name=config$legend_label)
  
  p
}


add_qq_plot_data_ <- function(p, qq_additional_data=NULL, config=qq_plot_config()) {
  additional_data <- qq_additional_data$result
  additional_decoration <- qq_additional_data$decoration
  
  if ("AdditionalDecoration" %in% config$colour_labels) {
    p <- p + ggplot2::geom_abline(data=additional_decoration, mapping=ggplot2::aes(intercept=b, slope=0, colour='AdditionalDecoration'), show.legend = F)
  } else {
    stop("Undecided colour for additional qq decoration")
  }
  
  if (!is.null(config$colour_column)) {
    p <- p + ggplot2::geom_point(data = additional_data, mapping=ggplot2::aes_string(x="x", y="y", colour=config$colour_column), size=config$point_size)
  } else if ("AdditionalData" %in% names(config$colour_labels)) {
    p <- p + ggplot2::geom_point(data = additional_data, mapping=ggplot2::aes(x=x, y=y, colour='AdditionalData'), size=config$point_size)
  } else {
    stop("Undecided colour for additional qq data")
  }
  
  p
}

simple_qq_plot_2 <- function(data, additional_data=NULL, transform=NULL, config=qq_plot_config()) {
  qq_data_ <- prepare_qq_(data, facet_column=config$facet_column,  transform=transform, threshold=config$threshold, facet_order=config$facet_order)
  do_facet <- "facet_w" %in% names(qq_data_$result)
  
  qq_p_ <- build_qq_plot_(qq_data_, config)
  
  if (!is.null(additional_data)) {
    qq_additional_ <- prepare_qq_(additional_data, facet_column=config$facet_column, transform = transform, threshold=config$threshold, facet_order=config$facet_order )
    qq_p_ <- add_qq_plot_data_(qq_p_, qq_additional_, config=config)
    if ("facet_w" %in% names(qq_additional_$result)) {
      do_facet <- TRUE
    }
  }
  
  if (do_facet) {
    qq_p_ <- qq_p_ + ggplot2::facet_wrap(~facet_w, scales = config$scales, ncol=config$columns)
  }
  
  qq_p_
}

simple_qq_plot <- function(data, additional_data=NULL, transform=NULL, facet_column=NULL, facet_order=NULL, columns=2, scales="fixed", threshold=30) {
  config <- qq_plot_config(columns=columns, scales=scales, threshold=threshold, facet_column=facet_column, facet_order=facet_order)
  simple_qq_plot_2(data, additional_data, transform, config)
}



