

#' @export
prepare_qq_data_ <- function(data, label=NULL, threshold=30){
  if("pvalue" %in% colnames(data)) {
    p_val <- data$pvalue
  } else {
    p_val <- 2*pnorm(-abs(data$zscore))
  }

  y <- -sort(log10(p_val))
  y <- pmin(y, threshold) #upper threshold value
  nn <- length(y)
  x <- -log10((1:nn)/(nn+1))
  b <- -log10(0.05/nn) #bonferroni


  result <- data.frame(x=x, y=y, stringsAsFactors = FALSE)
  decoration <- data.frame(b=b,  stringsAsFactors = FALSE)
  if (!is.null(label)){
    result$label <- label
    decoration$label <- label
  }

  list("result"=result, "decoration"=decoration)
}

#' @export
prepare_qq_faceted_data_ <- function(data, facet_column=NULL, label=NULL, threshold=30, facet_order=NULL) {
  facets <- unique(data[[facet_column]])
  result <- data.frame()
  decoration <- data.frame()
  for (facet in facets) {
    d <- data %>% dplyr::filter_(paste0(facet_column,"=='", facet, "'"))
    qq_ <- prepare_qq_data_(d, label = label, threshold=threshold)

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
prepare_qq_faceted_data_2_ <- function(data, config, label=NULL) {
  prepare_qq_faceted_data_(data, facet_column=config$facet_column, threshold = config$threshold, facet_order = config$facet_order, label=label)
}

#' @export
prepare_qq_ <- function(data, facet_column=NULL, label=NULL, threshold=30, facet_order=NULL) {
  r <- if (!is.null(facet_column) && facet_column %in% names(data)) {
    prepare_qq_faceted_data_(data, facet_column, label, threshold, facet_order)
  } else {
    prepare_qq_data_(data, label, threshold)
  }
  r
}

#' @export
qq_plot_config <- function(columns=2,
                           point_size=1,
                           labels=c("Main", "Secondary"),
                           colours = c("black", "dodgerblue3"),
                           scales="free",
                           facet_column="tissue",
                           facet_order=NULL,
                           threhsold=30) {
  config <- list()
  config$point_size <- point_size
  config$scales <- scales
  config$columns <- columns
  config$facet_column <- facet_column
  config$facet_order <- facet_order
  config$threshold <- 30

  config$labels <- if(is.null(labels) || length(labels) != 2) {
    c("Main", "Secondary")
  } else {
    labels
  }

  config$colours <- if(is.null(colours) || length(colours) != 2) {
    c("black", "dodgerblue3")
  } else {
    colours
  }

  config
}

#' @export
build_qq_plot_ <- function(qq_data, qq_additional_data=NULL, config=qq_plot_config()) {
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

  if (!is.null(additional_data)) {
    additional_map <- ggplot2::aes(x=x, y=y, colour="secondary")
    p <- p +
      ggplot2::geom_abline(data=additional_decoration, mapping=ggplot2::aes(intercept=b, slope=0, colour='secondary')) +
      ggplot2::geom_point(data=additional_data, mapping=additional_map, size=config$point_size)

    p <- p +
      ggplot2::theme(legend.position="bottom",legend.direction="horizontal")
    #p <- p + ggplot2::guides(color=TRUE)
  } else {
    p <- p + ggplot2::guides(color=FALSE)
  }

  p <- p +
    ggplot2::geom_abline(data=decoration, mapping=ggplot2::aes(intercept=b, slope=0, colour='primary')) +
    ggplot2::geom_abline(data = data, mapping=ggplot2::aes(intercept=0, slope=1), colour='black') +
    ggplot2::geom_point(data = data, mapping=ggplot2::aes(x=x, y=y, colour='primary'), size=config$point_size)

  p <- p +
    ggplot2::scale_colour_manual(name = 'Colours:',
        values =c('primary'=config$colours[1],'secondary'=config$colours[2]), labels = config$labels)

  p
}

#' @export
simple_qq_plot <- function(data, additional_data=NULL, facet_column=NULL, columns=2, scales="fixed", threshold=30) {
  qq_data_ <- prepare_qq_(data, facet_column=facet_column, threshold=threshold )

  qq_additional_data_ <- if (!is.null(additional_data)) {
    prepare_qq_(additional_data, facet_column=facet_column, threshold=threshold )
  } else {
    NULL
  }
  config <- qq_plot_config(columns=columns, scales=scales)
  qq_p_ <- build_qq_plot_(qq_data_, qq_additional_data_, config)
  return(qq_p_)
}

#' @export
simple_qq_plot_2 <- function(data, additional_data=NULL, config=qq_plot_config()) {
  qq_data_ <- prepare_qq_(data, facet_column=config$facet_column, threshold=config$threshold, facet_order=config$facet_order)

  qq_additional_data_ <- if (!is.null(additional_data)) {
    prepare_qq_(additional_data, facet_column=config$facet_column, threshold=config$threshold, facet_order=config$facet_order)
  } else {
    NULL
  }
  qq_p_ <- build_qq_plot_(qq_data_, qq_additional_data_, config)
  return(qq_p_)
}

