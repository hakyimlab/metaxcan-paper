###############################################################################

load_predixcan_data <- function(template, prefix, tokens) {
  results <- list()
  for (t in tokens) {
    path <- sprintf(template, prefix, t)
    p <- read.delim(path, sep = " ")
    p$gene_key <- gsub("\\..*","", p$gene) 
    results[[t]] <- p
  }
  
  results
}

m_key <- function(b, c) {sprintf("b_%s_c_%s",b, c)}

load_metaxcan_data <- function(template,prefix, tokens) {
  results <- list()
  for (t1 in tokens) {
    for (t2 in tokens) {
      k <- m_key(t1, t2)
      p <- sprintf(template, prefix, t1, t2)
      results[[k]] <- read.csv(p, stringsAsFactors = FALSE)
    }
  }
  
  results
}

###############################################################################

lm_eqn <- function(df){
  m <- lm(predixcan_z ~ metaxcan_z, df)
  r2_ <- format(summary(m)$r.squared, digits = 3)
  eq <- if (r2_ == "1") { 
    substitute(italic(r)^2 %~~% r2, list(r2 = "1")) 
  } else { 
    substitute(italic(r)^2~"="~r2, list(r2 = r2_)) 
  } 
  as.character(as.expression(eq));
}

m_p_plot_config <- function(r_font_size=17, axis_title_font_size=40, axis_text_size=20, facet_font_size=35, columns=1, r_x=-4, r_y=3.8) {
  config <- list()
  config$r_font_size <- r_font_size
  config$axis_title_font_size <- axis_title_font_size
  config$axis_text_size <- axis_text_size
  config$facet_font_size <- facet_font_size
  config$r_x <- r_x
  config$r_y <- r_y
  config$columns <- columns
  config
}

m_p_plot <- function(data, r_support=NULL, config=m_p_plot_config()) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x=metaxcan_z, y=predixcan_z)) + 
    ggplot2::theme_bw()
  
  p <- p +  ggplot2::labs(y = "PrediXcan Association Result", x = "S-PrediXcan Association Result") +
    ggplot2::geom_abline(intercept=0, slope=1, colour="grey49") +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size=config$axis_text_size, face="bold"),
      axis.title.y = ggplot2::element_text(size = config$axis_title_font_size, face="bold", margin=ggplot2::margin(0,20,0,10)),
      axis.title.x = ggplot2::element_text(size = config$axis_title_font_size, face="bold", margin=ggplot2::margin(20,0,10,0)),
      strip.text = ggplot2::element_text(size=config$facet_font_size, face ="bold"))
  
  if (!is.null(r_support)) {
    p <- p + 
      ggplot2::geom_text(data=r_support, 
                   ggplot2::aes(x = x, y = y, label = R2), 
                  fontface="bold", size=config$r_font_size, hjust=0.0, parse = TRUE, colour = "#335533")
    
  }
  p
}

m_p_facet_plot <- function(data, config=m_p_plot_config()) {
  support <- data %>% dplyr::distinct(the_facet, R2)
  support$x <- config$r_x
  support$y <- config$r_y

  p <- m_p_plot(data, support, config) +
    ggplot2::facet_wrap(~the_facet,scales="fixed", ncol=config$columns) +
    ggplot2::theme(panel.spacing  = ggplot2::unit(2, "lines"))
  p
}

m_p_grid_plot <- function(data, config=m_p_plot_config()) {
  support <- data %>% dplyr::distinct(study,reference, R2)
  support$x <- config$r_x
  support$y <- config$r_y

  p <- m_p_plot(data, support, config) +
    ggplot2::facet_grid(reference ~ study) 
  p
}

merge_p_m <- function(p, m, by_p="gene_key", by_m="gene"){
  by_c <- c()
  by_c[[by_p]] <- by_m
  r <- p %>% dplyr::inner_join(m, by = by_c)
  r <- r %>% dplyr::select(gene, metaxcan_z=zscore, predixcan_z=z.stat)
  r$R2 <- lm_eqn(r)
  r
}

build_predixcan_metaxcan_plot_data <- function(predixcan_results, metaxcan_results, tokens) {
  data <- data.frame()
  for (population in tokens) {
    for (target in tokens) {
      p <- predixcan_results[[target]]
      #Careful of this. We want to compare Predixcan(population=target) with
      #Metaxcan(GWAS=target, covariance=reference)
      m <- metaxcan_results[[m_key(b=target, c=population)]] 
      r <- p %>% dplyr::inner_join(m, by = c("gene_key" = "gene"))
      r <- r %>% dplyr::select(gene, metaxcan_z=zscore, predixcan_z=z.stat)
      r$R2 <- lm_eqn(r)
      r$study <- sprintf("%s GWAS Study", toupper(target))
      r$reference <- sprintf("%s Reference", toupper(population))
      data <- rbind(data,r)
    }
  }
  data
}


