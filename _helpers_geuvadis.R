
format_geuvadis_expression__ <- function(geuvadis_expression, individuals) {
  geuvadis_expression %>% 
    reshape2::melt(id.vars="TargetID", measure.vars=individuals, variable.name="id", value.name="expression") %>%
    tidyr::spread(TargetID, expression) %>% 
    dplyr::select(-id)
}

#Produces same result as function above but is faster
format_geuvadis_expression_ <- function(geuvadis_expression, individuals) {
  genes <- geuvadis_expression$TargetID %>% as.character()
  filtered <- geuvadis_expression %>% dplyr::select_(.dots=individuals)
  rownames(filtered) <- genes
  filtered <- t(filtered)
  filtered %>% data.frame()
}

get_cor_stats <- function(observed, model, gene, model_name, group) {
  obs <- scale(observed[[gene]])
  mod <- scale(model[[gene]])
  skip <- all(is.nan(obs)) || all(is.nan(mod))
  if (skip) {
    k <- list(estimate=NA, p.value=NA)
  } else {
    k <- cor.test(obs, mod) 
  }
  data.frame(model_name=model_name, gene=gene, group=group, cor=k$estimate, pvalue=k$p.value)
}

plot_ <- function(observed, model, gene) {
  obs <- scale(observed[[gene]])
  mod <- scale(model[[gene]])
  plot(obs, mod)
}

###############################################################################

process_geuvadis_expression <- function(geuvadis_expression, geuvadis_predicted_expression, geuvadis_samples) {

  #build filters that help speed up stuff
  samples <- geuvadis_samples %>% dplyr::filter(id %in% names(geuvadis_expression))
  individuals_ <- samples$id %>% as.character()
  african_ <- samples %>% dplyr::filter(group == "AFR") %>% .$id %>% as.character()
  european_ <- samples %>% dplyr::filter(group == "EUR") %>% .$id %>% as.character()
  africans <- individuals_ %in% african_
  europeans <- individuals_ %in% european_
  
  #Not sure about removing ensemble id
  observed <- geuvadis_expression %>% 
    dplyr::mutate(TargetID = ratools::remove_id_from_ensemble(TargetID)) %>%
    format_geuvadis_expression_(individuals_)
  
  observed_european <- observed %>% dplyr::filter(europeans)
  observed_african <- observed %>% dplyr::filter(africans)
  
  print("Processing geuvadis correlations")
  results <- data.frame()
  for (model_name in names(geuvadis_predicted_expression)) {
    print(model_name)
    model <- geuvadis_predicted_expression[[model_name]] %>% dplyr::filter(geuvadis_samples$id %in% individuals_)
    colnames(model) <- colnames(model) %>% ratools::remove_id_from_ensemble()
    available <- colnames(model)[colnames(model) %in% colnames(observed)]
    
    model_european <- model %>% dplyr::filter(europeans)
    model_african <- model %>% dplyr::filter(africans)
    
    for (gene in available) {
      results <- rbind(results, get_cor_stats(observed, model, gene, model_name, "ALL"))
      results <- rbind(results, get_cor_stats(observed_european, model_european, gene, model_name, "EUR"))
      results <- rbind(results, get_cor_stats(observed_african, model_african, gene, model_name, "AFR"))
    }
  }
  results
}

geuvadis_qq_plot <- function(geuvadis_correlation) {
  o_ <- geuvadis_correlation$model_name %>% unique() %>% sort() %>% as.character()
  all <- geuvadis_correlation %>% dplyr::filter(group == "ALL") %>% ratools::prepare_qq_(facet_column = "model_name", label="ALL", facet_order = o_)
  eur <- geuvadis_correlation %>% dplyr::filter(group == "EUR") %>% ratools::prepare_qq_(facet_column = "model_name", label="EUR", facet_order = o_)
  afr <- geuvadis_correlation %>% dplyr::filter(group == "AFR") %>% ratools::prepare_qq_(facet_column = "model_name", label="AFR",facet_order = o_)
 
  ncols <- ceiling(sqrt(length(o_)))
  
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::xlab(expression(Expected~~-log[10](italic(p)))) +
    ggplot2::ylab(expression(Observed~~-log[10](italic(p)))) +
    ggplot2::theme(legend.position="bottom",legend.direction="horizontal") +
    ggplot2::geom_abline(data=all$decoration, mapping=ggplot2::aes(intercept=b, slope=0), colour='black') +
    ggplot2::geom_abline(data=all$decoration, mapping=ggplot2::aes(intercept=0, slope=1), colour='black') +
    ggplot2::geom_point(data=all$result, mapping=ggplot2::aes(x=x, y=y, colour=label), size=2) +
    ggplot2::geom_point(data=eur$result, mapping=ggplot2::aes(x=x, y=y, colour=label), size=2) +
    ggplot2::geom_point(data=afr$result, mapping=ggplot2::aes(x=x, y=y, colour=label), size=2)
  
  if (ncols > 1) { 
    p <- p + ggplot2::facet_wrap(~facet_w, scales="fixed",ncol=ncols)
  }
  
  p <- p +
    ggplot2::scale_colour_manual(name = 'Colours:',
                                 values =c('ALL'="black",'EUR'="dodgerblue3", 'AFR'="red"))
  
  p 
}

###############################################################################

process_geuvadis_expression_2 <- function(geuvadis_expression, geuvadis_predicted_expression, geuvadis_samples) {
  #build filters that help speed up stuff
  samples <- geuvadis_samples %>% dplyr::filter(id %in% names(geuvadis_expression))
  individuals_ <- samples$id %>% as.character()
  african_ <- samples %>% dplyr::filter(group == "AFR") %>% .$id %>% as.character()
  european_ <- samples %>% dplyr::filter(group == "EUR") %>% .$id %>% sample(length(african_)) %>% as.character()
  africans <- individuals_ %in% african_
  europeans <- individuals_ %in% european_
  
  #Not sure about removing ensemble id
  observed <- geuvadis_expression %>% 
    dplyr::mutate(TargetID = ratools::remove_id_from_ensemble(TargetID)) %>%
    format_geuvadis_expression_(individuals_)
  
  observed_european <- observed %>% dplyr::filter(europeans)
  observed_african <- observed %>% dplyr::filter(africans)
  
  print("Processing geuvadis correlations")
  results <- data.frame()
  for (model_name in names(geuvadis_predicted_expression)) {
    print(model_name)
    model <- geuvadis_predicted_expression[[model_name]] %>% dplyr::filter(geuvadis_samples$id %in% individuals_)
    colnames(model) <- colnames(model) %>% ratools::remove_id_from_ensemble()
    available <- colnames(model)[colnames(model) %in% colnames(observed)]
    
    model_european <- model %>% dplyr::filter(europeans)
    model_african <- model %>% dplyr::filter(africans)
    
    for (gene in available) {
      results <- rbind(results, get_cor_stats(observed_european, model_european, gene, model_name, "EUR"))
      results <- rbind(results, get_cor_stats(observed_african, model_african, gene, model_name, "AFR"))
    }
  }
  list(results=results, selected_europeans=european_)
}

geuvadis_qq_plot_2 <- function(geuvadis_correlation) {
  o_ <- geuvadis_correlation$model_name %>% unique() %>% sort() %>% as.character()
  eur <- geuvadis_correlation %>% dplyr::filter(group == "EUR") %>% ratools::prepare_qq_(facet_column = "model_name", label="EUR", facet_order = o_)
  afr <- geuvadis_correlation %>% dplyr::filter(group == "AFR") %>% ratools::prepare_qq_(facet_column = "model_name", label="AFR",facet_order = o_)
  
  ncols <- ceiling(sqrt(length(o_)))
  
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::xlab(expression(Expected~~-log[10](italic(p)))) +
    ggplot2::ylab(expression(Observed~~-log[10](italic(p)))) +
    ggplot2::theme(legend.position="bottom",legend.direction="horizontal") +
    ggplot2::geom_abline(data=eur$decoration, mapping=ggplot2::aes(intercept=b, slope=0), colour='black') +
    ggplot2::geom_abline(data=eur$decoration, mapping=ggplot2::aes(intercept=0, slope=1), colour='black') +
    ggplot2::geom_point(data=eur$result, mapping=ggplot2::aes(x=x, y=y, colour=label), size=2) +
    ggplot2::geom_point(data=afr$result, mapping=ggplot2::aes(x=x, y=y, colour=label), size=2)
  
  if (ncols > 1) { 
    p <- p + ggplot2::facet_wrap(~facet_w, scales="fixed",ncol=ncols)
  }
  
  p <- p +
    ggplot2::scale_colour_manual(name = '',
                                 values =c('ALL'="black",'EUR'="dodgerblue3", 'AFR'="red"))
  
  p 
}
