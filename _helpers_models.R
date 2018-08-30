###############################################################################

load_model_extra_u <- function(file_path) {
  con <- DBI::dbConnect(RSQLite::SQLite(),file_path)
  
  query <- "SELECT gene, genename, R2 as `pred.perf.R2`, `n.snps` as `n.snps.in.model`, `pval` as `pred.perf.pval` FROM extra"
  results <- DBI::dbGetQuery(con, query)
  
  DBI::dbDisconnect(con)
  return(results)
}

load_models_extra_u <- function(folder, pattern=".db") {
  names <- list.files(folder)
  names <- names[grepl(pattern, names)]
  paths <- file.path(folder,names)
  dbs <- list()
  for(i in 1:length(paths)) {
    path <- paths[i]
    results <- load_model_extra_u(path)
    #match naming convention
    name <- gsub("_0.5.db", "", names[i])
    name <- gsub("TW_", "", name)
    name <- gsub("-", "_", name)
    dbs[[name]] <- results
  }
  return(dbs)
}

flatten_models <- function(models) {
  flat <- data.frame()
  for (name in names(models)) {
    model <- data.frame(models[[name]])
    model$tissue <- name
    flat <- rbind(flat, model)
  }
  flat
}

plot_model_extra_histograms <- function(models, column="pred.perf.pval", xlabel="Prediction Performance P-value", title="Prediction Models' Performance") {
  flat <- flatten_models(models)
  
  bins <- 30
  n <- nrow(flat)
  hy <- n/bins
  
  ggplot2::ggplot(data=flat, mapping=ggplot2::aes_string(x=column)) + 
    ggplot2::geom_histogram(colour="black", bins=bins) +
    ggplot2::geom_hline(yintercept=hy) +
    ggplot2::annotate("text", x = 0.7, y = hy+10e3, label = "(total count)/(number of bins)", size=7)+
    ggplot2::xlab(xlabel) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=20),
                   axis.title = ggplot2::element_text(size=15),
                   axis.text = ggplot2::element_text(size=12))
}