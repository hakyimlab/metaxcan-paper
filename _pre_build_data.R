suppressWarnings(source("_postgre_utilities.R"))
suppressWarnings(source("_helpers.R"))
suppressWarnings(library(readr))

data.dir <- 'data'
pheno_selected <- read_csv(file.path(data.dir,'selected-phenotypes.txt'))
replication_relationship <- read.delim("data/replicate_relationship.txt")
twas_phenos <- read.delim("data/twas_phenos.txt")

conn <- db_v6p_hapmap_advanced_4_()
metaxcan_results <- build_data_2(conn, pheno_selected$pheno)
metaxcan_results_replication <- build_data_2(conn, replication_relationship$pheno_repl)
pheno_info <- get_pheno_info(conn)
twas <- get_twas_for_metaxcan(conn)
DBI::dbDisconnect(conn)

ra_save_delim(metaxcan_results, "data/metaxcan_results_replication.txt")
ra_save_delim(metaxcan_results, "data/metaxcan_results.txt")
ra_save_delim(metaxcan_results_replication, "data/metaxcan_results_replication.txt")
ra_save_delim(pheno_info, "data/pheno_info.txt")
ra_save_delim(twas, "data/twas.txt")
           