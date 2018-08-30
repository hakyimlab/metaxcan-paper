library(RPostgreSQL)
library(dplyr)

###############################################################################
drv <- RPostgreSQL::PostgreSQL()


db_v6p_hapmap_advanced_4_ <- function() {
  DBI::dbConnect(drv, host="development.ccsriudvs1y0.us-east-1.rds.amazonaws.com",
            port="5432",
            dbname="prototype_2017_07_03",
            user="metaxcan_ro",
            password="M3t4xc4n")
}

get_all_pheno_names <- function(connection) {
  p <- DBI::dbGetQuery(connection, "SELECT tag FROM pheno")
  return(p$tag)
}

get_pheno_info <- function(connection) {
  DBI::dbGetQuery(connection, "select * from pheno as p inner join pheno_info as pi on p.id = pi.pheno_id")
}

## dbListTables(db)
## dbListFields(db,'metaxcan_result')
## dbDisconnect(db)

###############################################################################
query_base <- paste0( "SELECT ",
                      " g.gene,",
                      " g.gene_name,",
                      " m.zscore,",
                      " m.effect_size,",
                      " m.pval,",
                      " p.tag as phenotype,",
                      " t.tissue as tissue,",
                      " m.pred_perf_R2,",
                      " m.pred_perf_pval,",
                      " m.n_snps_used,",
                      " m.n_snps_model",
                      " FROM gene AS g ",
                      " INNER JOIN metaxcan_result AS m ON g.id = m.gene_id ",
                      " INNER JOIN tissue AS t ON t.id = m.tissue_id ",
                      " INNER JOIN pheno AS p ON p.id = m.pheno_id " )

query_base_2 <- paste0(
  "SELECT",
  " g.gene_name,",
  " m.zscore,",
  " m.effect_size,",
  " m.pval,", 
  " p.tag as phenotype,",
  " t.tissue as tissue,",
  " m.pred_perf_R2,",
  " m.pred_perf_pval,",
  " m.pred_perf_qval,",
  " m.n_snps_used,",
  " m.n_snps_model,",
  " mi.p_smr,",
  " mi.p_heidi,",
  " mi.coloc_p0 as P_H0,",
  " mi.coloc_p1 as P_H1,",
  " mi.coloc_p2 as P_H2,",
  " mi.coloc_p3 as P_H3,",
  " mi.coloc_p4 as P_H4,",
  " tr.p as twas_pvalue,",
  " g.gene",
  " FROM gene AS g",
  " INNER JOIN metaxcan_result AS m ON g.id = m.gene_id ",
  " INNER JOIN tissue AS t ON t.id = m.tissue_id ",
  " INNER JOIN pheno AS p ON p.id = m.pheno_id ",
  " LEFT JOIN metaxcan_result_info as mi on m.id = mi.metaxcan_result_id ",
  " LEFT JOIN twas_relationship as trr ON m.id = trr.metaxcan_result_id ",
  " LEFT JOIN twas_result as tr ON trr.twas_result_id = tr.id "    
)

query_twas_with_metaxcan <- paste0(
  "SELECT",
  " g.gene,",
  " g.gene_name,",
  " p.tag as phenotype,",
  " t.tissue as tissue,",
  " tr.z as twas_zcore,",
  " tr.p as twas_pvalue",
  " FROM gene AS g",
  " INNER JOIN metaxcan_result AS m ON g.id = m.gene_id ",
  " INNER JOIN tissue AS t ON t.id = m.tissue_id ",
  " INNER JOIN pheno AS p ON p.id = m.pheno_id ",
  " INNER JOIN twas_relationship as trr ON m.id = trr.metaxcan_result_id ",
  " INNER JOIN twas_result as tr ON trr.twas_result_id = tr.id "    
)

get_twas_for_metaxcan <- function(connection) {
  DBI::dbGetQuery(connection, query_twas_with_metaxcan)
}

query.pheno.tissue <- function(connection, pheno, tissue)
{
  query <- paste0(query_base,
                  "WHERE p.tag = '", pheno,"'",
                  " and m.pred_perf_R2 > ", Rthres,
                  "and t.tissue = '" %&% tissue %&% "'")
  DBI::dbGetQuery(connection, query)
}

query_2_ <- function(connection, pheno=NULL, tissue=NULL, Rthres=0.) {
  query <- query_base_2
  where_clause <- ""
  if (!is.null(pheno)){
    where_clause <- paste0(" WHERE ", " p.tag = '" , pheno  ,"' ")
  }
  if (!is.null(tissue)){
    if (nchar(where_clause) == 0) { 
      where_clause <- " WHERE " 
    } else {
      where_clause <- paste0(where_clause, " AND ") 
    }
    where_clause <- paste0(where_clause, " t.tissue = '" , tissue  ,"' ")
  }
  if(Rthres>0) {
    if (nchar(where_clause) == 0) { 
      where_clause <- " WHERE " 
    } else {
      where_clause <- paste0(where_clause, " AND ") 
    }
    where_clause <- paste0(query, " and m.pred_perf_R2 >= " %&% Rthres)
  }
  query <- paste0(query, where_clause)
  DBI::dbGetQuery(connection, query)
}

build_data_2 <- function(connection, phenos, verbose=FALSE) {
  d <- data.frame()
  for(phenoname in sort(phenos)) {
    if (verbose) { print(paste0("Acquiring ", phenoname)) }
    data <- query_2_(connection, phenoname)
    data$phenotype <- as.character(data$phenotype)
    data$tissue < as.character(data$tissue)
    data$gene <- as.character(data$gene)
    data$gene_name <- as.character(data$gene_name)
    d <- rbind(d, data)
  }
  d <- d %>% arrange(phenotype, tissue)
  return(d)
}

build_data_3 <- function(connection, spec) {
  d <- data.frame()
  for(phenoname in sort(unique(spec$phenotype))) {
    for (tissuename in sort(unique(spec$tissue))) {
      data <- query_2_(connection, phenoname, tissuename)
      data$phenotype <- as.character(data$phenotype)
      data$tissue < as.character(data$tissue)
      data$gene <- as.character(data$gene)
      data$gene_name <- as.character(data$gene_name)
      d <- rbind(d, data)
    }
  }
  d <- d %>% arrange(phenotype, tissue)
  return(d)
}

###############################################################################
query.pheno <- function(connection, pheno,Rthres=NA)
{
  query <- paste0(query_base,
                  "WHERE p.tag = '", pheno, "'")
  if (!is.na(Rthres)) {
    query <- paste0(query," and m.pred_perf_R2 > ", Rthres)
  }
  DBI::dbGetQuery(connection, query)
}

query.all <- function(Rthres=0.)
{
  query <- paste0(query_base,
                  " and m.pred_perf_R2 > ", Rthres)
  DBI::dbGetQuery(db,query)
}

build_data <- function(connection, phenos) {
  d <- data.frame(pred_perf_r2=numeric(), 
                  pred_perf_pval = numeric(),
                  zscore = numeric(),
                  phenotype = character(),
                  tissue = character(),
                  gene = character(),
                  gene_name = character())
  for(phenoname in sort(phenos))
  {
    print(paste0("Acquiring ", phenoname))
    data <- query.pheno(phenoname, connection=connection)
    append <- data.frame(pred_perf_r2 = data$pred_perf_r2,
                         pred_perf_pval = data$pred_perf_pval,
                         zscore = data$zscore,
                         pval = data$pval,
                         phenotype = data$phenotype,
                         tissue = data$tissue,
                         gene = data$gene,
                         gene_name = data$gene_name,
                         stringsAsFactors=FALSE)
    d <- rbind(d, append)
  }
  d <- d %>% arrange(phenotype, tissue)
  return(d)
}

