
###############################################################################
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


###############################################################################

#metaxcan_results <- build_data_2(db_v6p_hapmap_advanced_4_(), pheno_selected$pheno)
#gencode <- read.delim("data/gencode.v19.txt.gz")
#gencode$gene_key <- gsub("\\..*","",gencode$gene) 

WINDOW <- 2e6
HLA <- select_around_genes(gencode, c("C4A", "C4B"), WINDOW)
hla_genes <- unique(HLA$gene_name)

metaxcan_without_hla <- metaxcan_results %>%  dplyr::filter(!(gene_name %in% hla_genes))