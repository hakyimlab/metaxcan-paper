###############################################################################
# Clinvar
find_clinvar_genes <- function(clinvar, disease_name) {
  c <- clinvar %>% dplyr::filter(grepl(disease_name, clinvar$DiseaseName, ignore.case = TRUE))
  unique(c$GeneSymbol) 
}

build_clinvar_filter <- function(clinvar, clinvar_pheno_selected) {
  filter <- data.frame()
  for (pheno_name in clinvar_pheno_selected$pheno) {
    disname <- (clinvar_pheno_selected %>% dplyr::filter(pheno == pheno_name))$search.term
    genes <- find_clinvar_genes(clinvar, disname)
    f <- data.frame(phenotype=pheno_name, gene_name=genes)
    filter <- rbind(filter, f)
  }
  filter
}

p4_filter_ <- function(data) {
  data %>% dplyr::filter(p_h4>0.5)
}

p3_filter_ <- function(data) {
  data %>% dplyr::filter(p_h3>0.5)
}

p_rest_filter_ <- function(data) {
  data %>% dplyr::filter(p_h0+p_h1+p_h2>0.5)
}

clinvar_pheno_selected_ <- function(pheno_selected, chosen_clinvar_tags) {
  pheno_selected %>% 
    dplyr::filter(pheno.short %in% chosen_clinvar_tags) %>%
    dplyr::mutate(pheno.short = factor(pheno.short, chosen_clinvar_tags)) %>%
    dplyr::arrange(pheno.short) %>%
    dplyr::mutate(pheno = factor(pheno, pheno), label = factor(label, label))
}

###############################################################################

clinvar_for_genes <- function(metaxcan_results, pheno_selected, clinvar,  chosen_clinvar_tags, genes)  {
  clinvar_pheno_selected <- clinvar_pheno_selected_(pheno_selected, chosen_clinvar_tags)
  
  m <- metaxcan_results %>%
    dplyr::inner_join(clinvar_pheno_selected %>% dplyr::select(phenotype=pheno, label=label), by="phenotype")
  
  clinvar_filter <- build_clinvar_filter(clinvar, clinvar_pheno_selected)
  m %>% 
    dplyr::inner_join(clinvar_filter, by=c("phenotype", "gene_name")) %>% 
    dplyr::filter(gene_name %in% genes) %>%
    dplyr::select(gene_name, phenotype, pval) %>% group_by(phenotype, gene_name) %>% summarize(n=n())

}
