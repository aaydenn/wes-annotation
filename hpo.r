#' @author Atakan Ayden
#' @title HPO enrichment analysis

home_path = "C://Users/atakan.ayden/"

setwd(home_path)

library(readr)
library(tidyverse)


#' OMIM ----------------------------

col_names_omim <- c("Chromosome","Genomic Position Start","Genomic Position End",
                    "Cyto Location","Computed Cyto Location","MIM Number",
                    "Gene Symbols","Gene Name","Approved Gene Symbol",
                    "Entrez Gene ID","Ensembl Gene ID","Comments","Phenotypes",
                    "Mouse Gene Symbol/ID")

omim <- read_delim("omim.txt", "\t", col_names = col_names_omim, comment = "#", col_types = c())


#' ORPHA ---------------------------

orpha <- read_delim("orpha.txt", delim = "\t")


#' HPO -----------------------------

col_names_hpo <- c("HPO_ID","HPO_label","Entrez_ID","Entrez_symbol",
                   "Additional Info from G-D source","G-D source","DatabaseID")

col_names_hpoa <- c("DatabaseID", "DiseaseName", "Qualifier", 
                    "HPO_ID", "Reference", "Evidence", "Onset", "Frequency", 
                    "Sex", "Modifier", "Aspect", "Biocuration")

hpo <- read_delim("phenotype_to_genes.txt", delim = "\t", 
                  col_names = col_names_hpo, comment = "#") %>% 
  select(c(1:4, 7))

hpoa <- read_delim("phenotype.hpoa", delim = "\t", 
                   comment = "#", col_names = col_names_hpoa) %>% 
  select(c(1,2,4))

hpo_db <-  hpo %>% left_join(hpoa, by = c("HPO_ID", "DatabaseID")) %>% separate(DatabaseID, c("Database", "ID"))

#' ----------------------------

# write_tsv(hpo_db, "hpo.tsv")
# save(hpo_db, file = "hpo.rda")

hpo2gene <- hpo_db %>% select(HPO_label, Entrez_symbol) %>% unique()


hpo2geneli <- hpo2gene %>% group_by(Entrez_symbol) %>% 
  group_map(~.x)

names(hpo2geneli) <- hpo2gene %>% group_by(Entrez_symbol) %>%
  group_map(~.y) %>%
  unlist()



hpo2disease <- hpo_db %>% select(HPO_label, DiseaseName) %>% unique()

hpo2diseaseli <- hpo2disease %>% group_by(DiseaseName) %>% 
  group_map(~.x)

names(hpo2diseaseli) <- hpo2disease %>% group_by(DiseaseName) %>%
  group_map(~.y) %>%
  unlist()

save(hpo_db, hpo2geneli, hpo2diseaseli, file = "hpo.rda")


#' ORA ----------------------------------------

my_endications <- sample(hpo_db$HPO_label,5) # fixed

# number_of_endications <- length(my_endications) # fixed
# total <- length(unique(hpo_db$HPO_label)) # fixed
# disease <- nrow(hpo2diseaseli$`12q14 microdeletion syndrome`)
# overlap <- length(intersect(my_endications, unlist(hpo2diseaseli$`12q14 microdeletion syndrome`)))

# ora <- phyper(q = overlap-1, m = disease, n = total-disease, k = number_of_endications, lower.tail = FALSE)
# ura <- phyper(q = overlap, m = disease, n = total-disease, k = number_of_endications, lower.tail = TRUE)

# over-representation analysis
ora_hpo_sub <- function(cond, endications, hpo, N){
  N <- length(unique(hpo))
  
  if (!is.na(cond)) {
    number_of_endications <- length(endications) # fixed
    
    if (length(hpo[[cond]])>0) {
      cond_m <- nrow(hpo[[cond]])
      overlap <- length(intersect(endications, unlist(hpo[cond])))
      
      ora <- phyper(q = overlap-1, 
                    m = cond_m, 
                    n = N-cond_m, 
                    k = number_of_endications, 
                    lower.tail = FALSE)
      
      ptable <- tibble(term = cond, p.val = ora, 
                       `e/g ratio` = paste0(overlap,"/",cond_m))
      return(ptable)
      
    } else { 
      cli::cat_bullet("no gene-related indication found for ", cond) 
    }
  }
  else {
    cli::cat_bullet("Invalid term!")
  }
}

# ora_hpo_sub("Acrodysostosis", my_endications, hpo = hpo2diseaseli)
# ora_hpo_sub("CDC4", my_endications, hpo = hpo2geneli)
# ora_hpo_sub("AARS1", my_endications, hpo = hpo2geneli)

ora_hpo <- function(endications, hpo) {
  
  if (hpo == "disease") {
    result <- pbapply::pblapply(names(hpo2diseaseli), ora_hpo_sub, endications = endications, hpo = hpo2diseaseli)
  } else if (hpo == "gene") {
    result <- pbapply::pblapply(names(hpo2geneli), ora_hpo_sub, endications = endications, hpo = hpo2geneli)
  }
  
  adj.tbl <- do.call(rbind, result)
  
  return(adj.tbl %>% mutate(adj.p.val = p.adjust(p.val, method = "fdr")) %>% arrange(adj.p.val))
}

ora_hpo(my_endications, hpo = "disease")

my_gene_enrichment <- ora_hpo(my_endications, hpo = "gene")









